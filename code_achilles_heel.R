## ----setup, include=FALSE--------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_packages, results="hide", message=FALSE, warning=FALSE----

# LOAD PACKAGES ---------------------------------------------------------------

# Function to read in all required packages in one go:
loadPackages <- function(x) {
  for(i in x) {
    if(!require(i, character.only = TRUE)) {
      install.packages(i, dependencies = TRUE)
      library(i, character.only = TRUE)
    }
  }
}

loadPackages(c("data.table", "sensobol", "ncdf4", 
               "rworldmap", "sp", "countrycode", 
               "IDPmisc", "boot", "parallel", "scales",
               "MASS", "doParallel", "complmrob", 
               "mvoutlier", "sandwich", "lmtest", "mice", 
               "ggridges", "broom", "naniar", "cowplot", 
               "benchmarkme", "tidyverse"))

# SET CHECKPOINT --------------------------------------------------------------

dir.create(".checkpoint")

library("checkpoint")

checkpoint("2020-05-17", 
           R.version ="3.6.3", 
           checkpointLocation = getwd())

# CUSTOM FUNCTION TO DEFINE THE PLOT THEMES -----------------------------------

theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
}


## ----functions_data, cache=TRUE--------------------------------

# CREATE FUNCTIONS -----------------------------------------------------------------

# Function to obtain UN code, Continent and Country names
country_code <- function(dt) {
  dt[, `:=` (Code = countrycode(dt[, Country], 
                                origin = "country.name", 
                                destination = "un"), 
             Continent = countrycode(dt[, Country], 
                                     origin = "country.name", 
                                     destination = "continent"))]
  dt[, Country:= countrycode(dt[, Code], 
                             origin = "un", 
                             destination = "country.name")]
  setcolorder(dt, c("Country", "Continent", "Code", "Water.Withdrawn"))
  return(dt)
}

## Function to transform longitude and latitude to country.
# It is borrowed from Andy:
# https://stackoverflow.com/questions/14334970/convert-latitude-and-longitude-coordinates-to-country-name-in-r)
coords2country = function(points) {  
  countriesSP <- rworldmap::getMap(resolution = 'low')
  pointsSP = sp::SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  indices = sp::over(pointsSP, countriesSP)
  indices$ADMIN  
  #indices$ISO3 # returns the ISO3 code 
  #indices$continent   # returns the continent (6 continent model)
  #indices$REGION   # returns the continent (7 continent model)
}

# Function to load and extract data from .nc files produced
# by Huang et al. 2018
get_nc_data <- function(nc_file) {
  nc <- nc_open(nc_file)
  ww <- ncvar_get(nc, "withd_irr")
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  water <- rowSums(ww[, 469:ncol(ww)]) # Obtain year values for 2010 only
  ww.df <- data.frame(cbind(lon, lat, water)) 
  countries <- coords2country(ww.df[1:nrow(ww.df), 1:2])
  df <- cbind(countries, ww.df)
  setDT(df)
  final <- df[, .(Water.Withdrawn = sum(water)), countries]
  setnames(final, "countries", "Country")
  country_code(final)
  out <- na.omit(final[order(Continent)])
  out[, Water.Withdrawn:= Water.Withdrawn / 1000] # From mm to m
  return(out)
}

# Function to load and extract data from .nc files from ISIMIP
open_nc_files <- function(file, dname) {
  ncin <- nc_open(file)
  # get longitude, latitude, time
  lon <- ncvar_get(ncin, "lon")
  lat <- ncvar_get(ncin, "lat")
  # Get variable
  tmp_array <- ncvar_get(ncin, dname)
  # Retrieve last 12 months
  # get last year
  m <- (dim(tmp_array)[3] - 11):dim(tmp_array)[3]
  tmp_slice <- lapply(m, function(m) tmp_array[, , m])
  # create dataframe -- reshape data
  # matrix (nlon*nlat rows by 2 cols) of lons and lats
  lonlat <- as.matrix(expand.grid(lon,lat))
  # vector of `tmp` values
  tmp_vec <- lapply(tmp_slice, function(x) as.vector(x))
  # create dataframe and add names
  tmp_df01 <- lapply(tmp_vec, function(x) data.frame(cbind(lonlat, x)))
  names(tmp_df01) <- m
  da <- lapply(tmp_df01, data.table) %>%
    rbindlist(., idcol = "month") %>%
    na.omit()
  # Convert coordinates to country
  countries <- coords2country(da[1:nrow(da), 2:3])
  df <- cbind(countries, da)
  setDT(df)
  out <- na.omit(df)[, .(Water.Withdrawn = sum(x)), countries]
  out[, Water.Withdrawn:= Water.Withdrawn * 10000]
  return(out)
}


## ----water_with_dataset, cache=TRUE, warning=FALSE-------------

# READ IN DATASETS ON IRRIGATION WATER WITHDRAWAL ----------------------------------

# Define vector to reorder columns
cols_order <- c("Continent", "Country", "Code", "Water.Dataset", "Water.Withdrawn")

# FAO data (Table 4) ----------------------------
# UNIT IS KM3/YEAR
table4.tmp <- fread("table_4.csv", skip = 3, nrows = 167) %>%
  .[, .(Country, Year, Water.withdrawal)] %>%
  setnames(., "Water.withdrawal", "Water.Withdrawn") 

# Extract the selected years
table4.dt <- country_code(table4.tmp[Year %in% 1999:2012])[
  , Water.Dataset:= "Aquastat"][
    , Year:= NULL] %>%
  setcolorder(., cols_order)

# Liu et al. dataset ----------------------------
#UNIT IS 10^9 m3/year = km3/year
liu.dt <- fread("liu.csv")[, .(country, irr)] %>%
  setnames(., c("country", "irr"), c("Country", "Water.Withdrawn")) %>%
  country_code(.) %>%
  .[, Water.Dataset:= "Liu et al. 2016"] %>%
  setcolorder(., cols_order)
  

# Huang et al datasets --------------------------
names_nc_files <- c("withd_irr_lpjml.nc", "withd_irr_pcrglobwb.nc", 
                    "withd_irr_h08.nc", "withd_irr_watergap.nc")
out.nc <- mclapply(names_nc_files, function(x) get_nc_data(x), mc.cores = detectCores() * 0.75)
names(out.nc) <- c("LPJmL", "PCR-GLOBWB", "H08", "WaterGap")
out.final <- rbindlist(out.nc, idcol = "Water.Dataset")

# ISIMIP datasets --------------------------------
files <- list("dbh_wfdei_nobc_hist_varsoc_co2_airrww_global_monthly_1971_2010.nc", 
              "mpi-hm_miroc5_ewembi_picontrol_histsoc_co2_airrww_global_monthly_1861_2005.nc")

dname <- "airrww"

isimip.dt <- mclapply(files, function(x) open_nc_files(file = x, dname = dname), 
                      mc.cores = detectCores() * 0.75)

names(isimip.dt) <- c("DBHM", "MPI-HM")

isimip.final <- rbindlist(isimip.dt, idcol = "Water.Dataset") %>%
  setnames(., "countries","Country") %>%
  country_code(.) %>%
  na.omit()

GHM.dt <- rbind(out.final, isimip.final) %>%
  .[order(Country)]

# The value for Gabon by MPI-HM is a clear outlier
GHM.dt[Country == "Gabon"]

# So what we do is to give it a NA
GHM.dt <- GHM.dt[, Water.Withdrawn:= ifelse(Country == "Gabon" & 
                                              Water.Dataset == "MPI-HM", NA, 
                                            Water.Withdrawn)] %>%
  setcolorder(., cols_order)


## ----final_water_dataset, cache=TRUE, dependson=c("water_with_dataset" ,"arrange_total_countries")----

# CREATE THE FINAL IRRIGATION WATER WITHDRAWAL DATASET -----------------------------

water.dt <- rbind(liu.dt, table4.dt, GHM.dt)

# Check if there are duplicated countries
duplicated.countries <- water.dt[, .N, .(Country)][N > 8][, Country]

# Show duplicated countries
water.dt[Country %in% duplicated.countries]

# Compute mean
water.dt <- water.dt[, .(Water.Withdrawn = mean(Water.Withdrawn)),
                     .(Water.Dataset, Country, Code, Continent)]


## ----area_dataset, cache=TRUE----------------------------------

# READ IN IRRIGATED AREA DATASET ---------------------------------------------------

meier.dt <- fread("meier.csv") %>%
  setnames(., "Codes", "Code") %>%
  na.omit() %>%
  .[, .(Country, Continent, Code, `FAO-GMIA`)] %>%
  setnames(., "FAO-GMIA", "Irrigated.Area")


## ----merge_with_area, cache=TRUE, dependson=c("area_dataset", "water_with_dataset", "final_water_dataset")----

# MERGE DATASETS -------------------------------------------------------------------

full.dt <- water.dt[, merge(.SD, meier.dt, by = c("Country", "Code", "Continent"), 
                            all.y = TRUE), Water.Dataset] %>%
  .[!Continent == "Oceania"]

# Show countries with missing values in Water Withdrawal
full.dt[is.na(Water.Withdrawn), ] %>%
  .[, unique(Country)]


## ----plot_merged, cache=TRUE, dependson="merge_with_area", dev="tikz", fig.height=3, fig.width=4.7, fig.cap="Scatterplots of irrigated areas reported by FAO-GMIA against irrigation water withdrawals. Each dot is a country."----

# PLOT -----------------------------------------------------------------------------

full.dt %>%
  na.omit() %>%
  ggplot(., aes(Irrigated.Area, Water.Withdrawn, 
                color = Continent)) +
  geom_point(size = 0.6) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (4 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Irrigated area (ha)", 
       y = expression(paste("Water withdrawal ", " ", "(", 10^9, m^3/year, "", ")"))) +
  facet_wrap(~Water.Dataset, ncol = 4) +
  theme_AP() +
  theme(legend.position = "top")


## ----log10, cache=TRUE, dependson="merge_with_area"------------

# TRANSFORM DATASET ----------------------------------------------------------------

cols <- c("Water.Withdrawn", "Irrigated.Area")
col_names <- c("Continent", "Water.Dataset", "Area.Dataset", "Regression", 
               "Imputation.Method", "Iteration")
cols_group <- c("Continent", "Water.Dataset")
full.dt <- full.dt[, (cols):= lapply(.SD, log10), .SDcols = (cols)]


## ----export_dataset_log10, cache=TRUE, dependson="log10"-------

# EXPORT FULL DATASET WITH MISSING VALUES --------------------------------------------

fwrite(full.dt, "full.dt.csv")
fwrite(water.dt, "water.dt.csv")


## ----plot_missing2, cache=TRUE, dependson="log10", fig.height=3, fig.width=4.5, dev="tikz", fig.cap="Proportion of missing values in irrigation water withdrawal per dataset."----

# PLOT PERCENTAGE OF MISSING 2 -----------------------------------------------------

full.dt[, sum(is.na(.SD) == TRUE) / .N, 
        .(Continent, Water.Dataset), 
        .SDcols = "Water.Withdrawn"] %>%
  ggplot(., aes(Continent, V1, fill = Water.Dataset)) +
  geom_bar(stat = "identity", 
           position = position_dodge(0.7), 
           color = "black") +
  labs(x = "", 
       y = "Proportion of missing values") +
  scale_fill_discrete(name = "") +
  theme_AP() +
  theme(legend.position = "top")


## ----missing, cache=TRUE, dependson="log10", message=FALSE, warning=FALSE----

# IMPUTATION OF MISSING VALUES -----------------------------------------------------

# Substitute Inf values for NA
for (j in 1:ncol(full.dt)) set(full.dt, which(is.infinite(full.dt[[j]])), j, NA)
  
full.dt[, lapply(.SD, function(x) sum(is.infinite(x)))] # Check
  
# Imputation settings
m.iterations <- 40
imputation.methods <- c("norm", "norm.boot", "norm.nob")
  
# Run
imput <- full.dt[, .(Group = lapply(imputation.methods, function(x) 
  mice(.SD, m = m.iterations, maxit = m.iterations, method = x, seed = 500, 
       print = FALSE))), 
  cols_group]
  
imput <- imput[, Imputation.Method:= rep(imputation.methods, .N /
                                           length(imputation.methods))]
  
# Extract iterations
imput <- imput[, Datasets:= lapply(Group, function(x) 
  lapply(1:m.iterations, function(y) data.table(mice::complete(x, y))))] %>%
  .[, Data:= lapply(Datasets, function(x) rbindlist(x, idcol = "Iteration"))]
  
# Vector to loop onto
columns_add <- c("Country", "Iteration", "Irrigated.Area", "Water.Withdrawn")
tmp <- as.list(columns_add)
names(tmp) <- columns_add
  
# Extract columns
for(i in names(tmp)) {
  imput <- imput[, tmp[[i]]:= lapply(.SD, function(x) 
    lapply(x, function(y) y[, ..i])), .SDcols = "Data"]
}
  
# Unlist
full.imput <- imput[, lapply(.SD, unlist), 
                    .SDcols = columns_add, 
                    .(Continent,Water.Dataset, Imputation.Method)]


## ----conduct_lm, cache=TRUE, dependson=c("missing_values", "missing")----

# COMPUTE LINEAR REGRESSIONS --------------------------------------------------------

# Compute regressions in each combination
regressions <- full.imput %>%
  group_by(Continent, Water.Dataset, Imputation.Method, Iteration) %>%
  nest() %>%
  mutate(fit = map(.x = data, .f = ~lm(Water.Withdrawn ~ Irrigated.Area, data = .)), 
         results = map(fit, glance), 
         residuals = map(fit, augment))

# Extract r squared
results <- regressions %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, results) %>%
  unnest(results) %>%
  data.table() %>%
  .[, index:= paste(Continent, Water.Dataset, 
                    Imputation.Method, Iteration, sep = "_")]

# Extract residuals
residuals <- regressions %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, residuals) %>%
  unnest(residuals) %>%
  data.table()


## ----predict, cache=TRUE, dependson="conduct_lm"---------------

# PREDICT WATER WITHDRAWALS --------------------------------------------------------
size.gmia <- meier.dt[, .(Country, Continent, Irrigated.Area)] %>%
  .[!Irrigated.Area == 0] %>%
  .[, Irrigated.Area:= log10(Irrigated.Area)] %>%
  .[!Continent == "Oceania"] %>%
  .[order(Continent)]

countries <- split(size.gmia, size.gmia$Continent) %>%
  lapply(., function(x) x[, Country]) %>%
  lapply(., data.table)

areas <- split(size.gmia, size.gmia$Continent) %>%
  lapply(., function(x) x[, Irrigated.Area]) %>%
  lapply(., data.frame) %>%
  lapply(., function(x) setnames(x, "X..i..", "Irrigated.Area"))

tmp.regressions <- regressions %>%
  split(., .$Continent)

out <- list()
for(i in names(tmp.regressions)) {
  out[[i]] <- mutate(tmp.regressions[[i]], 
                     pred = map(fit, .f = ~predict(., areas[[i]])))
}

water.predicted <- lapply(out, function(x) {
  select(x, Continent, Water.Dataset, Imputation.Method, Iteration, pred) %>%
  unnest(pred) %>%
  data.table()
})

out <- list()
for(i in names(water.predicted)) {
  out[[i]] <- water.predicted[[i]][, Country:= rep(countries[[i]][, V1], 
                                                   times = nrow(water.predicted[[i]]) / 
                                                                  nrow(countries[[i]]))]
}

water.predicted <- rbindlist(water.predicted) %>%
  .[, pred:= 10 ^ pred]

# Compute quantiles
water.quantiles <- water.predicted[, .(min = min(pred), 
                                       max = max(pred), 
                                       q0.025 = quantile(pred, 0.025), 
                                       q0.1 = quantile(pred, 0.1), 
                                       q0.25 = quantile(pred, 0.25),
                                       q0.5 = quantile(pred, 0.5), 
                                       q0.75 = quantile(pred, 0.75), 
                                       q0.975 = quantile(pred, 0.975),
                                       q0.99 = quantile(pred, 0.99), 
                                       q1 = quantile(pred, 1), 
                                       mean = mean(pred), 
                                       median = median(pred)), 
                                   .(Continent, Country)] %>%
  .[order(Country, Continent)]


## ----plot_predicted, cache=TRUE, dependson=c("predict", "conduct_lm", "final_water_dataset"), dev="tikz", fig.height=6, fig.width=5.7, fig.cap="Validation of our approach. The black dots and the error bars show the range of irrigation water withdrawal values predicted from irrigated areas only. The colored dots show the irrigation water withdrawal values outputted by Global Hydrological Models (DBHM, Ho8, LPJmL, MPI-HM, PCR-GLOBWB, WaterGap) and FAO-based datasets (Aquastat, Liu et al. 2016)."----

# PLOT PREDICTIONS AGAINST GHM AND FAO OUTPUTS -------------------------------------

water.tmp <- water.dt[Country %in% water.quantiles[, Country]]
Cont <- c("Africa", "Americas", "Asia", "Europe")
gg <- list()
for(i in Cont) {
  gg[[i]] <- water.quantiles[Continent == i] %>%
    ggplot(., aes(reorder(Country, median), median)) +
    geom_point() +
    geom_point(data = water.tmp[Continent == i], 
               aes(Country, Water.Withdrawn, color = Water.Dataset)) +
    geom_errorbar(aes(ymin = min, 
                      ymax = max)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    scale_color_discrete(name = "Water dataset") +
    labs(y = expression(paste("Water withdrawal ", " ", "(", 10^9, m^3/year, "", ")")), 
         x = "") +
    coord_flip() +
    theme_AP() 
}
gg


## ----lookup, cache=TRUE, dependson="conduct_lm"----------------

# CREATE LOOKUP TABLE --------------------------------------------------------------

lookup <- setkey(results, index)


## ----export_datasets, cache=TRUE, dependson=c("missing_values", "conduct_lm", "lookup", "predict")----

# EXPORT DATASETS -------------------------------------------------------------------

fwrite(full.imput, "full.imput.csv")
fwrite(results, "results.csv")
fwrite(residuals, "residuals.csv")
fwrite(lookup, "lookup.csv")
fwrite(water.quantiles, "water.quantiles.csv")


## ----set_sample_matrix, cache=TRUE-----------------------------

# DEFINE THE SETTINGS OF THE SAMPLE MATRIX -----------------------------------------

Continents <- c("Africa", "Americas", "Asia", "Europe")

# Create a vector with the name of the columns
parameters <- paste("X", 1:3, sep = "")

# Select sample size
n <- 2 ^ 13

# Define order
order <- "third"


## ----sample_matrix, cache=TRUE, dependson="set_sample_matrix"----

# CREATE THE SAMPLE MATRIX ---------------------------------------------------------

# Create an A, B and AB matrices for each continent
sample.matrix <- lapply(Continents, function(Continents) 
  sobol_matrices(N = n,
                 params = parameters,
                 order = order) %>%
    data.table())

# Name the slots, each is a continent
names(sample.matrix) <- Continents

# Name the columns
sample.matrix <- lapply(sample.matrix, setnames, parameters)


## ----add_diagram, echo=FALSE, fig.align="center", fig.cap="Tree diagram coding the discrete probability distributions of each trigger into each uncertainty level.", out.width = '100%'----

# INSERT THE TREE DIAGRAM ----------------------------------------------------------

knitr::include_graphics("./tree_diagram.pdf")


## ----transform_sample_matrix, cache=TRUE, dependson=c("sample_matrix", "set_sample_matrix", "missing_values")----

# TRANSFORM THE SAMPLE MATRIX ------------------------------------------------------

# Function to transform sample matrix to appropriate distributions
transform_sample_matrix <- function(dt) {
  dt[, X1:= floor(X1 * (8 - 1 + 1)) + 1] %>%
    .[, X1:= ifelse(X1 == 1, "LPJmL", 
                    ifelse(X1 == 2, "H08", 
                           ifelse(X1 == 3, "PCR-GLOBWB", 
                                  ifelse(X1 == 4, "WaterGap", 
                                         ifelse(X1 == 5, "Aquastat", 
                                                ifelse(X1 == 6, "Liu et al. 2016", 
                                                       ifelse(X1 == 7, "DBHM", "MPI-HM")))))))] %>%
    .[, X2:= floor(X2 * (length(imputation.methods) - 1 + 1)) + 1] %>%
    .[, X2:= ifelse(X2 == 1, imputation.methods[1], 
                    ifelse(X2 == 2, imputation.methods[2], imputation.methods[3]))] %>%
    .[, X3:= floor(X3 * (m.iterations - 1 + 1)) + 1] 
}

sample.matrix <- lapply(sample.matrix, transform_sample_matrix)
sample.matrix.dt <- rbindlist(sample.matrix, idcol = "Continent")


## ----print_matrix)---------------------------------------------

# PRINT SAMPLE MATRIX --------------------------------------------------------------

print(sample.matrix.dt)


## ----define_model, cache=TRUE----------------------------------

# THE MODEL -------------------------------------------------------------------

model <- function(X) lookup[.(paste0(X[, 1:4], collapse = "_"))][, r.squared]


## ----run_model, cache=TRUE, dependson=c("define_model", "set_boot", "lookup", "transform_sample_matrix", "sample_matrix_gmia")----

# RUN THE MODEL---------------------------------------------------------------------

# Set number of cores at 75%
n_cores <- floor(detectCores() * 0.75)

# Create cluster
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Run model in parallel
r.squared <- foreach(i=1:nrow(sample.matrix.dt),
             .packages = "data.table", 
             .combine = "c") %dopar%
  {
    model(sample.matrix.dt[i])
  }

# Stop parallel cluster
stopCluster(cl)


## ----arrange_output, cache=TRUE, dependson=c("run_model", "transform_sample_matrix")----

# ARRANGE MODEL OUTPUT -------------------------------------------------------------

sample.matrix.dt <- cbind(sample.matrix.dt, r.squared) 

# Select the A and B matrix only (for uncertainty analysis)
AB.dt <- sample.matrix.dt[, .SD[1:(n * 2)], Continent]

# Export results
fwrite(sample.matrix.dt, "sample.matrix.dt.csv")
fwrite(AB.dt, "AB.dt.csv")


## ----quantiles, cache=TRUE, dependson="arrange_output"---------

# COMPUTE QUANTILES AND MEAN -------------------------------------------------------

AB.dt[, .(q0.025 = quantile(r.squared, 0.025), 
          q0.1 = quantile(r.squared, 0.1), 
          q0.25 = quantile(r.squared, 0.25),
          q0.5 = quantile(r.squared, 0.5), 
          q0.75 = quantile(r.squared, 0.75), 
          q0.975 = quantile(r.squared, 0.975),
          q0.99 = quantile(r.squared, 0.99), 
          q1 = quantile(r.squared, 1), 
          mean = mean(r.squared), 
          median = median(r.squared)), Continent]


## ----plot_uncertainty, cache=TRUE, dependson="arrange_output", dev="tikz", fig.height=5, fig.width=2.5, fig.cap="Uncertainty in the empirical distribution of $r^2$."----

# PLOT UNCERTAINTY -----------------------------------------------------------------

# Plot r2
unc.plot <- ggplot(AB.dt, aes(r.squared)) + 
      geom_histogram(color = "black", fill = "white") + 
      theme_AP() +
      labs(x = expression(italic(r) ^ 2), 
           y = "Density") +
      scale_y_continuous(breaks = pretty_breaks(n = 2)) +
      scale_x_continuous(breaks = pretty_breaks(n = 3)) +
      facet_wrap(~Continent, ncol = 1) +
      theme(panel.spacing.x = unit(4, "mm"))

unc.plot


## ----plot_uncertainty_GHM, cache=TRUE, dependson="arrange_output"----

# PLOT UNCERTAINTY IN EACH GHM AND FAO-BASED DATASET -------------------------------

unc.GHM <- AB.dt %>%
  ggplot(., aes(reorder(X1, r.squared), r.squared, fill = Continent)) +
  geom_boxplot(position = position_dodge(0.6), 
               outlier.size = 0.3) +
  theme_AP() +
  labs(y = expression(italic(r)^2), 
       x = "") +
  scale_x_discrete(labels = c("PCR-GLOBWB" = expression(bold(PCR-GLOBWB)), 
                              "DBHM" = expression(bold(DBHM)), 
                              "MPI-HM" = expression(bold(MPI-HM)), 
                              "LPJmL" = expression(bold(LPJmL)), 
                              "H08" = expression(bold(H08)), 
                              "WaterGap" = expression(bold(WaterGap)))) +
  theme(legend.position = "none") + 
  coord_flip()

# Get legend
legend <- get_legend(unc.GHM + theme(legend.position = "top"))


## ----merge_plots, cache=TRUE, dependson=c("plot_uncertainty", "plot_uncertainty_GHM", "run_model"), fig.height=4, fig.width=5.2, dev="tikz"----

# MERGE PLOTS ----------------------------------------------------------------------

bottom <- plot_grid(unc.plot, unc.GHM, ncol = 2, 
                    rel_widths = c(0.5, 1), labels = "auto")
all <- plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.1, 1))


## ----cumulative_r2, cache=TRUE, dependson="arrange_output", dev = "tikz", fig.height=2.2, fig.width=3.3, fig.cap="Cumulative empirical distribution for $r^2$."----


# PLOT CUMULATIVE EMPIRICAL DISTRIBUTION FOR R2 ------------------------------------

ggplot(AB.dt, aes(r.squared, colour = Continent)) + 
  stat_ecdf() +
  labs(x = "$r^2$", 
       y = "y") + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP() 


## ----scatterplots, cache=TRUE, dependson="arrange_output", fig.height=7, fig.width=5.3, fig.cap="Scatterplots of $r^2$ against the triggers' levels."----

# PLOT SCATTERPLOTS OF PARAMETERS VS MODEL OUTPUT ----------------------------------

AB.dt <- AB.dt[, X3:= factor(X3, levels = as.factor(1:m.iterations))]

scatter.dt <- melt(AB.dt[, .SD[1:n], Continent], 
                   measure.vars = paste("X", 1:3, sep = "")) 

# R squared
ggplot(scatter.dt, aes(r.squared, value)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_grid(variable ~ Continent,
             scales = "free_y") +
  labs(y = "", 
       x = expression(italic(r)^2)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_AP() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top")


## ----sensitivity, cache=TRUE, dependson=c("arrange_output", "set_sample_matrix")----

# SENSITIVITY ANALYSIS -------------------------------------------------------------

# Number of bootstrap replicas
R <- 1000

parameters.recoded <- c("$X_1$", "$X_2$", "$X_3$")

# Sobol' indices for r2
indices <- sample.matrix.dt[, sobol_indices(Y = r.squared, 
                                            N = n,
                                            params = parameters.recoded, 
                                            first = "jansen",
                                            R = R, 
                                            boot = TRUE,
                                            parallel = "multicore", 
                                            ncpus = n_cores, 
                                            order = order), 
                            Continent]


## ----print_sensitivity, cache=TRUE, dependson="sensitivity"----

# PRINT AND EXPORT SENSITIVITY INDICES ---------------------------------------------

print(indices[sensitivity %in% c("Si", "Ti")])
fwrite(indices, "indices.csv")


## ----plot_sobol, cache=TRUE, dependson="sensitivity", dev="tikz", fig.width = 4.7, fig.height=2, fig.cap="Sobol' indices. $S_i$ and $T_i$ refer respectively to Sobol' first and total order indices. $S_i$ measures the influence of a parameter in the model output, while $T_i$ measures the influence of a parameter jointly with its interactions."----

# PLOT UNCERTAINTY AND SOBOL' INDICES ----------------------------------------------

bottom <- indices[sensitivity %in% c("Si", "Ti")] %>%
  ggplot(., aes(parameters, original, fill = sensitivity)) +
    geom_bar(stat = "identity",
             position = position_dodge(0.6),
             color = "black") +
    geom_errorbar(aes(ymin = low.ci,
                      ymax = high.ci),
                  position = position_dodge(0.6)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    facet_wrap(~Continent,
               ncol = 4) +
    labs(x = "",
         y = "Sobol' index") +
    scale_fill_discrete(name = "Sobol' indices",
                        labels = c(expression(S[italic(i)]),
                                   expression(T[italic(i)]))) +
    theme_AP() +
    theme(legend.position = "top")

bottom


## ----merge_all_plots, cache=TRUE, dependson=c("plot_sobol", "merge_plots"), dev="tikz", fig.height=5.7, fig.width=5.2, fig.cap="Uncertainty and sensitivity analysis. a) Empirical distribution for $r^2$ at the continental level. b) Boxplots of $r^2$ values obtained when regressions where run with GHM (in bold) and FAO-based datasets. c) Sobol' indices. $S_i$ and $T_i$ refer respectively to Sobol' first and total order indices. $S_i$ measures the influence of a parameter in the model output, while $T_i$ measures the influence of a parameter jointly with its interactions."----

# MERGE UNCERTAINTY AND SENSITIVITY ANALYSIS PLOTS ---------------------------------

plot_grid(all, bottom, align = "hv", rel_heights = c(0.8, 0.4), 
          labels = c("", "c"), ncol = 1)


## ----sum_si, cache=TRUE, dependson="sensitivity"---------------

# CHECK SUM OF FIRST-ORDER INDICES -------------------------------------------------

indices[sensitivity == "Si",  sum(original), Continent]


## ----plot_sobol_second_third, cache=TRUE, dependson="sensitivity", dev = "tikz", fig.height = 2.3, fig.width=4.3, fig.cap="High-order interactions between the triggers. The dots and the errorbars show the mean and the 95\\% confidence intervals after bootstraping (R=1000)."----

# PLOT SOBOL' INDICES (SECOND AND THIRD ORDER) -------------------------------------

indices[sensitivity == "Sij" | sensitivity == "Sijk"] %>%
  .[low.ci > 0] %>%
  .[, sensitivity:= ifelse(sensitivity %in% "Sij", "$S_{ij}$", 
                           "$S_{ijk}$")] %>%
  ggplot(., aes(reorder(parameters, original), original, color = Continent)) +
  geom_point(position = position_dodge(0.6)) +
  geom_errorbar(aes(ymax = high.ci, ymin = low.ci),
                position = position_dodge(0.6)) +
  geom_hline(yintercept = 0,
             lty = 2,
             color = "red") +
  facet_grid(~sensitivity, 
             scales = "free_x", 
             space = "free_x") +
  scale_color_discrete(name = "Continent") +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  labs(x = "",
       y = "Sobol' index") +
  theme_AP()


## ----session_information---------------------------------------

# SESSION INFORMATION --------------------------------------------------------------

sessionInfo()

## Return the machine CPU
cat("Machine: "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores: "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = TRUE))

## Return the machine RAM
cat("RAM: "); print (get_ram()); cat("\n")

