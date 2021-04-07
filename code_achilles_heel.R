## ----setup, include=FALSE------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, dev = "tikz")


## ----load_packages, results="hide", message=FALSE, warning=FALSE---------------------------

# LOAD PACKAGES -----------------------------------------------------------------

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
               "benchmarkme", "tidyverse", "grid", "gridExtra",
               "robustbase", "rsample"))

# SET CHECKPOINT --------------------------------------------------------------

dir.create(".checkpoint")

library("checkpoint")

checkpoint("2021-04-07", 
           R.version ="4.0.3", 
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


## ----functions_data, cache=TRUE------------------------------------------------------------

# CREATE FUNCTIONS ---------------------------------------------------------------

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
  out <- lapply(selected.years, function(x) {
    water <- rowSums(ww[, vec[[x]]]) 
    ww.df <- data.frame(cbind(lon, lat, water)) 
    ww.df <- data.frame(cbind(lon, lat, water)) 
    countries <- coords2country(ww.df[1:nrow(ww.df), 1:2])
    df <- cbind(countries, ww.df)
    setDT(df)
    final <- df[, .(Water.Withdrawn = sum(water)), countries]
    setnames(final, "countries", "Country")
    country_code(final)
    out <- na.omit(final[order(Continent)])
    out[, Water.Withdrawn:= Water.Withdrawn / 1000] # From mm to m
  })
  names(out) <- selected.years
  return(out)
}

# Function to load and extract data from .nc files from ISIMIP
open_nc_files <- function(file, dname, selected.years, vec) {
  ncin <- nc_open(file)
  # get longitude, latitude, time
  lon <- ncvar_get(ncin, "lon")
  lat <- ncvar_get(ncin, "lat")
  # Get variable
  tmp_array <- ncvar_get(ncin, dname)
  m <- lapply(selected.years, function(x) vec[[x]])
  
  out <- lapply(m, function(x) {
    tmp_slice <- lapply(x, function(y) tmp_array[, , y])
    # create dataframe -- reshape data
    # matrix (nlon*nlat rows by 2 cols) of lons and lats
    lonlat <- as.matrix(expand.grid(lon,lat))
    # vector of `tmp` values
    tmp_vec <- lapply(tmp_slice, function(x) as.vector(x))
    # create dataframe and add names
    tmp_df01 <- lapply(tmp_vec, function(x) data.frame(cbind(lonlat, x)))
    names(tmp_df01) <- x
    da <- lapply(tmp_df01, data.table) %>%
      rbindlist(., idcol = "month") %>%
      na.omit()
    # Convert coordinates to country
    Country <- coords2country(da[1:nrow(da), 2:3])
    df <- cbind(Country, da)
    setDT(df)
    out <- na.omit(df)[, .(Water.Withdrawn = sum(x)), Country]
    out[, Water.Withdrawn:= Water.Withdrawn * 10000]
    out[, Continent:= countrycode(out[, Country], 
                                  origin = "country.name", 
                                  destination = "continent")] %>%
      .[, Code:= countrycode(out[, Country], 
                             origin = "country.name", 
                             destination = "un")] %>%
      .[, Country:= countrycode(out[, Code], 
                                origin = "un", 
                                destination = "country.name")] %>%
      .[!Continent == "Oceania"]
    setcolorder(out, c("Country", "Continent", "Code", "Water.Withdrawn"))
  })
  return(out)
}


## ----ghm_datasets, cache=TRUE, dependson="functions_data"----------------------------------

# READ IN GHM DATASETS -----------------------------------------------------------

# Read in Huang et al. datasets -------------------------------

# Define vector with the months and the years
vecs <- 1:((2010 - 1970) * 12)
vec <- split(vecs, ceiling(seq_along(vecs) / 12))
names(vec) <- 1971:2010
selected.years <- c("1971", "1980", "1990", "2000", "2005")

# Vector with the files
names_nc_files <- c("withd_irr_lpjml.nc", "withd_irr_pcrglobwb.nc", 
                    "withd_irr_h08.nc", "withd_irr_watergap.nc")

# Read in datasets
huang.dt <- mclapply(names_nc_files, function(x) get_nc_data(x), mc.cores = detectCores() * 0.75)
names.huang <- c("LPJmL", "PCR-GLOBWB", "H08", "WaterGap")
names(huang.dt) <- names.huang

# ISIMIP datasets ----------------------------------------------

files <- list("dbh_wfdei_nobc_hist_varsoc_co2_airrww_global_monthly_1971_2010.nc", 
              "mpi-hm_miroc5_ewembi_picontrol_histsoc_co2_airrww_global_monthly_1861_2005.nc", 
              "vic_wfdei_nobc_hist_pressoc_co2_airrww_global_monthly_1971_2010.nc", 
              "clm45_gfdl-esm2m_ewembi_historical_2005soc_co2_pirrww_global_monthly_1861_2005.nc")

isimip.dt <- mclapply(files, function(x) {
  if(x == "mpi-hm_miroc5_ewembi_picontrol_histsoc_co2_airrww_global_monthly_1861_2005.nc") {
    vecs <- 1:((2005 - 1860) * 12)
    vec <- split(vecs, ceiling(seq_along(vecs) / 12))
    names(vec) <- 1861:2005
    dname <- "airrww"
  } else if (x == "clm45_gfdl-esm2m_ewembi_historical_2005soc_co2_pirrww_global_monthly_1861_2005.nc") {
    vecs <- 1:((2005 - 1860) * 12)
    vec <- split(vecs, ceiling(seq_along(vecs) / 12))
    names(vec) <- 1861:2005
    dname <- "pirrww"
  } else {
    vecs <- 1:((2010 - 1970) * 12)
    vec <- split(vecs, ceiling(seq_along(vecs) / 12))
    names(vec) <- 1971:2010
    dname <- "airrww"
  }
  open_nc_files(file = x, dname = dname, selected.years = selected.years, 
                vec = vec) 
}, 
mc.cores = detectCores() * 0.75)

names.isimip <- c("DBHM", "MPI-HM", "VIC", "CLM45")
names(isimip.dt) <- names.isimip

for(i in names(isimip.dt)) {
  names(isimip.dt[[i]]) <- selected.years
}  


## ----historic_fao_gmia, cache=TRUE, dependson="ghm_datasets"-------------------------------

# READ IN HISTORIC FAO-GMIA ------------------------------------------------------

gmia.hist <- fread("FAO_GMIA_historic.csv") 
gmia.hist <- gmia.hist[, Continent:= countrycode(gmia.hist[, Country], 
                                                 origin = "country.name", 
                                                 destination = "continent")] %>%
  .[, Code:= countrycode(gmia.hist[, Country], 
                         origin = "country.name", 
                         destination = "un")] %>%
  .[, Country:= countrycode(gmia.hist[, Code], 
                            origin = "un", 
                            destination = "country.name")] %>%
  .[!Continent == "Oceania"]

# Create list with the columns selected

dt.tmp <- lapply(c("AEI_1970", paste("AEI_", c(seq(1980, 2000, 10), 2005), sep = "")), function(x)
  gmia.hist[, .SD, .SDcols = c("Country", "Code", "Continent", x)])
names(dt.tmp) <- selected.years

# MERGE GHM AND HISTORIC FAO-GMIA ------------------------------------------------

all.ghm <- c(huang.dt, isimip.dt)
ghm.dt <- lapply(names(all.ghm), function(x)
  lapply(selected.years, function(y)
    all.ghm[[x]][[y]][, merge(.SD, dt.tmp[[y]], by = c("Country", "Code", "Continent"), 
                              all.y = TRUE)])) %>%
  lapply(., function(x) lapply(x, function(y) setnames(y, 5, "AEI") %>%
                                 na.omit()))

names(ghm.dt) <- c(names.huang, names.isimip)

for(i in names(ghm.dt)) {
  names(ghm.dt[[i]]) <- selected.years
}     


## ----plot_historic_fao_gmia, cache=TRUE, dependson="historic_fao_gmia", fig.height=4.3, fig.width=5.4----

# PLOT GHM WATER WITHDRAWALS THROUGH TIME (1970-2005) ------------------------------

tmp <- lapply(ghm.dt, function(x) rbindlist(x, idcol = "Year"))

# Plot
gg <- list()
for(i in 1:length(tmp)) {
  gg[[i]] <- ggplot(tmp[[i]], aes(AEI, Water.Withdrawn, 
                                  color = Continent)) +
    geom_point(size = 0.4) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    labs(x = "Irrigated area (ha)", 
         y = expression(paste("Water withdrawal ", " ", "(", 10^9, m^3/year, "", ")"))) +
    facet_wrap(~Year) +
    theme_AP() +
    theme(legend.position = "top", 
          strip.background = element_rect(fill = "white")) +
    ggtitle(names(tmp[i]))
}

gg


## ----long_term_trends, cache=TRUE, dependson=c("functions_data", "historic_fao_gmia")------

# ANALYZE LONG-TERM HISTORICAL TRENDS --------------------------------------------

files <- list("pcr-globwb_gfdl-esm2m_ewembi_historical_histsoc_co2_airrww_global_monthly_1861_2005.nc", 
              "mpi-hm_gfdl-esm2m_ewembi_historical_histsoc_co2_airrww_global_monthly_1861_2005.nc", 
              "lpjml_miroc5_ewembi_historical_histsoc_co2_airrww_global_monthly_1861_2005.nc", 
              "h08_miroc5_ewembi_historical_histsoc_co2_airrww_global_monthly_1861_2005.nc")

vecs <- 1:((2005 - 1860) * 12)
vec <- split(vecs, ceiling(seq_along(vecs) / 12))
names(vec) <- 1861:2005
dname <- "airrww"
selected.years <- as.character(seq(1900, 2000, 10))

isimip.hist <- mclapply(files, function(x)
  open_nc_files(file = x, dname = dname, selected.years = selected.years, 
                vec = vec), mc.cores = detectCores() * 0.75)

names.isimip.hist <- c("PCR-GLOBWB", "MPI-HM", "LPJmL", "H08")
names(isimip.hist) <- names.isimip.hist

for(i in names(isimip.hist)) {
  names(isimip.hist[[i]]) <- selected.years
}  

# Create list with the columns selected

dt.tmp <- lapply(paste("AEI_", seq(1900, 2000, 10), sep = ""), function(x)
  gmia.hist[, .SD, .SDcols = c("Country", "Code", "Continent", x)])
names(dt.tmp) <- selected.years


out <- lapply(names(isimip.hist), function(x)
  lapply(selected.years, function(y)
    isimip.hist[[x]][[y]][, merge(.SD, dt.tmp[[y]], by = c("Country", "Code", "Continent"), 
                              all.y = TRUE)])) %>%
  lapply(., function(x) lapply(x, function(y) setnames(y, 5, "AEI") %>%
                                 na.omit()))

names(out) <- c("PCR-GLOBWB", "MPI-HM", "LPJmL", "H08")

for(i in names(out)) {
  names(out[[i]]) <- selected.years
}     


## ----plot_long_term_trends, cache=TRUE, dependson="long_term_trends", fig.width=5.4, fig.height=4.3----

# PLOT GHM WATER WITHDRAWALS THROUGH TIME (1900-1960) ------------------------------

tmp <- lapply(out, function(x) rbindlist(x, idcol = "Year")) %>%
  na.omit()

# Plot
gg <- list()
for(i in 1:length(tmp)) {
  gg[[i]] <- ggplot(tmp[[i]], aes(AEI, Water.Withdrawn, 
                                  color = Continent)) +
    geom_point(size = 0.4) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    labs(x = "Irrigated area (ha)", 
         y = expression(paste("Water withdrawal ", " ", "(", 10^9, m^3/year, "", ")"))) +
    facet_wrap(~Year, ncol = 5) +
    theme_AP() +
    theme(legend.position = "top", 
          strip.background = element_rect(fill = "white")) +
    ggtitle(names(tmp[i]))
}

gg


## ----data_water, cache=TRUE, dependson="historic_fao_gmia"---------------------------------

# READ IN DATASETS ON IRRIGATION WATER WITHDRAWAL AND MEIER DATASET --------------

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
  country_code(.) %>%
  setcolorder(., cols_order)

# Liu et al. dataset ----------------------------
#UNIT IS 10^9 m3/year = km3/year
liu.dt <- fread("liu.csv")[, .(country, irr)] %>%
  setnames(., c("country", "irr"), c("Country", "Water.Withdrawn")) %>%
  country_code(.) %>%
  .[, Water.Dataset:= "Liu et al. 2016"] %>%
  country_code(.) %>%
  setcolorder(., cols_order) 

# READ IN FAO-GMIA BY MEIER  ------------------------------------------------------

meier.dt <- fread("meier.csv") %>%
  setnames(., "Codes", "Code") %>%
  na.omit() %>%
  .[, .(Country, Continent, Code, `FAO-GMIA`)] %>%
  setnames(., "FAO-GMIA", "Irrigated.Area")

# MERGE WITH MEIER ET AL. FAO-GMIA -----------------------------------------------

water.dt <- lapply(all.ghm, function(x) rbindlist(x, idcol = "Year")) %>%
  rbindlist(., idcol = "Water.Dataset") %>%
  .[Year == 2005] %>%
  .[, Year:= NULL] %>%
  setcolorder(., cols_order) %>%
  rbind(., table4.dt, liu.dt)

# Check if there are duplicated countries
water.dt[, .N, .(Country)][N > 10][, Country]

# Compute the mean
water.dt <- water.dt[, .(Water.Withdrawn = mean(Water.Withdrawn)), 
                     .(Continent, Country, Code, Water.Dataset)]

# And check again
water.dt[, .N, .(Country)][N > 10][, Country]

full.dt <- water.dt[, merge(.SD, meier.dt, by = c("Country", "Code", "Continent"), 
                       all.y = TRUE), Water.Dataset] %>%
  .[!Continent == "Oceania"]

# Show countries with missing values in Water Withdrawal
full.dt[is.na(Water.Withdrawn), ] %>%
  .[, unique(Country), Water.Dataset] %>%
  print(n = Inf)

# Show unique countries missing
full.dt[is.na(Water.Withdrawn), ] %>%
  .[, unique(Country)]


## ----plot_merged, cache=TRUE, dependson=c("data_water"), fig.height=3, fig.width=5.4, fig.cap="Scatterplots of irrigated areas reported by FAO-GMIA against irrigation water withdrawals. Each dot is a country."----

# PLOT ---------------------------------------------------------------------------

full.dt %>%
  na.omit() %>%
  ggplot(., aes(Irrigated.Area, Water.Withdrawn, 
                color = Continent)) +
  geom_point(size = 0.4) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (4 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Irrigated area (ha)", 
       y = expression(paste("Water withdrawal ", " ", "(", 10^9, m^3/year, "", ")"))) +
  facet_wrap(~Water.Dataset, ncol = 5) +
  theme_AP() +
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "white")) 


## ----plot_divided, cache=TRUE, dependson="data_water"--------------------------------------

# CHECK THE INFLUENCE OF POTENTIAL EVAPOTRANSPIRATION ----------------------------

full.dt <- full.dt[, div:= Water.Withdrawn / Irrigated.Area] 

# Create function to access the evapotranspiration data from ISIMIP
open_nc_totevap <- function(file, dname) {
ncin <- nc_open(file)
# get longitude, latitude, time
lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
# Get variable
tmp_array <- ncvar_get(ncin, dname)
m <- (dim(tmp_array)[3] - 11):dim(tmp_array)[3]
tmp_slice <- lapply(m, function(m) tmp_array[, , m])
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
dt <- cbind(countries, da)
setDT(dt)
dt <- na.omit(dt)[, .(totevap = sum(x)), countries]
setnames(dt, "countries", "Country")
dt[, `:=` (Code = countrycode(dt[, Country], 
                              origin = "country.name", 
                              destination = "un"), 
           Continent = countrycode(dt[, Country], 
                                   origin = "country.name", 
                                   destination = "continent"))]
dt[, Country:= countrycode(dt[, Code], 
                           origin = "un", 
                           destination = "country.name")]
return(dt)
}

# List of files to access
files <- list("watergap2_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc",
              "pcr-globwb_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc",
              "mpi-hm_wfdei_nobc_hist_pressoc_co2_evap_global_monthly_1971_2010.nc",
              "lpjml_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc",
              "h08_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc",
              "vic_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc", 
              "dbh_wfdei_nobc_hist_varsoc_co2_evap_global_monthly_1971_2010.nc")

# Retrieve and arrange the data
dname <- "evap"
totevap.dt <- mclapply(files, function(x) 
  open_nc_totevap(file = x, 
                  dname = dname), 
  mc.cores = detectCores() * 0.75)
gm_list <- c("WaterGap", "PCR-GLOBWB", "MPI-HM", "LPJmL", "H08", "VIC", "DBHM")
names(totevap.dt) <- gm_list


## ----plot_evapotranspiration, cache=TRUE, dependson="plot_divided", fig.height=3, fig.width=5.4----

# PLOT EVAPOTRANSPIRATION -------------------------------------------------------

rbindlist(totevap.dt, idcol = "Water.Dataset") %>%
  na.omit() %>%
  merge(full.dt[Water.Dataset %in% gm_list], ., by = c("Continent", "Country", 
                                                       "Code", "Water.Dataset")) %>%
  ggplot(., aes(totevap, div, 
                color = Continent)) +
  geom_point(size = 0.6) +
  scale_x_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Total evapotranspiration (kg $\\textrm{m}^{-2} \\quad \\textrm{s}^{-1}$)", 
       y = "$\\frac{\\textrm{Irrigation Water Withdrawal}}{\\textrm{Irrigated area}}$") +
  facet_wrap(~Water.Dataset, ncol = 4) +
  theme_AP() +
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "white")) 


## ----check_potevap, cache=TRUE, dependson = "plot_divided"---------------------------------

# CHECK THE INFLUENCE OF POTENTIAL EVAPORATION ----------------------------------

files <- c("watergap2_wfdei_nobc_hist_pressoc_co2_potevap_global_monthly_1971_2010.nc", 
           "h08_wfdei_nobc_hist_pressoc_co2_potevap_global_monthly_1971_2010.nc", 
           "pcr-globwb_wfdei_nobc_hist_pressoc_co2_potevap_global_monthly_1971_2010.nc", 
           "mpi-hm_wfdei_nobc_hist_pressoc_co2_potevap_global_monthly_1971_2010.nc", 
           "lpjml_wfdei_nobc_hist_pressoc_co2_potevap_global_monthly_1971_2010.nc")

# Retrieve and arrange the data
dname <- "potevap"
potevap.dt <- mclapply(files, function(x) 
  open_nc_totevap(file = x, 
                  dname = dname), 
  mc.cores = detectCores() * 0.75)
gm_list <- c("WaterGap", "H08", "PCR-GLOBWB", "MPI-HM", "LPJmL")
names(potevap.dt) <- gm_list


## ----plot_potevap, cache=TRUE, dependson="check_potevap", fig.height=2, fig.width=5.4------

# PLOT POTENTIAL EVAPORATION ----------------------------------------------------

rbindlist(potevap.dt, idcol = "Water.Dataset") %>%
  na.omit() %>%
  merge(full.dt[Water.Dataset %in% gm_list], ., by = c("Continent", "Country", 
                                                       "Code", "Water.Dataset")) %>%
  ggplot(., aes(totevap, div, 
                color = Continent)) +
  geom_point(size = 0.6) +
  scale_x_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Potential evaporation (kg $\\textrm{m}^{-2} \\quad  \\textrm{s}^{-1}$)", 
       y = "$\\frac{\\textrm{Irrigation Water Withdrawal}}{\\textrm{Irrigated area}}$") +
  facet_wrap(~Water.Dataset, ncol = 5) +
  theme_AP() +
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "white")) 


## ----influence_efficiency, cache=TRUE, dependson=c("plot_divided", "data_water"), fig.height=2.3, fig.width=5.4----

# CHECK INFLUENCE OF WATER EFFICIENCY --------------------------------------------

# Retrieve the data from Rohwer
rohwer_data <- fread("rohwer_data.csv")

# Plot
merge(full.dt[Water.Dataset %in% c("WaterGap", "H08", "LPJmL", 
                                   "PCR-GLOBWB")], 
      rohwer_data, by = "Country") %>%
  na.omit() %>%
  ggplot(., aes(Project.efficiency, div, 
                color = Continent)) +
  geom_point(size = 0.6) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_x_log10(breaks = pretty_breaks(n = 3)) +
  labs(x = "Water efficiency", 
       y = "$\\frac{\\textrm{Irrigation Water Withdrawal}}{\\textrm{Irrigated area}}$") +
  facet_wrap(~Water.Dataset, ncol = 4) +
  theme_AP() +
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "white")) 


## ----log10, cache=TRUE, dependson="data_water"---------------------------------------------

# TRANSFORM DATASET -------------------------------------------------------------

cols <- c("Water.Withdrawn", "Irrigated.Area")
col_names <- c("Continent", "Water.Dataset", "Area.Dataset", "Regression", 
               "Imputation.Method", "Iteration")
cols_group <- c("Continent", "Water.Dataset")
full.dt <- full.dt[, (cols):= lapply(.SD, log10), .SDcols = (cols)]
full.dt <- full.dt[, Country:= str_replace(Country, "\\&", "and")]


## ----diagnostic, cache=TRUE, dependson=c("data_water", "log10"), fig.width=4.5, fig.height=7----

# REGRESSION DIAGNOSTICS --------------------------------------------------------

# Conduct regressions
regression.diag <- full.dt %>%
  NaRV.omit() %>%
  group_by(Continent, Water.Dataset) %>%
  nest() %>%
  mutate(fit = map(.x = data, .f = ~lm(Water.Withdrawn ~ Irrigated.Area, data = .)), 
         results = map(fit, glance), 
         residuals = map(fit, augment))

# Prepare dataset
diagnostics.dt <- regression.diag %>%
  dplyr::select(Continent, Water.Dataset, residuals) %>%
  unnest(residuals) %>%
  data.table() 

# Residuals versus fitted
ggplot(diagnostics.dt, aes(.fitted, .resid)) +
  geom_point(size = 0.5) + 
  geom_hline(yintercept = 0, lty = 2, color = "red") +
  geom_smooth(size = 0.5) +
  labs(x = "Fitted value", y = "Residuals") +
  facet_grid(Water.Dataset ~ Continent, 
             scales = "free_y") +
  theme_AP() +
  theme(strip.text.y = element_text(size = 5.5), 
        strip.background = element_rect(fill = "white"))

# QQ plot
ggplot(diagnostics.dt, aes(sample = .std.resid)) +
  stat_qq(size = 0.5) + 
  stat_qq_line() +
  facet_grid(Water.Dataset ~ Continent, 
             scales = "free_y") +
  labs(x = "Theoretical quantiles", 
       y = "Standardized Residuals") +
  theme_AP() +
  theme(strip.text.y = element_text(size = 5.5), 
        strip.background = element_rect(fill = "white"))

# Scale location plot
ggplot(diagnostics.dt, aes(.fitted, sqrt(abs(.std.resid)))) +
  geom_point(size = 0.5) + 
  geom_smooth() +
  facet_grid(Water.Dataset ~ Continent, 
             scales = "free_y") + 
  labs(x = "Fitted value", 
       y = "$\\sqrt{|\\textrm{Standardized residuals}|}$") +
  theme_AP() +
  theme(strip.text.y = element_text(size = 5.5), 
        strip.background = element_rect(fill = "white"))

# Cook's distance
ID <- diagnostics.dt[, .N, .(Continent, Water.Dataset)]
diagnostics.dt <- diagnostics.dt[, ID:= unlist(sapply(ID[, N], function(x) 1:x))]
ggplot(diagnostics.dt, aes(ID, .cooksd)) +
  geom_col() +
  facet_grid(Water.Dataset ~ Continent, 
             scales = "free_x") +
  theme_AP() +
  labs(x = "Observation number", y = "Cook's distance") +
  theme(strip.text.y = element_text(size = 5.5), 
        strip.background = element_rect(fill = "white"))


## ----export_dataset_log10, cache=TRUE, dependson="log10"-----------------------------------

# EXPORT FULL DATASET WITH MISSING VALUES ---------------------------------------

fwrite(full.dt, "full.dt.csv")
fwrite(water.dt, "water.dt.csv")


## ----plot_missing2, cache=TRUE, dependson="log10", fig.height=3, fig.width=4.5, fig.cap="Proportion of missing values in irrigation water withdrawal per dataset."----

# PLOT PERCENTAGE OF MISSING 2 --------------------------------------------------

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


## ----missing, cache=TRUE, dependson="log10", message=FALSE, warning=FALSE------------------

# IMPUTATION OF MISSING VALUES ----------------------------------------------------

# Substitute Inf values for NA
for (j in 1:ncol(full.dt)) set(full.dt, which(is.infinite(full.dt[[j]])), j, NA)

full.dt[, lapply(.SD, function(x) sum(is.infinite(x)))] # Check

# Imputation settings
m.iterations <- 40
imputation.methods <- c("norm.boot", "norm", "norm.nob")

# Run
full.dt <- full.dt[, div:= NULL]

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


## ----conduct_lm, cache=TRUE, dependson=c("missing_values", "missing")----------------------

# COMPUTE LINEAR REGRESSIONS -----------------------------------------------------

# Compute regressions in each combination
# Settings for robust regressions
a1<-lmrob.control()
a1$k.max <- 1000

# Compute regressions
regressions <- full.imput %>%
  group_by(Continent, Water.Dataset, Imputation.Method, Iteration) %>%
  nest() %>%
  mutate(fit = map(.x = data, .f = ~lm(Water.Withdrawn ~ Irrigated.Area, data = .)), 
         fit.robust = map(.x = data, .f = ~lmrob(Water.Withdrawn ~ Irrigated.Area, data = ., 
                                                 control = a1)),
         results = map(fit, glance), 
         results.robust = map(fit.robust, glance),
         residuals = map(fit, augment))

# EXTRACT R SQUARED ----------------------
# Regular r squared
results <- regressions %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, results,) %>%
  unnest(results) %>%
  data.table() %>%
  .[, Regression:= "Normal"] %>%
  .[, index:= paste(Continent, Water.Dataset, Imputation.Method, 
                    Iteration, Regression, sep = "_")]

# Robust r squared
results.robust <- regressions %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, results.robust) %>%
  unnest(results.robust) %>%
  data.table() %>%
  .[, Regression:= "Robust"] %>%
  .[, index:= paste(Continent, Water.Dataset, Imputation.Method, 
                    Iteration, Regression, sep = "_")]

cols_extract <- c("Continent", "Water.Dataset", "Imputation.Method", 
                  "Iteration", "Regression", "r.squared", 
                  "index")

# Bind regular and robust
all.results <- rbind(results[, ..cols_extract], 
                     results.robust[, ..cols_extract])


## ----obtain_beta, cache=TRUE, dependson="conduct_lm", fig.width=5, fig.height=6.5----------

# EXTRACT SLOPE (BETA) -----------------------------------------------------------

# Extract regular slope
slope <- regressions %>%
  mutate(slope = map(fit, tidy)) %>%
  unnest(slope) %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, term, estimate) %>%
  data.table() %>%
  .[term == "Irrigated.Area"]

# Extract robust slope
slope.robust <- regressions %>%
  mutate(slope = map(fit.robust, tidy)) %>%
  unnest(slope) %>%
  dplyr::select(Continent, Water.Dataset, 
                Imputation.Method, Iteration, term, estimate) %>%
  data.table() %>%
  .[term == "Irrigated.Area"]

# Plot
rbind(slope, slope.robust) %>%
  ggplot(., aes(estimate)) +
  geom_histogram(color = "black", aes(fill = estimate > 1)) +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "$\\beta$", 
       y = "Count") +
  scale_fill_discrete(name = "$\\beta$", 
                      labels = c("$<1$", "$>1$")) +
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  facet_grid(Water.Dataset~Continent) +
  theme_AP() + 
  theme(strip.text.y = element_text(size = 5),
        strip.background = element_rect(fill = "white"))


## ----predict, cache=TRUE, dependson=c("conduct_lm", "data_water")--------------------------

# PREDICT WATER WITHDRAWALS ----------------------------------------------------

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

out <- out.robust <- list()
for(i in names(tmp.regressions)) {
  out[[i]] <- mutate(tmp.regressions[[i]], 
                     pred = map(fit, .f = ~predict(., areas[[i]])))
  out.robust[[i]] <- mutate(tmp.regressions[[i]], 
                            pred = map(fit.robust, .f = ~predict(., areas[[i]])))
}

water.predicted <- lapply(out, function(x) {
  select(x, Continent, Water.Dataset, Imputation.Method, Iteration, pred) %>%
  unnest(pred) %>%
  data.table()
})

water.predicted.rob <- lapply(out.robust, function(x) {
  select(x, Continent, Water.Dataset, Imputation.Method, Iteration, pred) %>%
    unnest(pred) %>%
    data.table()
})

out <- out.robust <- list()
for(i in names(water.predicted)) {
  out[[i]] <- water.predicted[[i]][, Country:= rep(countries[[i]][, V1], 
                                                   times = nrow(water.predicted[[i]]) / 
                                                                  nrow(countries[[i]]))]
  out.robust[[i]] <- water.predicted.rob[[i]][, Country:= rep(countries[[i]][, V1], 
                                                          times = nrow(water.predicted[[i]]) / 
                                                            nrow(countries[[i]]))]
}

water.predicted <- rbindlist(water.predicted) %>%
  .[, pred:= 10 ^ pred]

water.predicted.rob <- rbindlist(water.predicted.rob) %>%
  .[, pred:= 10 ^ pred]

# Compute quantiles
water.quantiles <- rbind(water.predicted, water.predicted.rob) %>%
  .[, .(min = min(pred), 
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

water.quantiles <- water.quantiles[, Country:= str_replace(Country, "\\&", "and")]


## ----plot_predicted, cache=TRUE, dependson=c("predict", "conduct_lm", "data_water"), fig.height=6.5, fig.width=5.7, fig.cap="Validation of our approach. The black dots and the error bars show the range of irrigation water withdrawal values predicted from irrigated areas only. The colored dots show the irrigation water withdrawal values outputted by Global Hydrological Models (DBHM, Ho8, LPJmL, MPI-HM, PCR-GLOBWB, WaterGap) and FAO-based datasets (Aquastat, Liu et al. 2016)."----

# PLOT PREDICTIONS AGAINST GHM AND FAO OUTPUTS ----------------------------------

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

all.plots <- lapply(1:4, function(x) 
  gg[[x]] + 
    theme(legend.position = "none") +
    labs(x = "", y = ""))

legend <- get_legend(gg[[1]] + theme(legend.position = "top", 
                                     legend.key.size = unit(0.5, 'lines')))
bottom1 <- plot_grid(all.plots[[1]], all.plots[[2]], ncol = 2, labels = "auto", align = "hv")

all1 <- plot_grid(legend, bottom1, ncol = 1, rel_heights = c(0.1, 1))

grid.arrange(arrangeGrob(all1, bottom = grid::textGrob(
  label = expression(paste("Irrigation water withdrawal ", " ", "(", 10^9, m^3/year, "", ")")))))

bottom2 <- plot_grid(all.plots[[3]], all.plots[[4]], ncol = 2, labels = "auto", align = "hv")
all <- plot_grid(legend, bottom2, ncol = 1, rel_heights = c(0.1, 1))

grid.arrange(arrangeGrob(all, bottom = grid::textGrob(
  label = expression(paste("Irrigation water withdrawal ", " ", "(", 10^9, m^3/year, "", ")")))))


## ----check_estimates, cache=TRUE, dependson=c("predict", "plot_predicted"),fig.height=3, fig.width=4----

# CHECK HOW MANY ESTIMATES FIT WITHIN THE REGRESSION BOUNDS ---------------------

da <- merge(water.tmp, water.quantiles[, .(Continent, Country, min, max)], 
            by = c("Continent", "Country")) %>%
  .[, fit:= ifelse(Water.Withdrawn >= min & Water.Withdrawn <= max, TRUE, FALSE)]

# Check which GHM or FAO-based dataset shows more estimates
# beyond or below our predictions
da[, .(N = sum(Water.Withdrawn < min | Water.Withdrawn > max), 
       prop = sum(Water.Withdrawn < min | Water.Withdrawn > max) / .N),
   .(Water.Dataset, Continent)]

# Plot bars
prove <- da[, sum(fit), .(Country, Continent)]

ggplot(prove, aes(V1)) +
  geom_bar() +
  facet_wrap(~Continent) + 
  theme_AP() + 
  labs(x = "Number of estimates", 
       y = "Number of countries") +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  theme(strip.background = element_rect(fill = "white"))

# Check how many countries show more than 7 or 10 estimates
# bounded by our predictions
lapply(c(7, 10), function(x) 
  prove[, .(Total.countries = .N, 
            Bounded = sum(V1 >= x), 
            Proportion = sum(V1 >= x) / .N)])

# Which countries do not show any framed estimate, or only 1, 2 or 3
lapply(c(0:3), function(x) prove[V1 == x])


## ----all_together, cache=TRUE, dependson="check_estimates", fig.height=2, fig.width=2.3----

# PLOT ALL ESTIMATES TOGETHER ---------------------------------------------------

data.table(table(prove$V1)) %>%
  .[, V1:= factor(V1, levels = 0:10)] %>%
  ggplot(., aes(V1, N)) +
  geom_point() +
  labs(x = "Number of point estimates", 
       y = "Number of countries") +
  theme_AP()


## ----check_gh_fao, cache=TRUE, dependson="check_estimates", fig.height=5, fig.width=4.7----

# CHECK BIASED POINT ESTIMATES BY DATASET ---------------------------------------

da[, .(N = sum(Water.Withdrawn < min | Water.Withdrawn > max), 
       prop = sum((Water.Withdrawn < min | Water.Withdrawn > max) / .N) * 100),
   .(Water.Dataset, Continent)] %>%
  ggplot(., aes(Water.Dataset, prop, fill = Continent)) +
  geom_bar(stat = "identity", 
           position = position_dodge(0.7), 
           color = "black") +
  labs(x = "", 
       y = "Percent (\\%)") +
  coord_flip() + 
  theme_AP() +
  theme(legend.position = "top")


## ----future_irrigation, cache=TRUE---------------------------------------------------------

# FUTURE IRRIGATION WATER WITHDRAWALS --------------------------------------------

files <- list(
  "pcr-globwb_miroc5_ewembi_rcp60_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "pcr-globwb_miroc5_ewembi_rcp26_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "lpjml_miroc5_ewembi_rcp60_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "lpjml_miroc5_ewembi_rcp26_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "lpjml_miroc5_ewembi_rcp60_rcp60soc_co2_airrww_global_monthly_2006_2099.nc",
  "lpjml_miroc5_ewembi_rcp26_rcp26soc_co2_airrww_global_monthly_2006_2099.nc",
  "h08_miroc5_ewembi_rcp60_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "h08_miroc5_ewembi_rcp26_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "h08_miroc5_ewembi_rcp60_rcp60soc_co2_airrww_global_monthly_2006_2099.nc",
  "h08_miroc5_ewembi_rcp26_rcp26soc_co2_airrww_global_monthly_2006_2099.nc",
  "mpi-hm_miroc5_ewembi_rcp60_2005soc_co2_airrww_global_monthly_2006_2099.nc",
  "mpi-hm_miroc5_ewembi_rcp26_2005soc_co2_airrww_global_monthly_2006_2099.nc"
)

vecs <- 1:((2099 - 2005) * 12)
vec <- split(vecs, ceiling(seq_along(vecs) / 12))
names(vec) <- 2006:2099
dname <- "airrww"
selected.years <- as.character(seq(2030, 2050, 10))

# Read in datasets
isimip.future <- mclapply(
  files, function(x)
    open_nc_files(file = x, dname = dname, selected.years = selected.years, 
                  vec = vec), mc.cores = detectCores() * 0.75
)


## ----plot_future_irrigation, cache=TRUE, dependson="future_irrigation", fig.height=7, fig.width=5.7----

# PLOT ---------------------------------------------------------------------------

# Create vector to name the slots
ghms <- c(c("PCR-GLOBWB", "PCR-GLOBWB"), rep(c("LPJmL", "H08"), each = 4), c("MPI-HM", "MPI-HM"))
climate_scenario <- rep(c(60, 26), 6)
climate2 <- c(rep(2005, 4), rep(c(60, 26, 2005, 2005), 2))
names.isimip.future <- paste(paste(ghms, climate_scenario, sep = "/"), climate2, sep = ".")

# Name the slots
names(isimip.future) <- names.isimip.future

for(i in names(isimip.future)) {
  names(isimip.future[[i]]) <- selected.years
}  

# Arrange data
isimip.future.dt <- lapply(isimip.future, function(x) rbindlist(x, idcol = "Year")) %>%
  rbindlist(., idcol = "Model") %>%
  .[!Continent == "Oceania"] %>%
  separate(., "Model", c("Model", "Climate scenario"), "/") %>%
  na.omit()

Cont <- unique(isimip.future.dt$Continent)
water.tmp <- isimip.future.dt[Country %in% data.table(water.quantiles)[, Country]]

gg <- list()
for(j in Cont) {
  gg[[j]] <- water.quantiles[Continent == j] %>%
    ggplot(., aes(reorder(Country, median), median)) +
    geom_point() +
    geom_point(data = water.tmp[Year == "2050" & Continent == j], 
               aes(reorder(Country, Water.Withdrawn), Water.Withdrawn, 
                   shape = `Climate scenario`, 
                   color = Model)) +
    geom_errorbar(aes(ymin = min, 
                      ymax = max)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                  labels = trans_format("log10", math_format(10 ^ .x))) +
    scale_color_discrete(name = "Water dataset") +
    labs(y = "", x = "") +
    coord_flip() +
    theme_AP() +
    theme(legend.position = "none")
}

legend <- get_legend(gg[[1]] +
                       theme(legend.position = "top", 
                             legend.box = "vertical"))

bottom1 <- plot_grid(gg[[1]] + theme(legend.position = "none"), 
                     gg[[2]] + theme(legend.position = "none"), ncol = 2)
all1 <- plot_grid(legend, bottom1, rel_heights = c(0.2, 1), ncol = 1)
grid.arrange(arrangeGrob(all1, bottom = grid::textGrob(
  label = expression(paste("Irrigation water withdrawal ", " ", "(", 10^9, m^3/year, "", ")")))))

bottom2 <- plot_grid(gg[[3]] + theme(legend.position = "none"), 
                     gg[[4]] + theme(legend.position = "none"), ncol = 2)
all2 <- plot_grid(legend, bottom2, rel_heights = c(0.2, 1), ncol = 1)
grid.arrange(arrangeGrob(all2, bottom = grid::textGrob(
  label = expression(paste("Irrigation water withdrawal ", " ", "(", 10^9, m^3/year, "", ")")))))


## ----lookup, cache=TRUE, dependson="conduct_lm"--------------------------------------------

# CREATE LOOKUP TABLE -----------------------------------------------------------

lookup <- setkey(all.results, index)


## ----export_datasets, cache=TRUE, dependson=c("missing_values", "conduct_lm", "lookup", "predict")----

# EXPORT DATASETS ----------------------------------------------------------------

fwrite(full.imput, "full.imput.csv")
fwrite(results, "results.csv")
fwrite(lookup, "lookup.csv")
fwrite(water.quantiles, "water.quantiles.csv")


## ----set_sample_matrix, cache=TRUE---------------------------------------------------------

# DEFINE THE SETTINGS OF THE SAMPLE MATRIX ---------------------------------------

Continents <- c("Africa", "Americas", "Asia", "Europe")

# Create a vector with the name of the columns
parameters <- paste("X", 1:4, sep = "")

# Select sample size
n <- 2 ^ 13

# Define order
order <- "third"


## ----sample_matrix, cache=TRUE, dependson="set_sample_matrix"------------------------------

# CREATE THE SAMPLE MATRIX -------------------------------------------------------

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

# INSERT THE TREE DIAGRAM --------------------------------------------------------

knitr::include_graphics("./tree_diagram.pdf")


## ----transform_sample_matrix, cache=TRUE, dependson=c("sample_matrix", "set_sample_matrix", "missing_values")----

# TRANSFORM THE SAMPLE MATRIX ----------------------------------------------------

# Function to transform sample matrix to appropriate distributions
transform_sample_matrix <- function(dt) {
  dt[, X1:= floor(X1 * (10 - 1 + 1)) + 1] %>%
    .[, X1:= ifelse(X1 == 1, "LPJmL", 
                    ifelse(X1 == 2, "H08", 
                           ifelse(X1 == 3, "PCR-GLOBWB", 
                                  ifelse(X1 == 4, "WaterGap", 
                                         ifelse(X1 == 5, "Aquastat", 
                                                ifelse(X1 == 6, "Liu et al. 2016", 
                                                       ifelse(X1 == 7, "DBHM", 
                                                              ifelse(X1 == 8, "VIC", 
                                                                     ifelse(X1 == 9, "MPI-HM", "CLM45")))))))))] %>%
    .[, X2:= floor(X2 * (length(imputation.methods) - 1 + 1)) + 1] %>%
    .[, X2:= ifelse(X2 == 1, imputation.methods[1], 
                    ifelse(X2 == 2, imputation.methods[2], imputation.methods[3]))] %>%
    .[, X3:= floor(X3 * (m.iterations - 1 + 1)) + 1] %>%
    .[, X4:= floor(X4 * (2 - 1 + 1)) + 1] %>%
    .[, X4:= ifelse(X4 == 1, "Normal", "Robust")]
}

sample.matrix <- lapply(sample.matrix, transform_sample_matrix)
sample.matrix.dt <- rbindlist(sample.matrix, idcol = "Continent")


## ----print_matrix)-------------------------------------------------------------------------

# PRINT SAMPLE MATRIX ------------------------------------------------------------

print(sample.matrix.dt)


## ----define_model, cache=TRUE--------------------------------------------------------------

# THE MODEL ----------------------------------------------------------------------

model <- function(X) lookup[.(paste0(X[, 1:5], collapse = "_"))][, r.squared]


## ----run_model, cache=TRUE, dependson=c("define_model", "set_boot", "lookup", "transform_sample_matrix", "sample_matrix")----

# RUN THE MODEL-------------------------------------------------------------------

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


## ----arrange_output, cache=TRUE, dependson=c("run_model", "transform_sample_matrix")-------

# ARRANGE MODEL OUTPUT ------------------------------------------------------------

sample.matrix.dt <- cbind(sample.matrix.dt, r.squared) 

# Select the A and B matrix only (for uncertainty analysis)
AB.dt <- sample.matrix.dt[, .SD[1:(n * 2)], Continent]

# Export results
fwrite(sample.matrix.dt, "sample.matrix.dt.csv")
fwrite(AB.dt, "AB.dt.csv")


## ----quantiles, cache=TRUE, dependson="arrange_output"-------------------------------------

# COMPUTE QUANTILES AND MEAN ------------------------------------------------------

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


## ----plot_uncertainty, cache=TRUE, dependson="arrange_output", fig.height=5, fig.width=2.5, fig.cap="Uncertainty in the empirical distribution of $r^2$."----

# PLOT UNCERTAINTY --------------------------------------------------------------

# Plot r2
unc.plot <- ggplot(AB.dt, aes(r.squared)) + 
      geom_histogram(color = "black", fill = "white") + 
      theme_AP() +
      labs(x = expression(italic(r) ^ 2), 
           y = "Count") +
      scale_y_continuous(breaks = pretty_breaks(n = 2)) +
      scale_x_continuous(breaks = pretty_breaks(n = 3)) +
      facet_wrap(~Continent, ncol = 1) +
      theme(panel.spacing.x = unit(4, "mm"), 
            strip.background = element_rect(fill = "white"))

unc.plot


## ----plot_uncertainty_GHM, cache=TRUE, dependson="arrange_output"--------------------------

# PLOT UNCERTAINTY IN EACH GHM AND FAO-BASED DATASET ----------------------------

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
                              "WaterGap" = expression(bold(WaterGap)), 
                              "CLM45" = expression(bold(CLM45)), 
                              "VIC" = expression(bold(VIC)))) +
  theme(legend.position = "none") + 
  coord_flip()

# Get legend
legend <- get_legend(unc.GHM + theme(legend.position = "top"))


## ----merge_plots, cache=TRUE, dependson=c("plot_uncertainty", "plot_uncertainty_GHM", "run_model"), fig.height=4, fig.width=5.2----

# MERGE PLOTS -------------------------------------------------------------------

bottom <- plot_grid(unc.plot, unc.GHM, ncol = 2, 
                    rel_widths = c(0.5, 1), labels = "auto")
all <- plot_grid(legend, bottom, ncol = 1, rel_heights = c(0.1, 1))


## ----cumulative_r2, cache=TRUE, dependson="arrange_output", fig.height=2.2, fig.width=3.3, fig.cap="Cumulative empirical distribution for $r^2$."----

# PLOT CUMULATIVE EMPIRICAL DISTRIBUTION FOR R2 -----------------------------------

ggplot(AB.dt, aes(r.squared, colour = Continent)) + 
  stat_ecdf() +
  labs(x = "$r^2$", 
       y = "y") + 
  scale_y_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP() 


## ----scatterplots, cache=TRUE, dependson="arrange_output", fig.height=7, fig.width=5.3, fig.cap="Scatterplots of $r^2$ against the triggers' levels.", dev="pdf"----

# PLOT SCATTERPLOTS OF PARAMETERS VS MODEL OUTPUT ---------------------------------

AB.dt <- AB.dt[, X3:= factor(X3, levels = as.factor(1:m.iterations))]

scatter.dt <- melt(AB.dt[, .SD[1:n], Continent], 
                   measure.vars = paste("X", 1:4, sep = "")) 

# R squared
ggplot(scatter.dt, aes(r.squared, value)) +
  geom_point(alpha = 0.05, size = 0.5) +
  facet_grid(variable ~ Continent,
             scales = "free_y", 
             space = "free_y") +
  labs(y = "", 
       x = expression(italic(r)^2)) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_color_manual(values = c("#00BFC4", "#F8766D")) +
  theme_AP() +
  theme(axis.text.y = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top")


## ----sensitivity, cache=TRUE, dependson=c("arrange_output", "set_sample_matrix")-----------

# SENSITIVITY ANALYSIS ----------------------------------------------------------

# Number of bootstrap replicas
R <- 10^3

parameters.recoded <- c("$X_1$", "$X_2$", "$X_3$", "$X_4$")

# Sobol' indices for r2
indices <- sample.matrix.dt[, sobol_indices(Y = r.squared, 
                                            N = n,
                                            params = parameters.recoded, 
                                            first = "jansen",
                                            R = R, 
                                            boot = TRUE,
                                            parallel = "multicore", 
                                            ncpus = n_cores, 
                                            order = order)$results, 
                            Continent]


## ----print_sensitivity, cache=TRUE, dependson="sensitivity"--------------------------------

# EXPORT SENSITIVITY INDICES -----------------------------------------------------

fwrite(indices, "indices.csv")


## ----plot_sobol, cache=TRUE, dependson="sensitivity", dev="tikz", fig.width = 4.7, fig.height=2, fig.cap="Sobol' indices. $S_i$ and $T_i$ refer respectively to Sobol' first and total order indices. $S_i$ measures the influence of a parameter in the model output, while $T_i$ measures the influence of a parameter jointly with its interactions."----

# PLOT UNCERTAINTY AND SOBOL' INDICES -------------------------------------------

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
    theme(legend.position = "top", 
          strip.background = element_rect(fill = "white"))

bottom


## ----merge_all_plots, cache=TRUE, dependson=c("plot_sobol", "merge_plots"), dev="tikz", fig.height=5.7, fig.width=5.2, fig.cap="Uncertainty and sensitivity analysis. a) Empirical distribution for $r^2$ at the continental level. b) Boxplots of $r^2$ values obtained when regressions where run with GHM (in bold) and FAO-based datasets. c) Sobol' indices. $S_i$ and $T_i$ refer respectively to Sobol' first and total order indices. $S_i$ measures the influence of a parameter in the model output, while $T_i$ measures the influence of a parameter jointly with its interactions."----

# MERGE UNCERTAINTY AND SENSITIVITY ANALYSIS PLOTS --------------------------------

plot_grid(all, bottom, align = "hv", rel_heights = c(0.8, 0.4), 
          labels = c("", "c"), ncol = 1)


## ----sum_si, cache=TRUE, dependson="sensitivity"-------------------------------------------

# CHECK SUM OF FIRST-ORDER INDICES ----------------------------------------------

indices[sensitivity == "Si",  sum(original), Continent]


## ----plot_sobol_second_third, cache=TRUE, dependson="sensitivity", dev = "tikz", fig.height = 2.3, fig.width=4.3, fig.cap="High-order interactions between the triggers. The dots and the errorbars show the mean and the 95\\% confidence intervals after bootstraping (R=1000)."----

# PLOT SOBOL' INDICES (SECOND AND THIRD ORDER) ------------------------------------

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


## ----australia.dt, cache=TRUE--------------------------------------------------------------

# RETRIEVE DATA FROM DIFFERENT SCALES ---------------------------------------------

# Australia irrigation schemes -----------------------------------------

australia.scheme <- fread("australia_scheme.csv")

# Filter out NA and rows with 0
australia.dt <- australia.scheme[, lapply(.SD, function(x) 
  ifelse(x == 0, NA, x))] %>%
  na.omit() %>%
  .[, Year:= factor(Year, levels = c("97.98", "2002.2003"))] %>%
  .[, Year:= gsub("\\.", "-", Year)] %>%
  .[, Scale:= "Australian \n irrigation systems"]

# USA level ------------------------------------------------------------

# Irrigated area is in thousand acres
# Water withdrawals in thousand acres feet

# To convert from acre feet to cubic meters, multiply IWW by 1233.48.
# To convert from acre to ha, divide the area value by 2.471.

col_transf <- c("Irrigated.Area", "Water.Withdrawal")

solley.dt <- fread("solley.dt.csv") %>%
  .[, Scale:= "USA states"] # state level

colorado.dt <- fread("colorado_data.csv") # colorado
colorado.dt <- colorado.dt[, `:=`(Irrigated.Area = Flood + Sprinkler + Micro, 
                                  Water.Withdrawal = Groundwater + Surface.water)] %>%
  .[, .SD, .SDcols = (col_transf)] %>%
  .[, Scale:= "Colorado counties"]

# Merge both data sets
usa.dt <- rbind(solley.dt, colorado.dt) %>%
  .[, Irrigated.Area:= (Irrigated.Area * 1000) / 2.471] %>%
  .[, Water.Withdrawal:= Water.Withdrawal * 1000 * 1233.48]


## ----plot_australian, cache=TRUE, dependson=australia.dt, fig.height=2, fig.width=4--------

# PLOT SCATTERPLOT ---------------------------------------------------------------

scatter_usa <- ggplot(usa.dt, aes(Irrigated.Area, Water.Withdrawal)) +
  geom_point(size = 0.6) +
  scale_x_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  labs(x = "Irrigated area (ha)", 
       y = "Irrigation water \n Withdrawals (m$^3$)") +
  facet_wrap(~Scale) +
  theme_AP() +
  theme(strip.background = element_rect(fill = "white"))

scatter_australia <- ggplot(australia.dt, aes(Irrigated.Area, Water.Withdrawal)) +
  geom_point(size = 0.6) +
  scale_x_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10 ^ .x))) +
  facet_grid(Scale~Year) +
  labs(x = "Irrigated area (ha)", 
       y = "Irrigation deliveries \n (ML/year)") +
  theme_AP() +
  theme(strip.background = element_rect(fill = "white")) 

all_scatters <- plot_grid(scatter_australia, scatter_usa, ncol = 2)
all_scatters


## ----australian_regressions, cache=TRUE, dependson="australia.dt", fig.height=2, fig.width=3----

# COMPUTE REGRESSIONS AND BOOTSTRAP ----------------------------------------------

# AUSTRALIA IRRIGATION SYSTEMS -----------------------------------

australia.dt <- australia.dt[, (col_transf):= lapply(.SD, log10), .SDcols = col_transf]

australia.regressions <- australia.dt %>%
  group_by(Year) %>%
  nest() %>%
  mutate(fit = map(.x = data, .f = ~lm(Water.Withdrawal~ Irrigated.Area, 
                                      data = .)), 
         results = map(fit, glance), 
         residuals = map(fit, augment))

# Bootstrap
foo <- boot(australia.dt, function(data,indices)
  summary(lm(Water.Withdrawal ~ Irrigated.Area,
             data[indices, ] ))$r.squared, R = 5000)

# Confidence intervals
ci_plot <- boot.ci(foo, type = "all") 

# USA LEVEL (COLORADO COUNTY AND COUNTY LEVEL) -------------------

usa.dt.lm <- usa.dt[, (col_transf):= lapply(.SD, log10), .SDcols = col_transf] %>%
  NaRV.omit() 

usa.regressions <- usa.dt.lm %>%
  group_by(Scale) %>%
  nest() %>%
  mutate(fit = map(.x = data, .f = ~lm(Water.Withdrawal~ Irrigated.Area, 
                                       data = .)), 
         results = map(fit, glance), 
         residuals = map(fit, augment))

# Bootstrap
usa.loop <- c("USA states", "Colorado counties")
foo.usa <- lapply(usa.loop, function(x)
  boot(usa.dt.lm[Scale == x], function(data,indices)
    summary(lm(Water.Withdrawal ~ Irrigated.Area,
               data[indices, ] ))$r.squared, R = 5000))

names(foo.usa) <- usa.loop

# Confidence intervals
ci_plot2 <- lapply(foo.usa, function(x) boot.ci(x, type = "all"))

ci.dt <- data.table(Scale = c("USA states", "Colorado counties", "Australian irrigation systems"), 
                    low.ci = c(ci_plot2$`USA states`$bca[[4]], 
                               ci_plot2$`Colorado counties`$bca[[4]], 
                               ci_plot$bca[[4]]), 
                    high.ci = c(ci_plot2$`USA states`$bca[[5]], 
                                ci_plot2$`Colorado counties`$bca[[5]], 
                                ci_plot$bca[[5]])) %>%
  melt(., measure.vars = c("low.ci", "high.ci"))


## ----plot_all_australia, cache=TRUE, dependson=c("australian_regressions", "plot_australian"), fig.height=4.5, fig.width=5.6----

# PLOT ALL -----------------------------------------------------------------------

tmp <- lapply(usa.loop, function(x) foo.usa[[x]]$t) %>%
  lapply(., data.table)

names(tmp) <- usa.loop

all.histo <- rbindlist(tmp, idcol = "Scale") %>%
  rbind(., data.table(foo$t) %>%
          .[, Scale:= "Australian irrigation systems"] %>% 
          setcolorder(., c("Scale", "V1"))) %>%
  .[, Scale:= factor(Scale, levels = c("Colorado counties", "USA states", 
                                       "Australian irrigation systems"))] %>%
  ggplot(., aes(V1, group = Scale)) +
  geom_histogram(color = "black", fill = "white") +
  geom_vline(data = ci.dt, aes(xintercept = value), color = "blue") +
  labs(x = "$r^2$", 
       y = "Count") +
  facet_wrap(~Scale) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP() +
  theme(strip.background = element_rect(fill = "white")) 

plot_grid(all_scatters, all.histo, ncol = 1, labels = "auto")


## ----australian_diagnostics, cache=TRUE, dependson="australian_regressions", fig.height=5, fig.width=3.5----

# REGRESSION DIAGNOSTICS ---------------------------------------------------------

# AUSTRALIA ---------------------------------------
australia.residuals <- australia.regressions  %>%
  dplyr::select(Year, residuals) %>%
  unnest(residuals) %>%
  data.table() 

# USA  -------------------------------------------
usa.residuals <- usa.regressions  %>%
  dplyr::select(Scale, residuals) %>%
  unnest(residuals) %>%
  data.table() 

# functions to plot:

# residuals versus fitted
plot_res <- function(dt) {
  gg <- ggplot(dt, aes(.fitted, .resid)) +
    geom_point(size = 0.5) + 
    geom_hline(yintercept = 0, lty = 2, color = "red") +
    geom_smooth(size = 0.5) +
    labs(x = "Fitted value", y = "Residuals") +
    theme_AP() +
    theme(strip.text.y = element_text(size = 6.4), 
          strip.background = element_rect(fill = "white"))
  return(gg)
}

# qq plot
plot_qq <- function(dt) {
  gg <- ggplot(dt, aes(sample = .std.resid)) +
    stat_qq(size = 0.5) + 
    stat_qq_line() +
    labs(x = "Theoretical quantiles", 
         y = "Standardized Residuals") +
    theme_AP() +
    theme(strip.text.y = element_text(size = 6.4), 
          strip.background = element_rect(fill = "white"))
  return(gg)
}

# PLOT FOR AUSTRALIA ------------------

a <- plot_res(australia.residuals) +
  facet_grid(~Year, scales = "free_y") 
  
b <- plot_qq(australia.residuals) +
  facet_grid(~Year, scales = "free_y")

ID <- australia.residuals[, .N, Year]
numb <- unlist(sapply(ID[, N], function(x) 1:x))
c <- cbind(australia.residuals, numb) %>%
  ggplot(., aes(numb, .cooksd)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(~Year, scales = "free_x") +
  theme_AP() +
  labs(x = "Observation number", y = "Cook's distance") +
  theme(strip.text.y = element_text(size = 6.4), 
        strip.background = element_rect(fill = "white"))

plot_grid(a, b, c, labels = "auto", ncol = 1, align = "hv", hjust = -2.5)

# PLOT FOR USA ------------------------

a <- plot_res(usa.residuals) +
  facet_grid(~Scale, scales = "free_y") 

b <- plot_qq(usa.residuals) +
  facet_grid(~Scale, scales = "free_y")

ID <- usa.residuals[, .N, Scale]
numb <- unlist(sapply(ID[, N], function(x) 1:x))
c <- cbind(usa.residuals, numb) %>%
  ggplot(., aes(numb, .cooksd)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 1)) +
  facet_grid(~Scale, scales = "free_x") +
  theme_AP() +
  labs(x = "Observation number", y = "Cook's distance") +
  theme(strip.text.y = element_text(size = 6.4), 
        strip.background = element_rect(fill = "white"))

plot_grid(a, b, c, labels = "auto", ncol = 1, align = "hv", hjust = -2.5)


## ----boot_beta, cache=TRUE, dependson=c("australia.dt", "australian_regressions"), fig.height=2, fig.width=2.5----

# BOOTSTRAP BETA  ---------------------------------------------------------------

# Australia
foo.beta <- boot(australia.dt, function(data,indices)
  coef(lm(Water.Withdrawal ~ Irrigated.Area,
          data[indices, ]))[[2]], R = 5000)

# USA
foo.states <- boot(usa.dt.lm[Scale == "USA states"], function(data,indices)
  coef(lm(Water.Withdrawal ~ Irrigated.Area,
          data[indices, ]))[[2]], R = 5000)

foo.counties <- boot(usa.dt.lm[Scale == "Colorado counties"], function(data,indices)
  coef(lm(Water.Withdrawal ~ Irrigated.Area,
          data[indices, ]))[[2]], R = 5000)

# Confidence intervals
ci_plot2 <- lapply(list(foo.states,
                        foo.counties,
                        foo.beta), function(x) boot.ci(x, type = "all"))

names(ci_plot2) <- c("USA states", "Colorado counties", "Australian irrigation systems")

ci.dt <- data.table(Scale = c("USA states", "Colorado counties", "Australian irrigation systems"), 
                    low.ci = c(ci_plot2$`USA states`$bca[[4]], 
                               ci_plot2$`Colorado counties`$bca[[4]], 
                               ci_plot2$`Australian irrigation systems`$bca[[4]]), 
                    high.ci = c(ci_plot2$`USA states`$bca[[5]], 
                                ci_plot2$`Colorado counties`$bca[[5]], 
                                ci_plot2$`Australian irrigation systems`$bca[[5]])) %>%
  melt(., measure.vars = c("low.ci", "high.ci"))

tmp <- rbind(data.table(foo.states$t)[, Scale:= "USA states"], 
             data.table(foo.counties$t)[, Scale:= "Colorado counties"], 
             data.table(foo.beta$t)[, Scale:= "Australian irrigation systems"])


## ----plot_beta, cache=TRUE, dependson="boot_beta", fig.height=2.5--------------------------

# PLOT BETA ----------------------------------------------------------------------

# Plot
  ggplot(tmp, aes(V1)) +
  geom_histogram(color = "black", fill = "white") +
  geom_vline(data = ci.dt, aes(xintercept = value), color = "blue") +
  geom_vline(xintercept = 1, lty = 2) +
  labs(x = "$\\beta$", y = "Count") +
  facet_wrap(~Scale) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  theme_AP() +
  theme(strip.background = element_rect(fill = "white")) 


## ----get_all_pixel, cache=TRUE-------------------------------------------------------------

# PIXEL LEVEL --------------------------------------------------------------------

ncin <- nc_open("histsoc_landuse-totals_annual_1861_2005.nc")
lon <- ncvar_get(ncin, "lon")
lat <- ncvar_get(ncin, "lat")
tmp_array <- ncvar_get(ncin, "cropland_irrigated")
m <- (dim(tmp_array)[3] - 11):dim(tmp_array)[3]
tmp_slice <- lapply(m, function(m) tmp_array[, , m])
lonlat <- as.matrix(expand.grid(lon,lat))
tmp_vec <- lapply(tmp_slice, function(x) as.vector(x))
tmp_df01 <- lapply(tmp_vec, function(x) data.frame(cbind(lonlat, x)))
names(tmp_df01) <- m
da <- lapply(tmp_df01, data.table) %>%
  rbindlist(., idcol = "month") %>%
  na.omit() %>%
  .[!x == 0]
Country <- coords2country(da[1:nrow(da), 2:3])
df <- cbind(Country, da) 
setDT(df)
df <- setnames(df, c("Var1", "Var2"), c("lon", "lat"))
mirca <- na.omit(df)[, .(area = mean(x)), .(Country, lon, lat)] %>%
  .[order(Country)]

# Function to open the Huang datasets
open_nc_mirca <- function(file, dname) {
  nc <- nc_open(file)
  ww <- ncvar_get(nc, dname)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  water <- rowSums(ww[, 469:ncol(ww)]) 
  ww.df <- data.frame(cbind(lon, lat, water)) 
  setDT(ww.df)
  ww.df <- ww.df[!water == 0]
  Country <- coords2country(ww.df[1:nrow(ww.df), 1:2])
  dt <- cbind(Country, ww.df) %>%
    na.omit() %>%
    .[, water:= water / 1000]
  return(dt)
}

# Function to open ISIMIP datasets
open_nc_mirca_isimip <- function(file, dname) {
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
  Country <- coords2country(da[1:nrow(da), 2:3])
  df <- cbind(Country, da)
  setDT(df)
  out <- na.omit(df)[, .(water = sum(x)), .(Country, Var1, Var2)] %>%
    .[, water:= water * 10000]
  out <- out[order(Country)] %>%
    .[!water == 0] %>%
    setnames(., c("Var1", "Var2"), c("lon", "lat"))
  return(out)
}

# HUANG ET AL DATASETS
names_nc_files <- c("withd_irr_lpjml.nc", "withd_irr_pcrglobwb.nc", 
                    "withd_irr_h08.nc", "withd_irr_watergap.nc")
out.nc.mirca <- mclapply(names_nc_files, function(x) 
  open_nc_mirca(x, "withd_irr"), mc.cores = detectCores() * 0.75)
names(out.nc.mirca) <- c("LPJmL", "PCR-GLOBWB", "H08", "WaterGap")
out.final.mirca <- rbindlist(out.nc.mirca, idcol = "Water.Dataset")

# ISIMIP DATASETS
files <- list("dbh_wfdei_nobc_hist_varsoc_co2_airrww_global_monthly_1971_2010.nc", 
              "mpi-hm_miroc5_ewembi_picontrol_histsoc_co2_airrww_global_monthly_1861_2005.nc", 
              "vic_wfdei_nobc_hist_pressoc_co2_airrww_global_monthly_1971_2010.nc")
dname <- "airrww"

isimip.dt.mirca <- mclapply(files, function(x) 
  open_nc_mirca_isimip(file = x, dname = dname), mc.cores = detectCores() * 0.75)

names(isimip.dt.mirca) <- c("DBHM", "MPI-HM", "VIC")

# ADD CLM45
CLM45.mirca <- open_nc_mirca_isimip(file = "clm45_gfdl-esm2m_ewembi_historical_2005soc_co2_pirrww_global_monthly_1861_2005.nc", 
                                    dname = "pirrww")

CLM45.mirca <- CLM45.mirca[, Water.Dataset:= "CLM45"] %>%
  setcolorder(., c("Water.Dataset", "Country", "lon", "lat", "water"))

final.isimip.mirca <- rbindlist(isimip.dt.mirca, idcol = "Water.Dataset") %>%
  rbind(CLM45.mirca) %>%
  rbind(out.final.mirca) %>%
  na.omit()

full.isimip <- merge(final.isimip.mirca, mirca, by = c("Country", "lon", "lat"), all.y = TRUE) 

full.isimip[, `:=` (Code = countrycode(full.isimip[, Country], 
                                       origin = "country.name", 
                                       destination = "un"), 
                    Continent = countrycode(full.isimip[, Country], 
                                            origin = "country.name", 
                                            destination = "continent"))]

full.isimip <- full.isimip[!Continent == "Oceania"]


## ----plot_pixel, cache=TRUE, dependson="get_all_pixel", fig.height=7, fig.width=6----------

# PLOT ---------------------------------------------------------------------------

full.isimip.split <- split(full.isimip, full.isimip$Continent) %>%
  lapply(., function(x) split(x, x$Water.Dataset))

gg <- list()
for(i in names(full.isimip.split)) {
  for(j in names(full.isimip.split[[i]])) {
    gg[[i]][[j]] <- ggplot(full.isimip.split[[i]][[j]], aes(area, water)) +
      geom_point(size = 0.1, alpha = 0.5) +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                    labels = trans_format("log10", math_format(10 ^ .x))) +
      scale_y_log10(breaks = trans_breaks("log10", function(x) 10 ^ (2 * x)),
                    labels = trans_format("log10", math_format(10 ^ .x))) +
      facet_wrap(~Country, 
                 ncol = 5) +
      labs(x = "Irrigated area (fraction)", 
           y = expression(paste("Irrigation water withdrawal ", " ", "(", 10^9, m^3/year, "", ")"))) +
      theme_AP() +
      theme(strip.text.x = element_text(size = 6)) +
      ggtitle(label = names(full.isimip.split[i]), 
              subtitle = names(full.isimip.split[[i]][j]))
  }
}

gg


## ----session_information-------------------------------------------------------------------

# SESSION INFORMATION ------------------------------------------------------------

sessionInfo()

## Return the machine CPU
cat("Machine: "); print(get_cpu()$model_name)

## Return number of true cores
cat("Num cores: "); print(detectCores(logical = FALSE))

## Return number of threads
cat("Num threads: "); print(detectCores(logical = TRUE))

## Return the machine RAM
cat("RAM: "); print (get_ram()); cat("\n")

