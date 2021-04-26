
# Irrigated areas drive irrigation water withdrawals

Arnald Puy, Emanuele Borgonovo, Samuele Lo Piano, Simon Levin and Andrea Saltelli

This is the R code of the paper, whose abstract is the following:

*A sustainable management of global freshwater resources requires producing reliable estimates of the water demanded by irrigated agriculture. This has been attempted by FAO through country surveys and censuses, or through Global Models (GM), which compute irrigation water withdrawals with sub-models on crop types and calendars, evapotranspiration, irrigation efficiencies, weather data and irrigated areas, among others. Here we demonstrate that these strategies err on the side of excess complexity, as the values reported by FAO and outputted by GM are largely driven by irrigated areas only. Modeling irrigation water withdrawals as a function of irrigated areas yields almost the same results in a much cheaper and parsimonious way, while permitting the exploration of all model uncertainties. Our work offers a more robust and transparent approach to compute one of the most important indicators guiding our policies on water security worldwide.*

## Information

We provide the code in `.R`, `.rmd` and `.pdf` along with a detailed narration of the steps needed to reproduce our work and the software requirements. 

## Datasets included

Here you can find country-level estimates of irrigation water withdrawal produced as
a function of irrigated areas (`global_water_withdrawals.csv`). 

We also include all `.csv` datasets needed to replicate our work. The `.nc` files weight too much and should be downloaded from either ISIMIP (https://www.isimip.org) or from Huang et al. https://zenodo.org/record/1209296#.YIblVy2cb1J. All links below were functional as of 26 April 2021:

* `FAO_GMIA_historic.csv` https://mygeohub.org/publications/8/2.
* `table_4.csv` http://www.fao.org/3/bc824e/bc824e.pdf.
* `liu.csv` https://pubs.acs.org/doi/10.1021/acs.est.6b01065.
* `meier.csv` https://hess.copernicus.org/articles/22/1119/2018/.
* `rohwer_data.csv` https://www.pik-potsdam.de/en/output/publications/pikreports/.files/pr104.pdf.
* `australia_scheme.csv` https://www.icid.org/BMGuidelines.pdf.
* `solley.dt.csv` https://pubs.er.usgs.gov/publication/cir1200.
* `colorado_data.csv` https://pubs.usgs.gov/sir/2010/5002/pdf/SIR10-5002.pdf.
* `water_productivity.csv` https://data.worldbank.org/indicator/ER.GDP.FWTL.M3.KD.
