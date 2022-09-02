## Exploring inequalities in life expectancy and lifespan variation by race/ethnicity and urbanicity in the United States: 1990-2019.

Under Review at Social Science and Medicine - Population Health

**NCHS MORTALITY DATA folder contains:**

Originally, the NCHS vital resgitration data. However, since this data cannot be shared publicly, we just put a link to the data request. 

**POPULATION DATA folder contains:**
- Bridged-race post-censal population estimates obtained from Census Bureau (https://www.cdc.gov/nchs/nvss/bridged_race/data_documentation.htm)
- _LE_pop_denoms.rdata_: fips, year (1990-2019), age 5 yr group (0,1,5,10...), gender (M, W), race (NHW, NHB, H, NHAPI), population denominators.

**ANALYSIS folder contains:**
- LE_data_prep.R
  - uses _LE_nchs.rdata_ and _LE_pop_denoms.rdata_ to generate _1990_2019_nchs_mortality.rdata_, master data for overall analysis that contains fips, year (1990-2019), age 5 yr group (0,1,5,10...), sex (M, F), race (NHW, NHB, H, NHAPI), death counts, population denominators, metropolitan status (metropolitan vs nonmetropolitan, disaggregated metropolitan status (large central metro, large fringe metro, medium metro, small metro, micro, noncore), and census region (northeast, midwest, south, west). Of note, _1990_2019_nchs_mortality.rdata_ utilizes 2013 NCHS Urban-Rural Classification Scheme. 
  - also generates _1990_2019_90class.rdata_, which contains same variables as _1990_2019_nchs_mortality.rdata_ except utilizes 1990 NCHS Urban-Rural Classification Scheme. 
  - see LE Data Prep Documentation.docx for further detail.
  - Please note that we cannot share the rdata files publicly since they are summaries of the NCHS data (see above) with cells with < 10 counts. 
- LE_calc_analysis.R 
  - main analyses code for life expectancy and lifespan variation calculations by race/ethnicity, gender, and urbanicity.
  - see LE Calc Analysis Documentation.docx for further detail.
- LifeTableFUN.R
  - Helper code/function to create lifetables in LE_calc_analysis.R
- Crosswalks
  - contains crosswalks from FIPS to 2013/2006/1990 NCHS Urban-Rural Classification.
  - contains crosswalks from FIPS to Census Region Classification.
- Tables & Figures
  - contains all ggplot figures generated from LE_calc_analysis.R
- Sensitivity Analysis 
  - contains code and data needed for analysis comparing absolute and relative measures of lifespan variation.
