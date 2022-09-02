###########################################################
#   Author: Isabel De Ramos                               #
#   Date Created: 13 February 2022                        #
#   Function: POP GRADUATION / MORTALITY SINGLE YR AGE    #
###########################################################

### loading libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)

#################### POPULATION GRADUATION #################### 
#  OBTAINING SINGLE YR ESTIMATES BY METRO6, YEAR, AGE, SEX RACE

# LOADING IN MASTER DTA 
load("../1990_2019_nchs_mortality.rdata")

# CREATING HELPER FXNS
sort_age_group <- c(0, 1, 5, seq(10,85, by=5)) # SORTS AGE GROUP 

# PLCM GRADUATION FXN, RETURNS DATAFRAME 
pclm_grad <- function(.x, .y) {
  .x <- .x %>% mutate(age_5yr_group=as.numeric(age_5yr_group)) %>% 
    arrange(match(age_5yr_group, sort_age_group))
  results <- graduate(Value=.x$pop_denom, 
                      Age=.x$age_5yr_group, 
                      AgeInt=age2int(.x$age_5yr_group), 
                      method="pclm",
                      keep0=TRUE,
                      OAG=TRUE,
                      OAnew=110)
  single_yr_age  <- names2age(results)
  as_tibble(single_yr_age) %>% mutate(single_yr_age=value, pop_denom=results) %>% select(-value)
}

# TESTING FXN ON ONE FIPS, YEAR, SEX, RACE
# RE-GROUPING AND SUMMARIZING BY YEAR, AGE GROUP, SEX, RACE
test <- master_dta %>% select(-fips, -metro2, -census_region) %>% 
  mutate(age_5yr_group=as.numeric(age_5yr_group)) %>% 
  group_by(age_5yr_group, sex, race, year, metro6) %>% 
  arrange(match(age_5yr_group, sort_age_group)) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  filter(metro6%in%"1", year%in%"1990", sex%in%"1", race%in%"NHB") %>% 
  mutate(age_5yr_group=as.numeric(age_5yr_group)) %>% 
  pclm_grad()


# COMPARING PLOTS, PCLM GRADUATED FIRST
ggplot(test, aes(single_yr_age, pop_denom)) + 
  geom_line() + 
  scale_y_continuous(labels=comma)+
  ggtitle("graduation with pclm, 1990 male metro NHB")+
  labs(x="age",
       y="population") +
  isabel_theme

# COMPARING PLOTS, NO GRADUATION  
test2 <- master_dta %>% select(-fips, -metro2, -census_region) %>% 
  mutate(age_5yr_group=as.numeric(age_5yr_group)) %>% 
  group_by(age_5yr_group, sex, race, year, metro6) %>% 
  arrange(match(age_5yr_group, sort_age_group)) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  filter(metro6%in%"1", year%in%"1990", sex%in%"1", race%in%"NHB") %>% 
  mutate(age_5yr_group=as.numeric(age_5yr_group))

ggplot(test2, aes(age_5yr_group, pop_denom)) + 
  geom_line() + 
  scale_y_continuous(labels=comma)+
  ggtitle("no graduation, 1990 male metro NHB")+
  labs(x="age",
       y="population") +
  isabel_theme


# TAKING SINGLE YR ESTIMATES, REAGGREGATING TO 5 YR, AND CHECKING TO SEE IF IT MATCHES 5 YR POP ESTIMATES 
test %>% mutate(pop_denom_singleyr=pop_denom,
                age_5yr_group=case_when(
                  single_yr_age == 0 ~ '0',
                  single_yr_age < 5 ~ '1',
                  between(single_yr_age, 5,84)~as.character(floor(single_yr_age/5)*5),
                  TRUE~"85")) %>% 
  group_by(age_5yr_group) %>% 
  summarise(pop_denom_singleyr=sum(pop_denom_singleyr, na.rm=T)) %>% 
  arrange(match(age_5yr_group, sort_age_group)) %>% 
  full_join(test2 %>% mutate(age_5yr_group=as.character(age_5yr_group)), by="age_5yr_group") %>% 
  mutate(diff=pop_denom_singleyr-pop_denom)
# after checking diff column, agreed that differences were overall negligble b/w 
# single yr reaggregated to 5 yr estimates and original 5 yr estimates 


# TURNING PCLM.GRAD INTO APPLICABLE FUNCTION TO RACE-URBANICITY GROUPS 
# RE-GROUPING AND SUMMARIZING BY YEAR, AGE GROUP, SEX, RACE
pclm_results <- master_dta %>% select(-fips, -metro2, -census_region) %>% 
  mutate(age_5yr_group=as.numeric(age_5yr_group)) %>% 
  group_by(age_5yr_group, sex, race, year, metro6) %>% 
  arrange(match(age_5yr_group, sort_age_group)) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  group_by(year, sex, race, metro6) %>% 
  group_modify(~pclm_grad(.))


# PCLM_RESULTS IS NOW SINGLE YR POPULATION ESTIMATES BY METRO6, YEAR, RACE, SEX UP TO AGE 110 


#################### OBTAINING MORTALITY COUNTS BY SINGLE YR AGE #################### 
### Get mortality file names
# NOTE: individual level mortality files obtained from NCHS were cleaned and exported as .rdata's beforehand 
df_mort_files = tibble(file = list.files(path='../../NCHS Mortality Data/Clean/Raw', pattern = "mort[1|2]", full.names = TRUE)) %>% 
  mutate(year = str_sub(file,41,44) %>% parse_integer() ) %>% 
  arrange(year) 


# CREATING FUNCTIONS TO CLEAN AND PREPARE NCHS MORTALITY AND POPULATION DATA 
# (1) clean_nchs 
# (2) clean_nchs1


#~~~~ (1) clean_nchs ~~~~#
# clean_nchs extracts variables of interest from mortality data - year, fips, age, sex, and race
clean_nchs <- function(file) {
  #file <- df_mort_files %>% filter(year%in%"1991") %>% pull(file) 
  load(file)
  # after loading, object loaded into global environment is mortXXXX
  
  # NCHS DATASET - find counts per age group by race
  mort_dta <- mort_dta %>% 
    transmute(year=as.character(death_year),
              fips=res_fips_effective18,
              age=as.numeric(age_red),
              single_yr_age=ifelse(age>=110, 110, age),
              sex=male,
              race=hispanic_bridged) %>%
    group_by(fips, year, single_yr_age, sex, race) %>% 
    count(single_yr_age) %>%
    rename(count=n) %>% 
    arrange(fips, year, single_yr_age, sex, race) %>%
    ungroup()
  print(mort_dta) %>% 
    return()
}

#~~~~ (2) clean_nchs1 ~~~~#
clean_nchs1 <- function(nchs_dta_tmp) {
  #nchs_dta_tmp <- nchs
  
  # (1) filters out non-continental FIPS 
  nonUS_FIPS=c("^(02)","^(60)","^(66)","^(69)","^(72)","^(78)") 
  # Alaska, American Samoa, Guam, Northern Mariana Islands, Puerto Rico, Virgin Islands
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!grepl(paste(nonUS_FIPS, collapse="|"), fips))
  
  # (2) filters out foreign residents (fips="00000", "00999)
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!fips%in%c("00000", "00999", "13999"))
  # 13999 (county of < 100k or Georgia HIV death)
  # NOTE: For the years 1988 through 1991, if there were three or fewer deaths for a given Georgia county of residence 
  # (of deaths occurring in Georgia) with HIV infection (ICD codes *042-*044, 795.8) 
  # cited as a cause-of-death (underlying or non-underlying cause), these records were assigned a "missing" place of 
  # residence code (location code (FIPS code 13999).
  # These deaths do not appear in county death rates, but these deaths are included in the state and national death rates.
    
  # (3) filters out missing fips (checked in RAW data, missing fips were res_fips from US territories (PR, GUA, VI , etc.))
  nchs_dta_tmp <- nchs_dta_tmp %>% filter(!is.na(fips)==TRUE)
  
  # (4) fixes FIPS that were renamed over time 
  nchs_dta_tmp <- nchs_dta_tmp %>% 
    mutate(fips=ifelse(fips=="12025", "12086", fips)) %>% # fixing 12025
    mutate(fips=ifelse(fips=="51560", "51005", fips)) %>% # fixing 51560
    mutate(fips=ifelse(fips=="51780", "51083", fips)) %>% # fixing 51780
    # special case: 30113 
    # pop_denoms %>% filter(fips%in%"30031") %>% count(fips, wt=pop_denom) >> population of 1863635
    # pop_denoms %>% filter(fips%in%"30067") %>% count(fips, wt=pop_denom) >> population of 332532
    # since FIPS 30031 has largest population, 30113, 30067 and 30031 will all merge into 30031 
    mutate(fips=ifelse(fips%in%c("30113", "30067"), "30031", fips)) %>%
    arrange(year, single_yr_age, race) %>% 
    group_by(fips, year, single_yr_age, sex, race) %>%
    summarise(count=sum(count)) %>% 
    ungroup()
  print(nchs_dta_tmp) %>% 
    return()
}

# RUNNING CLEAN_NCHS AND CLEAN_NCHS1 ON MORT FILES 1990-2019
mort_files_tmp = df_mort_files %>% filter(between(year, 1990, 2019)) %>% pull(file)
nchs <- map_dfr(mort_files_tmp, ~clean_nchs(.x))
nchs <- clean_nchs1(nchs) 

# 7/11/22: R&R REVSION: DROP ALL NHAIAN, DROP HISP IN FOLLOWING STATES AND YEARS:
# Louisiana: 1990
# New Hampshire: 1990-1992
# Oklahoma: 1990-1996
hisp_FIPS_LA=c("^(22)")
hisp_FIPS_NH=c("^(33)")
hisp_FIPS_OK=c("^(40)")

nchs <- nchs %>% filter(!race%in%"NHAIAN")  # DROPPING NHAIAN ALL TOGETHER
nchs <- nchs %>% filter(!((year%in%"1990" & race%in%"H" & grepl(paste(hisp_FIPS_LA, collapse="|"), fips)) | # DROPPING HISP IN LA 1990
                            (year%in%c("1990", "1991", "1992") & race%in%"H" & grepl(paste(hisp_FIPS_NH, collapse="|"), fips)) | # DROPPING HISP IN NH 1990-1992
                            (year%in%c("1990", "1991", "1992", "1993", "1994", "1995", "1996") & race%in%"H" & grepl(paste(hisp_FIPS_OK, collapse="|"), fips)) )) # DROPPING HISP IN OK 1990-1996


# ADDING URBANICITY CLASSIFICATION TO NCHS
# create function that pads all fips (adds back leading zeroes)
format_5dig_fips <- function(fp) {
  ifelse(str_length(fp)==4,
         str_pad(fp, width=5, side="left", pad="0"),
         fp)
}
format_5digfips_v <- Vectorize(format_5dig_fips)
# reading in 2013 NCHS urban-rural crosswalk 
xwalk_fips_urbanrural <- read.csv('../Crosswalks/US/Clean/fips_urbanrural_xwalk.csv',  
                                  sep = ',', header =  T, colClasses = "character") %>% as_tibble() %>%
  transmute(fips=format_5digfips_v(fips),
            `metro6`=`X2013_code`)
# merging urban-rural codes to data
nchs_tmp <- nchs %>% full_join(xwalk_fips_urbanrural, by="fips") %>% 
  # creating `metro_class` column to delineate metropolitan or nonmetropolitan
  mutate(metro2=ifelse(metro6%in%c("1", "2", "3", "4"), "1", "0"))

# NCHS SHOULD NOW HAVE SINGLE YR MORTALITY COUNTS UP TO AGE 110 BY METRO2, YEAR, SEX, RACE 



#################### COMBINING SINGLE YR MORTALITY DTA AND POP ESTIMATES #################### 
# DROPPING FIPS AND RESUMMARIZING MORT COUNTS 
nchs_tmp <- nchs_tmp %>% select(-fips) %>% 
  group_by(year, sex, race, single_yr_age, metro6) %>% 
  summarise(count=sum(count, na.rm=T)) %>% 
  ungroup()

# JOINING MORT (NCHS) TO POP ESTIMATES (PCLM_RESULTS)
single_yr_mortdta <- nchs_tmp %>% full_join(pclm_results, by=c("year", "single_yr_age", "sex", "race", "metro6")) %>% 
  arrange(year, single_yr_age, sex, race) %>% 
  filter(!is.na(year)==TRUE)
# CREATING TEMPLATE TO SEE FULL COMBO OF AGE GROUPS / SEXES / RACES
template <- expand_grid(year=unique(single_yr_mortdta$year),
                        single_yr_age=unique(single_yr_mortdta$single_yr_age),
                        sex=unique(single_yr_mortdta$sex),
                        race=unique(single_yr_mortdta$race),
                        metro6=unique(single_yr_mortdta$metro6)) %>% 
  arrange(year, single_yr_age, sex, race, metro6)
  # REPLACING NA DEATH COUNTS W/ 0
single_yr_mortdta <- template %>% full_join(single_yr_mortdta, by=c("year", "single_yr_age", "sex", "race", "metro6")) %>% 
    mutate(count=replace_na(count,0)) %>% 
    arrange(year, single_yr_age, sex, race, metro6)

single_yr_mortdta <- single_yr_mortdta %>% 
  mutate(single_yr_age=ifelse(single_yr_age>=85, 85, single_yr_age)) %>% 
           group_by(single_yr_age, race, metro6, year, sex) %>% 
           summarise(pop_denom=sum(pop_denom), 
                     count=sum(count))


#################### SAVING SINGLE YR DTA FILE #################### 
save(single_yr_mortdta, file='single_yr_mortdta.rdata')
  




#################### LIFESPAN VARIATION SENSITIVITY ANALYSIS #################### 

# WORKING OUT FUNCTIONS THAT CALCULATE THE FOLLOWING LIFESPAN VARIATION MEASURES
# VARIANCE 
# STANDARD DEVIATION UNCONDITIONAL(le_lv())
# COEFFICIENT OF VARIANCE
# LIFE DISPARITY
# GINI (HAMADA METHOD)
# IQR
# AID 

## FOR THIS SENSITIVITY ANALYSIS: USED GRADUATED (SINGLE YR DTA) FOR ALL LIFESPAN VARIATION MEASURES 

# LOADING IN GRADUATED (SINGLE YR DTA) (CALCLATED IN SINGLE_YR_DTA.R) 
load("single_yr_mortdta.rdata")
single_yr_mortdta <- single_yr_mortdta %>%   
  mutate(year5=case_when(
    year%in%c(2015:2019) ~ as.character('2015-2019'),
    year%in%c(2010:2014) ~ as.character('2010-2014'),
    year%in%c(2005:2009) ~ as.character('2005-2009'),
    year%in%c(2000:2004) ~ as.character('2000-2004'),
    year%in%c(1995:1999) ~ as.character('1995-1999'),
    year%in%c(1990:1994) ~ as.character('1990-1994'))) %>% 
  mutate(metro2=ifelse(metro6%in%c("1", "2", "3", "4"), "1", "0"),
         metro2=factor(metro2, levels=c(1,0), labels=c("Metropolitan", "Non-metro")),
         gender=factor(sex, levels=c(1,0), labels=c("Men", "Women"))) %>% select(-sex, -year)

# TRYING OUT GINI ON ONE METRO, RACE, YEAR CATEGORY
.x<- single_yr_mortdta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age))

lt <- lt_single_mx(nMx=.x$mx, Age=.x$single_yr_age, sex="b", OAG=TRUE)

#~~~~ STANDARD DEVIATION UNCONDITIONAL AND CONDITIONAL (le_lv()) ~~~~#
le_lv<-function(.x, .y, age_num){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994')
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate lifespan variation 
  lt <- lt %>% mutate(xbar=NA_real_,
                      noname=NA_real_,
                      v=NA_real_,
                      sd=NA_real_)
  lt <- lt %>% mutate(ax=ifelse(x==110, 1.29, ax),
                      xbar=x+ax,
                      noname=case_when(
                        x==0 ~ dx/lx*(xbar-ex)^2,
                        x!=0 ~ {
                          lx_forlv <- lt %>% filter(x==0) %>% pull(lx)
                          ex_forlv <- lt %>% filter(x==0) %>% pull(ex)
                          dx/lx_forlv*(xbar-ex_forlv)^2}),
                      v=rev(cumsum(rev(noname))),
                      sd=sqrt(v))
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting LV 
  lv <- lt %>% filter(x==age_num) %>% pull(sd)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv=lv)
}


# OBTAINING STANDARD DEVIATION 
e_sd_0 <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, race, single_yr_age) %>% 
  group_by(metro2, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(., age_num=0)) 

# CREATING OVERALL METRO TO NONMETRO DISPARITY 
e_sd_0_overall <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange( metro2, single_yr_age) %>% 
  group_by(metro2) %>%
  group_modify(~le_lv(., age_num=0)) 

# MERGING OVERALL WITH RACE RESULTS
e_sd_0 <- e_sd_0_overall %>% mutate(race=rep("Overall")) %>% 
  select(metro2, race, le, lci, uci, lv) %>% 
  bind_rows(., e_sd_0) %>% 
  arrange(metro2, race)



#~~~~ COEFFICIENT OF VARIANCE ~~~~#
coeff_var<-function(.x, .y, age_num){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994')
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate lifespan variation 
  lt <- lt %>% mutate(xbar=NA_real_,
                      noname=NA_real_,
                      v=NA_real_,
                      sd=NA_real_)
  lt <- lt %>% mutate(ax=ifelse(x==110, 1.29, ax),
                      xbar=x+ax,
                      noname=case_when(
                        x==0 ~ dx/lx*(xbar-ex)^2,
                        x!=0 ~ {
                          lx_forlv <- lt %>% filter(x==0) %>% pull(lx)
                          ex_forlv <- lt %>% filter(x==0) %>% pull(ex)
                          dx/lx_forlv*(xbar-ex_forlv)^2}),
                      v=rev(cumsum(rev(noname))),
                      sd=sqrt(v))
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting LV 
  sd <- lt %>% filter(x==age_num) %>% pull(sd)
  CoV <- sd/le
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv=CoV)
}


# OBTAINING STANDARD DEVIATION 
CoV_0 <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, race, single_yr_age) %>% 
  group_by(metro2, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~coeff_var(., age_num=0)) 

# CREATING OVERALL METRO TO NONMETRO DISPARITY 
CoV_0_overall <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange( metro2, single_yr_age) %>% 
  group_by(metro2) %>%
  group_modify(~coeff_var(., age_num=0)) 

# MERGING OVERALL WITH RACE RESULTS
CoV_0 <- CoV_0_overall %>% mutate(race=rep("Overall")) %>% 
  select(metro2, race, le, lci, uci, lv) %>% 
  bind_rows(., CoV_0) %>% 
  arrange(metro2, race)




#~~~~ LIFE DISPARITY (E-DAGGER) ~~~~#
e_dagger <-function(.x, .y, age_num){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994')
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate gini
  lt <- lt %>% mutate(lost_length=NA_real_,
                      e_dagger=NA_real_)
  ex_85 <- lt %>% filter(x==85) %>% pull(ex)
  lx_85 <- lt %>% filter(x==85) %>% pull(lx)
  lt <- lt %>% mutate(ax=ifelse(x==85, 1.29, ax),
                      lost_length=case_when(
                        x==85 ~ 0,
                        x!=85 ~ lead(ex)+1-ax),
                      e_dagger=(rev(cumsum(dx*lost_length)+0.5*ex_85*lx_85))/lx)
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting e_dagger
  e_dagger <- lt %>% filter(x==age_num) %>% pull(e_dagger)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv=e_dagger)
}

# TESTING E_DAGGER
e_dag_results <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, race, single_yr_age) %>% 
  group_by(metro2, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~e_dagger(., age_num=0))

# CREATING OVERALL METRO TO NONMETRO DISPARITY 
e_dag_overall <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange( metro2, single_yr_age) %>% 
  group_by(metro2) %>%
  group_modify(~e_dagger(., age_num=0)) 

# MERGING OVERALL WITH RACE RESULTS
e_dag_results <- e_dag_overall %>% mutate(race=rep("Overall")) %>% 
  select(metro2, race, le, lci, uci, lv) %>% 
  bind_rows(., e_dag_results) %>% 
  arrange(metro2, race)


## NOTE: GINI AND IQR METHODS REQUIRE UNABRIDGED LIFE TABLES 
#~~~~  GINI using HANADA METHOD ~~~~#
gini <-function(.x, .y, age_num){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994')
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate gini
  # Tx_0 <- lt %>% filter(x==0) %>% pull(Tx)
  # lt <- lt %>% mutate(dx_l0=NA_real_,
  #                     dx_t0=NA_real_,
  #                     Fx=NA_real_,
  #                     FFx=NA_real_,
  #                     integral=NA_real_,
  #                     Gx=NA_real_, 
  #                     AID=NA_real_)
  # lt <- lt %>% mutate(dx_l0=dx/100000,
  #                     dx_t0=(dx*(x+ax))/Tx_0,
  #                     Fx=cumsum(dx_l0),
  #                     FFx=cumsum(dx_t0),
  #                     integral=(lead(Fx)-Fx)*(lead(FFx)+FFx),
  #                     integral=replace_na(integral,0),
  #                     Gx=1-sum(integral),
  #                     AID=ex*Gx)
  # adding variables needed to calculate gini (HANADA METHOD)
  lt <- lt %>% mutate(lx_sq=NA_real_,
                      integral=NA_real_,
                      Gx=NA_real_,
                      AID=NA_real_,)
  lt <- lt %>% mutate(lx_sq=lx^2/100000^2,
                      integral=lead(lx_sq)*(lead(x)-x)+ax*(lx_sq-lead(lx_sq))*(lead(x)-x),
                      integral=replace_na(integral,0),
                      Gx=1-rev(cumsum(integral))/ex*lx_sq,
                      AID=ex*Gx)
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting gini
  gini <- lt %>% filter(x==age_num) %>% pull(Gx)
  # extracting AID 
  aid <- lt %>% filter(x==age_num) %>% pull(AID)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             gini=gini, 
             aid=aid)
}



# TESTING GINI
gini_aid_results <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, race, single_yr_age) %>% 
  group_by(metro2, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~gini(., age_num=0)) 

# CREATING OVERALL METRO TO NONMETRO DISPARITY 
gini_aid_overall <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, single_yr_age) %>% 
  group_by(metro2) %>%
  group_modify(~gini(., age_num=0)) 

# MERGING OVERALL WITH RACE RESULTS
gini_aid_results <- gini_aid_overall %>% mutate(race=rep("Overall")) %>% 
  select(metro2, race, le, lci, uci, gini, aid) %>% 
  bind_rows(., gini_aid_results) %>% 
  arrange(metro2, race)



#~~~~ IQR ~~~~#
#QUARANTINED - IQR WITH ABRIDGED LIFE TABLES AND LAST AGE GROUP OF 85+ GENERATES IMPRECISE IQR
iqr <-function(.x, .y, age_num){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994')
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$single_yr_age, Nx=.x$pop_denom, Dx=.x$count, sex="B", ax=NULL, which.x=age_num, ns=1000, level=0.95)
  intFUN <- function(x, y, n, n.grid=10000){ 
    fit <- spline(x=x, y=y, n=n.grid)
    xi <- fit$x
    yi <- fit$y
    xn <- xi[which.min(abs(yi-n))] 
    return(xn) 
  }
  age25 <- intFUN(x=lt$x, y=lt$lx, n=75000) 
  age75 <- intFUN(x=lt$x, y=lt$lx, n=25000) 
  # extracting LE
  le <- lt %>% filter(x==age_num) %>% pull(ex)
  # extracting 95% CI
  ci <- ci_dta %>% pluck("CIex") %>% unname() 
  # extracting IQR
  iqr <- age75-age25
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             iqr=iqr,
             iqr_25=age25,
             iqr_75=age75)
}

# TESTING IQR
iqr_results <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange(metro2, race, single_yr_age) %>% 
  group_by(metro2, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~iqr(., age_num=0)) 

# CREATING OVERALL METRO TO NONMETRO DISPARITY 
iqr_overall <- single_yr_mortdta %>% filter(year5=='2015-2019') %>% 
  group_by(single_yr_age, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  mutate(mx=count/pop_denom, 
         single_yr_age=as.numeric(single_yr_age)) %>% 
  arrange( metro2, single_yr_age) %>% 
  group_by(metro2) %>%
  group_modify(~iqr(., age_num=0)) 

# MERGING OVERALL WITH RACE RESULTS
iqr_results <- iqr_overall %>% mutate(race=rep("Overall")) %>% 
  select(metro2, race, le, lci, uci, iqr, iqr_25, iqr_75) %>% 
  bind_rows(., iqr_results) %>% 
  arrange(metro2, race)


#### SENSITIVITY ANALYSIS FOR LV MEASURES ####
# VARIANCE (AT BIRTH)
# LIFESPAN VARIATION S(0), S(10), S(35), S(65),
# LIFE DISPARITY(E-DAGGER)
# IQR

# finding variance
variance <- e_sd_0 %>% mutate(lv=lv^2,
                              lv_measure=rep("var"))


# COMBINING LV MEASURES INTO ONE DATAFRAME
lv_dta <- variance %>% bind_rows(e_sd_0 %>% mutate(lv_measure=rep("sd"))) %>% 
  bind_rows(., CoV_0 %>% mutate(lv_measure=rep("CoV"))) %>% 
  bind_rows(., e_dag_results %>% mutate(lv_measure=rep("e_dagger"))) %>%
  bind_rows(., iqr_results %>% mutate(lv=iqr, lv_measure=rep("iqr"))) %>% 
  bind_rows(., gini_aid_results %>% select(-lci, -uci) %>% gather(key="lv_measure", value="lv", -metro2, -race, -le)) %>%
  select(-lci, -uci, -iqr_25, -iqr_75, -iqr) 
lv_dta_cor <- lv_dta %>% select(-le) %>% spread(lv_measure, lv)

# FIRST, DO ALL LV MEASURES FOLLOW NORMAL DISTRIBUTION 
# USE SHAPIRO-WILK NORMALITY TEST H0: NORMAL DIST, HA: NOT NORMAL DIST
shapiro.test(lv_dta_cor$e_dagger)
shapiro.test(lv_dta_cor$sd)
shapiro.test(lv_dta_cor$CoV)
shapiro.test(lv_dta_cor$var)
shapiro.test(lv_dta_cor$iqr)
shapiro.test(lv_dta_cor$gini)
shapiro.test(lv_dta_cor$aid)
# P-VALUES FOR ALL LV MEASURES ARE > 0.05, IMPLYING THAT WE CAN ASSUME NORMALITY 


# NEXT, FIND PEARSON CORRELATION COEFFICIENTS FOR LV MEASURES 
# install.packages("corrplot")
library(corrplot)
# install.packages("GGally")
library(GGally)

# CORR MATRIX
lv_dta_cor <- lv_dta_cor %>% ungroup() %>% select(-metro2, -race) 
lv_dta_cor %>% cor(., method = "pearson") 


# CORR MATRIX W/ VISUALIZATIONS 
sens_corrmatrix <- ggpairs(lv_dta_cor[, c(6,7,2,4,1,3,5)], 
        columnLabels = c("sd", "variance", "CoV", "gini", "AID", "e-dagger", "IQR"),
        upper = list(continuous = wrap('cor', size=4))) + 
  isabel_theme + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
ggsave("../Tables & Figures/sens_corrmatrix.pdf", sens_corrmatrix, width=15, height=8)


