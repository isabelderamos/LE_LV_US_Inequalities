###################################################
#   Author: Isabel De Ramos                       #
#   Date Created: 13 October 2021                 #
#   Function: Life Expectancy Calculations        #
###################################################

### loading libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)
library(grid)
library(gridExtra)
library(scales)
library(multcomp)
library(RColorBrewer)
library(scales)
library(ggpubr)
#install.packages("devtools")
#library(devtools)
#install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
#install.packages("remotes")
#remotes::install_github("timriffe/DemoTools")
#install.packages("DemoTools")
library(DemoTools)
library(RColorBrewer)
#remotes::install_github("coolbutuseless/ggpattern")
#install.packages("ggpattern")
library(ggpattern)
library(ggrepel)

#################### PREP WORK #################### 
# (1) ggplot helpers
# (2) Camarda's functions for lifetables, life expectancy 
# (3) loading in .rdata
# (4) create 5-year pooled periods and factors
# (5) create le_lv function

#~~~~ (1) ggplot helpers ~~~~#
### prep for generating plots ###
select<-dplyr::select
fontsize<-16
isabel_theme <-theme_bw()+
  theme(axis.text=element_text(color="black", size=fontsize),
        axis.title=element_text(color="black", size=fontsize, face="bold"),
        plot.title=element_text(color="black", size=fontsize, face="bold", hjust = 0.5),
        strip.background=element_blank(),
        strip.text=element_text(color="black", size=fontsize, face="bold"),
        legend.text=element_text(color="black", size=fontsize),
        legend.title=element_text(color="black", size=fontsize, face="bold"))
cols<-c(brewer_pal(type="qual", palette=2)(8),"blue", brewer_pal(type="qual", palette=2)(8), "blue", "red")
shapes<-c(rep(21, times=9), rep(22, times=9), 22)
# creating specific colors for racial groups
race_colors <- tibble(race=c("Overall", "NHW", "NHB", "NHAPI", "H")) %>% 
  mutate(color=c("Black", "#66A61E", "#E7298A", "#7570B3", "#1B9E77"))
race_cols <- race_colors$color
# creating specific colors for metro groups
metro_colors <- tibble(metro=c("Metropolitan", "Non-metro")) %>% 
  mutate(color=c("#366bcf", "#cf3636"))
metro_cols <- metro_colors$color
# creating sort order for races 
sort_order <- c("NHW", "NHB", "NHAPI", "H")
sort_order1 <- c("Overall", "NHW", "NHB", "NHAPI", "H")

#~~~~ (2) loading in Camarda's functions for lifetables, life expectancy ~~~~#
# see https://sites.google.com/site/carlogiovannicamarda/r-stuff/life-expectancy-confidence-interval 
## lifetable() calculates life expectancy
## CIex() calculates 95% CI 
source('LifeTableFUN.R')


#~~~~ (3) loading in master datafile (see LE_data_prep.R) ~~~~#
load('1990_2019_nchs_mortality.rdata')


#~~~~ (4) create year5 variable that pools 5-year periods from 2000-2019 ~~~~#
dta <- master_dta %>% select(-fips) %>% 
  mutate(gender=as.numeric(sex),
         age_5yr_group=as.numeric(age_5yr_group),
         year5=case_when(
          year%in%c(2015:2019) ~ as.character('2015-2019'),
          year%in%c(2010:2014) ~ as.character('2010-2014'),
          year%in%c(2005:2009) ~ as.character('2005-2009'),
          year%in%c(2000:2004) ~ as.character('2000-2004'),
          year%in%c(1995:1999) ~ as.character('1995-1999'),
          year%in%c(1990:1994) ~ as.character('1990-1994'))) %>% 
  # factoring metro2, metro6, gender, census_region to help in ggplot later
  mutate(metro2=factor(metro2, levels=c(1,0),
                       labels=c("Metropolitan", "Non-metro")),
         metro6=factor(metro6, levels=c(1:6),
                       labels=c("Large central", "Large fringe", "Medium", "Small", "Micro", "Noncore")),
         gender=factor(gender, levels=c(1,0),
                       labels=c("Men", "Women")),
         census_region=factor(census_region, levels=c(1:4),
                              labels=c("Northeast", "Midwest", "South", "West"))) %>% select(-sex)


#~~~~ (5) create le_lv function that calculates life expectancy and lifespan variation  ~~~~#
le_lv<-function(.x, .y, age_num, sex="B"){
  ## .x -> the data with one race and metro
  # example: .x<- dta %>% filter(metro2=="Metropolitan", race=="NHAIAN", year5%in%'1990-1994', gender=="Men")
  ## .y -> a dataframe wtih whichever metro or race category you are looking at
  # using lifetable() function from Camarda to create lifetable 
  lt <- lifetable(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex=sex, ax=NULL) 
  # using CIex function from Camarda to obtain confidence intervals 
  ci_dta <- CIex(x=.x$age_5yr_group, Nx=.x$pop_denom, Dx=.x$count, sex=sex, ax=NULL, which.x=age_num, ns=1000, level=0.95)
  # adding variables needed to calculate lifespan variation 
  lt <- lt %>% mutate(xbar=NA_real_,
                      noname=NA_real_,
                      v=NA_real_,
                      sd=NA_real_)
  lt <- lt %>% mutate(xbar=x+ax,
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
  lv <- sd/le # this is CoV = coefficient of variance (standard deviation / LE)
  # building df that summarizes results
  data.frame(le=le,
             lci=ci[1],
             uci=ci[2],
             lv=lv,
             sd=sd)
}

# ############ BASELINE ANALYSIS: 2015-2019, GENDER POOLED, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# baseline analysis examines disparities in life expectancy and lifespan variation by race/ethnicity
# and by urbanicity (2 categories) within the 5-year pooled period of 2015-2019
results_baseline <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, metro2, race) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  # creating age-specific death rates variable, Mx
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(metro2, race, age_5yr_group) %>%
  group_by(metro2, race) %>%
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(., age_num=0))

# CREATING OVERALL METRO TO NONMETRO DISPARITY
results_baseline_overall <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, metro2) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange( metro2, age_5yr_group) %>%
  group_by(metro2) %>%
  group_modify(~le_lv(., age_num=0))

# MERGING OVERALL WITH RACE RESULTS
results_baseline <- results_baseline_overall %>% mutate(race=rep("Overall")) %>%
  select(metro2, race, le, lci, uci, lv) %>%
  bind_rows(., results_baseline) %>%
  arrange(metro2, race)

# CREATING OVERALL OVERALL RACE DISPARITY
results_overall <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, race) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  # creating age-specific death rates variable, Mx
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(race, age_5yr_group) %>%
  group_by(race) %>%
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(., age_num=0))

# CREATING OVERALL OVERALL DISPARITY (basically LE for entire US pop)
results_overall_overall <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(age_5yr_group) %>%
  le_lv(., age_num=0)

results_overall <- results_overall_overall %>% mutate(race=rep("Overall")) %>%
  select(race, le, lci, uci, lv) %>%
  bind_rows(., results_overall) %>%
  arrange(race)




############ NEW BASELINE ANALYSIS: 2015-2019, MALE V FEMALE, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# analysis 2 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (2 categories), and by gender (men v women) within the 5-year pooled period of 2015-2019
# MEN BY RACE/URBANICITY LE/LV
results_gender <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, race, gender, age_5yr_group) %>% 
  group_by(race, metro2, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, metro2, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, race, gender, age_5yr_group) %>% 
      group_by(race, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%  
      group_by(age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%  
      group_by(age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )


############ BASELINE 0.5 ANALYSIS: 2015-2019, GENDER POOLED, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# life expectancy is e(10), e(35), e(65)
# lifespan variation is CoV(10), CoV(35), CoV(65)

# MEN FIRST 
BA_men_cond <- results_gender %>% filter(gender=="Men") %>% mutate(le_measure=rep("e0"), lv_measure=rep("CoV0")) %>% 
  bind_rows(
### LIFE EXP AND CoV 10 
dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, gender, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>%  
  filter(gender=="Men") %>% 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, gender, race, age_5yr_group) %>% 
  group_by(metro2, gender, race) %>% 
  # group_modify generates a data frame of outputs
  group_modify(~le_lv(., age_num=10, sex="M")) %>% mutate(le_measure=rep("e10"), lv_measure=rep("CoV10"))
) %>% 
  bind_rows(
# LIFE EXP AND CoV 10 OVERALL
dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, gender, metro2) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, gender, age_5yr_group) %>% 
  group_by(metro2, gender) %>%
  group_modify(~le_lv(., age_num=10, sex="M")) %>% mutate(race=rep("Overall"), 
                                                          le_measure=rep("e10"), lv_measure=rep("CoV10"))
) %>% 
  bind_rows(
    ### LIFE EXP AND CoV 35 
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2, race) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>%  
      filter(gender=="Men") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, race, age_5yr_group) %>% 
      group_by(metro2, gender, race) %>% 
      # group_modify generates a data frame of outputs
      group_modify(~le_lv(., age_num=35, sex="M")) %>% mutate(le_measure=rep("e35"), lv_measure=rep("CoV35"))
  ) %>% 
  bind_rows(
    # LIFE EXP AND CoV 35 OVERALL
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>%
      group_modify(~le_lv(., age_num=35, sex="M")) %>% mutate(race=rep("Overall"), 
                                                              le_measure=rep("e35"), lv_measure=rep("CoV35"))
  ) %>% 
  bind_rows(
    ### LIFE EXP AND CoV 65 
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2, race) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>%  
      filter(gender=="Men") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, race, age_5yr_group) %>% 
      group_by(metro2, gender, race) %>% 
      # group_modify generates a data frame of outputs
      group_modify(~le_lv(., age_num=65, sex="M")) %>% mutate(le_measure=rep("e65"), lv_measure=rep("CoV65"))
  ) %>% 
  bind_rows(
    # LIFE EXP AND CoV 65 OVERALL
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>%
      group_modify(~le_lv(., age_num=65, sex="M")) %>% mutate(race=rep("Overall"), 
                                                              le_measure=rep("e65"), lv_measure=rep("CoV65"))
  ) 


# NOW WOMEN 
BA_women_cond <- results_gender %>% filter(gender=="Women") %>% mutate(le_measure=rep("e0"), lv_measure=rep("CoV0")) %>% 
  bind_rows(
    ### LIFE EXP AND CoV 10 
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2, race) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>%  
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, race, age_5yr_group) %>% 
      group_by(metro2, gender, race) %>% 
      # group_modify generates a data frame of outputs
      group_modify(~le_lv(., age_num=10, sex="W")) %>% mutate(le_measure=rep("e10"), lv_measure=rep("CoV10"))
  ) %>% 
  bind_rows(
    # LIFE EXP AND CoV 10 OVERALL
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>%
      group_modify(~le_lv(., age_num=10, sex="W")) %>% mutate(race=rep("Overall"), 
                                                              le_measure=rep("e10"), lv_measure=rep("CoV10"))
  ) %>% 
  bind_rows(
    ### LIFE EXP AND CoV 35 
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2, race) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>%  
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, race, age_5yr_group) %>% 
      group_by(metro2, gender, race) %>% 
      # group_modify generates a data frame of outputs
      group_modify(~le_lv(., age_num=35, sex="W")) %>% mutate(le_measure=rep("e35"), lv_measure=rep("CoV35"))
  ) %>% 
  bind_rows(
    # LIFE EXP AND CoV 35 OVERALL
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>%
      group_modify(~le_lv(., age_num=35, sex="W")) %>% mutate(race=rep("Overall"), 
                                                              le_measure=rep("e35"), lv_measure=rep("CoV35"))
  ) %>% 
  bind_rows(
    ### LIFE EXP AND CoV 65 
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2, race) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>%  
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, race, age_5yr_group) %>% 
      group_by(metro2, gender, race) %>% 
      # group_modify generates a data frame of outputs
      group_modify(~le_lv(., age_num=65, sex="W")) %>% mutate(le_measure=rep("e65"), lv_measure=rep("CoV65"))
  ) %>% 
  bind_rows(
    # LIFE EXP AND CoV 65 OVERALL
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, gender, metro2) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>%
      group_modify(~le_lv(., age_num=65, sex="W")) %>% mutate(race=rep("Overall"), 
                                                              le_measure=rep("e65"), lv_measure=rep("CoV65"))
  )



############ ANALYSIS 1 LONGITUDINAL: 2000-2019, GENDER POOLED, GEOGRAPHY POOLED, 2 METRO CATEGORIES ######
# analysis 1 examines disparities in life expectancy and lifespan variation by race/ethnicity 
# and by urbanicity (2 categories) longitudinally across 5-year pooled period from 1990-2019
results_longitudinal <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% group_by(year5, age_5yr_group, metro2, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(year5, metro2, race, gender, age_5yr_group) %>% 
  group_by(year5, race, metro2, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, race, gender, age_5yr_group) %>% 
      group_by(year5, race, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, gender, age_5yr_group) %>% 
      group_by(year5, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, gender, age_5yr_group) %>% 
      group_by(year5, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )



############ ANALYSIS 3 URBANICITY: 2015-2019, GENDER POOLED, GEOGRAPHY POOLED, 6 METRO CATEGORIES ############ 
# analysis 3 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (6 categories) within the 5-year pooled period of 2015-2019
results_urbanicity <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, metro6, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro6, race, gender, age_5yr_group) %>% 
  group_by(race, metro6, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, race, gender, age_5yr_group) %>% 
      group_by(race, metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, gender, age_5yr_group) %>% 
      group_by(metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, gender, age_5yr_group) %>% 
      group_by(metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )

  
############ ANALYSIS 4 GEOGRAPHY: 2015-2019, GENDER POOLED, 4 CENSUS REGIONS, 2 METRO CATEGORIES  ######
# analysis 4 examines disparities in life expectancy and lifespan variation by race/ethnicity, 
# by urbanicity (2 categories), and by census region (NE, MW, S, W) within the 5-year pooled period of 2015-2019
# results_geography <- 
#   # MEN BY RACE/URBANICITY LE/LV
#   dta %>% filter(year5=='2015-2019') %>%
#   group_by(age_5yr_group, metro2, race, gender, census_region) %>% 
#   summarise(count=sum(count, na.rm=T),
#             pop_denom=sum(pop_denom, na.rm=T)) %>% 
#   ungroup() %>% 
#   filter(gender=="Men") %>% 
#   # creating age-specific death rates variable, Mx 
#   mutate(mx=count/pop_denom, 
#          age_5yr_group=as.numeric(age_5yr_group)) %>% 
#   arrange(metro2, race, gender, census_region, age_5yr_group) %>% 
#   group_by(race, metro2, gender, census_region) %>% 
#   group_modify(~le_lv(., age_num=0, sex="M")) %>% 
#   # WOMEN BY RACE/URBANICITY LE/LV
#   bind_rows(
#     dta %>% filter(year5=='2015-2019') %>%
#       group_by(age_5yr_group, metro2, race, gender, census_region) %>% 
#       summarise(count=sum(count, na.rm=T),
#                 pop_denom=sum(pop_denom, na.rm=T)) %>% 
#       ungroup() %>% 
#       filter(gender=="Women") %>% 
#       # creating age-specific death rates variable, Mx 
#       mutate(mx=count/pop_denom, 
#              age_5yr_group=as.numeric(age_5yr_group)) %>% 
#       arrange(metro2, race, gender, census_region, age_5yr_group) %>% 
#       group_by(race, metro2, gender, census_region) %>% 
#       group_modify(~le_lv(., age_num=0, sex="W"))
#   ) %>% 
#   # MEN BY URBANICITY OVERALL LE/LV
#   bind_rows(
#     dta %>% filter(year5=='2015-2019') %>%
#       group_by(age_5yr_group, metro2, gender, census_region) %>% 
#       summarise(count=sum(count, na.rm=T),
#                 pop_denom=sum(pop_denom, na.rm=T)) %>% 
#       ungroup() %>% 
#       filter(gender=="Men") %>%
#       mutate(mx=count/pop_denom, 
#              age_5yr_group=as.numeric(age_5yr_group)) %>% 
#       arrange(metro2, gender, census_region, age_5yr_group) %>% 
#       group_by(metro2, gender, census_region) %>% 
#       group_modify(~le_lv(., age_num=0, sex="B")) %>% 
#       mutate(race=rep("Overall"))
#   ) %>% 
#   # WOMEN BY URBANICITY OVERALL LE/LV
#   bind_rows(
#     dta %>% filter(year5=='2015-2019') %>%
#       group_by(age_5yr_group, metro2, gender, census_region) %>% 
#       summarise(count=sum(count, na.rm=T),
#                 pop_denom=sum(pop_denom, na.rm=T)) %>% 
#       ungroup() %>% 
#       filter(gender=="Women") %>%
#       mutate(mx=count/pop_denom, 
#              age_5yr_group=as.numeric(age_5yr_group)) %>% 
#       arrange(metro2, gender, census_region, age_5yr_group) %>% 
#       group_by(metro2, gender, census_region) %>% 
#       group_modify(~le_lv(., age_num=0, sex="W")) %>% 
#       mutate(race=rep("Overall"))
#   )




############ TABLES & FIGURES ############ 

#### BA_lelv_table : TABLE#### 
# Baseline Analysis table: le/lv and absolute/relative differences
# spreading LE and 95% CIs by metro2 
# le/lv and absolute/relative differences for each race-metro combination
BA_lelv_table <- results_baseline %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_baseline %>% select(race, metro2, le) %>% 
      spread(metro2, le) %>% 
      rename(metro_le=`Metropolitan`,
             nonmetro_le=`Non-metro`) %>% 
      full_join(
        results_baseline %>% select(race, metro2, lv) %>% 
          spread(metro2, lv) %>% 
          rename(metro_lv=`Metropolitan`,
                 nonmetro_lv=`Non-metro`)) %>% 
      mutate(abs_le=metro_le-nonmetro_le,
             rel_le=metro_le/nonmetro_le,
             abs_lv=metro_lv-nonmetro_lv,
             rel_lv=metro_lv/nonmetro_lv)) %>%
  ungroup() %>% 
  transmute(race, 
            metro_le_ci, 
            nonmetro_le_ci,
            abs_le=as.numeric(format(abs_le, digits=2, nsmall=2)), 
            rel_le=as.numeric(format(rel_le, digits=2, nsmall=2)), 
            metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1)), 
            nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1)), 
            abs_lv=as.numeric(format(abs_lv, digits=2, nsmall=2)), 
            rel_lv=as.numeric(format(rel_lv, digits=2, nsmall=2))) %>% 
  arrange(match(race, sort_order1))
write.csv(BA_lelv_table, "Tables & Figures/BA_lelv_table.csv", row.names=FALSE)


#### BA_relative2overall : TABLE#### 
# MEN FIRST

# finding overall race for men 
BA_rel2overall_men <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, gender, race) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  filter(gender=="Men") %>% 
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(race, gender, age_5yr_group) %>%
  group_by(race, gender) %>%
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # binding rows w/ overall-overall for men
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, gender) %>%
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>%
      ungroup() %>%
      filter(gender=="Men") %>% 
      mutate(mx=count/pop_denom,
             age_5yr_group=as.numeric(age_5yr_group)) %>%
      arrange(gender, age_5yr_group) %>%
      group_by(gender) %>%
      group_modify(~le_lv(., age_num=0, sex="M")) %>%
      mutate(race=rep("Overall"))
  ) %>%
  # finding abs and rel disparities for "overall men"
  select(-lci, -uci, -lv) %>% 
  arrange(match(race, sort_order1)) %>% 
  mutate(abs_disp=le-76.62171,
         rel_disp=le/76.62171,
         abs_disp=as.numeric(format(abs_disp, digits=2, nsmall=2)), 
         rel_disp=as.numeric(format(rel_disp, digits=2, nsmall=2))) %>% 
  # joining results_gender and finding abs and rel disparities for "overall race men"
    full_join(
      results_gender %>% filter(gender=="Men") %>% select(-lci, -uci, -lv) %>% 
        spread(metro2, le) %>% 
        arrange(match(race, sort_order1)) %>%
        mutate(abs_disp_met=Metropolitan-77.04266,
               rel_disp_met=Metropolitan/77.04266,
               abs_disp_met=as.numeric(format(abs_disp_met, digits=2, nsmall=2)), 
               rel_disp_met=as.numeric(format(rel_disp_met, digits=2, nsmall=2)),
               abs_disp_nonmet=`Non-metro`-74.47910,
               rel_disp_nonmet=`Non-metro`/74.47910,
               abs_disp_nonmet=as.numeric(format(abs_disp_nonmet, digits=2, nsmall=2)), 
               rel_disp_nonmet=as.numeric(format(rel_disp_nonmet, digits=2, nsmall=2)))) %>% 
    select(race, le, abs_disp, rel_disp,
           Metropolitan, abs_disp_met, rel_disp_met,
           `Non-metro`, abs_disp_nonmet, rel_disp_nonmet)
write.csv(BA_rel2overall_men, "Tables & Figures/BA_rel2overall_men.csv", row.names=FALSE)

# NOW WOMEN

BA_rel2overall_women <- dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, gender, race) %>%
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>%
  ungroup() %>%
  filter(gender=="Women") %>% 
  mutate(mx=count/pop_denom,
         age_5yr_group=as.numeric(age_5yr_group)) %>%
  arrange(race, gender, age_5yr_group) %>%
  group_by(race, gender) %>%
  group_modify(~le_lv(., age_num=0, sex="W")) %>% 
  # binding rows w/ overall-overall for women
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, gender) %>%
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>%
      ungroup() %>%
      filter(gender=="Women") %>% 
      mutate(mx=count/pop_denom,
             age_5yr_group=as.numeric(age_5yr_group)) %>%
      arrange(gender, age_5yr_group) %>%
      group_by(gender) %>%
      group_modify(~le_lv(., age_num=0, sex="W")) %>%
      mutate(race=rep("Overall"))
  ) %>% 
  # finding abs and rel disparities for "overall women"
  select(-lci, -uci, -lv) %>% 
  arrange(match(race, sort_order1)) %>% 
  mutate(abs_disp=le-81.64615,
         rel_disp=le/81.64615,
         abs_disp=as.numeric(format(abs_disp, digits=2, nsmall=2)), 
         rel_disp=as.numeric(format(rel_disp, digits=2, nsmall=2))) %>% 
  # joining results_gender and finding abs and rel disparities for "overall race women"
  full_join(
    results_gender %>% filter(gender=="Women") %>% select(-lci, -uci, -lv) %>% 
      spread(metro2, le) %>% 
      arrange(match(race, sort_order1)) %>%
      mutate(abs_disp_met=Metropolitan-82.03604,
             rel_disp_met=Metropolitan/82.03604,
             abs_disp_met=as.numeric(format(abs_disp_met, digits=2, nsmall=2)), 
             rel_disp_met=as.numeric(format(rel_disp_met, digits=2, nsmall=2)),
             abs_disp_nonmet=`Non-metro`-79.54948,
             rel_disp_nonmet=`Non-metro`/79.54948,
             abs_disp_nonmet=as.numeric(format(abs_disp_nonmet, digits=2, nsmall=2)), 
             rel_disp_nonmet=as.numeric(format(rel_disp_nonmet, digits=2, nsmall=2)))) %>% 
  select(race, le, abs_disp, rel_disp,
         Metropolitan, abs_disp_met, rel_disp_met,
         `Non-metro`, abs_disp_nonmet, rel_disp_nonmet)
write.csv(BA_rel2overall_women, "Tables & Figures/BA_rel2overall_women.csv", row.names=FALSE)


BA_relative2overall<- results_overall %>% select(-lci, -uci, -lv) %>% 
  arrange(match(race, sort_order1)) %>% 
  mutate(abs_disp=le-79.14423,
         rel_disp=le/79.14423,
         abs_disp=as.numeric(format(abs_disp, digits=2, nsmall=2)), 
         rel_disp=as.numeric(format(rel_disp, digits=2, nsmall=2))) %>% 
  full_join(
results_baseline %>% select(-lci, -uci, -lv) %>% 
  spread(metro2, le) %>% 
  arrange(match(race, sort_order1)) %>% 
  mutate(abs_disp_met=Metropolitan-79.56540,
         rel_disp_met=Metropolitan/79.56540,
         abs_disp_met=as.numeric(format(abs_disp_met, digits=2, nsmall=2)), 
         rel_disp_met=as.numeric(format(rel_disp_met, digits=2, nsmall=2)),
         abs_disp_nonmet=`Non-metro`-76.94276,
         rel_disp_nonmet=`Non-metro`/76.94276,
         abs_disp_nonmet=as.numeric(format(abs_disp_nonmet, digits=2, nsmall=2)), 
         rel_disp_nonmet=as.numeric(format(rel_disp_nonmet, digits=2, nsmall=2)))) %>% 
  select(race, le, abs_disp, rel_disp,
         Metropolitan, abs_disp_met, rel_disp_met,
         `Non-metro`, abs_disp_nonmet, rel_disp_nonmet)
write.csv(BA_relative2overall, "Tables & Figures/BA_relative2overall.csv", row.names=FALSE)



#### BA_le_bars : FIGURE #### 
# Baseline Analysis figure: le by race/urbanicity and gender
df <- results_baseline %>% ungroup() %>% mutate(gender=rep("Overall")) %>% select(race, metro2, gender, le, lci, uci) %>% 
  bind_rows(., results_gender %>% select(race, metro2, gender, le, lci, uci) %>% 
              mutate(gender=as.character(gender))) %>% 
  mutate(gender=case_when(
    gender=="Overall" ~ 0,
    gender=="Men" ~ 1,
    gender=="Women" ~ 2,
  )) %>% mutate(gender=factor(gender, levels=c(0,1,2), labels=c("Overall", "Men", "Women"))) %>% 
  arrange(race, metro2, gender)


BA_le_bars <- ggplot(df, aes(x=factor(race, level=sort_order1), y=le, fill=as.factor(gender))) + 
  geom_col(color="black", position=position_dodge(width=0.7), width=0.7) +
  geom_linerange(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.7))+
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-1) +
  scale_y_continuous(limits=c(df %>% ungroup() %>% select(le) %>% min,
                              df %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_brewer(palette="Set2")+
  labs(x="Race/Ethnicity",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  isabel_theme + 
  # emphasizes only y=lines (horizontal) 
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())
ggsave("Tables & Figures/BA_le_bars.pdf", BA_le_bars, width=15, height=8)



BA_le_bars1 <- ggplot(df , aes(x=factor(race, level=sort_order1), y=le, fill=as.factor(gender))) + 
  geom_col(color="black", position=position_dodge(width=0.7), width=0.7) +
  geom_linerange(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.7))+
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-1) +
  # scale_y_continuous(limits=c(df %>% ungroup() %>% select(le) %>% min,
  #                             df %>% ungroup() %>% select(le) %>% max),
  #                    oob=rescale_none)+
  scale_y_continuous(expand=expansion(mult=c(0, 0.1)))+
  scale_fill_brewer(palette="Set2")+
  labs(x="Race/Ethnicity",
       y="Life Expectancy",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  isabel_theme + 
  # emphasizes only y=lines (horizontal) 
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())
ggsave("Tables & Figures/BA_le_bars1.pdf", BA_le_bars1, width=15, height=8)




#### BA_cond: TABLES ####

# MEN FIRST - METRO
BA_metromen_cond <- BA_men_cond %>% filter(metro2=="Metropolitan") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  transmute(race, metro2, gender, le_measure, le_ci) %>% 
  spread(le_measure, le_ci) %>% 
  bind_cols(
    BA_men_cond %>% filter(metro2=="Metropolitan") %>% 
      mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                          " (",
                          format(lci, digits=1, nsmall=1),
                          ", ",
                          format(uci, digits=1, nsmall=1),
                          ")")) %>% # LE (95% CI) column 
      transmute(race, metro2, gender, lv_measure, lv) %>% 
      spread(lv_measure, lv)) %>% 
  transmute(race=`race...1`,
            metro2=`metro2...2`,
            gender=`gender...3`,
            e0, 
            CoV0=format(CoV0, digits=3, nsmall=3), 
            e10, 
            CoV10=format(CoV10, digits=3, nsmall=3), 
            e35, 
            CoV35=format(CoV35, digits=3, nsmall=3), 
            e65, 
            CoV65=format(CoV65, digits=3, nsmall=3)) %>% 
  arrange(match(race, sort_order1))
write.csv(BA_metromen_cond , "Tables & Figures/BA_metromen_cond.csv", row.names=FALSE)


# MEN FIRST - NONMETRO
BA_nonmetromen_cond <- BA_men_cond %>% filter(metro2=="Non-metro") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  transmute(race, metro2, gender, le_measure, le_ci) %>% 
  spread(le_measure, le_ci) %>% 
  bind_cols(
    BA_men_cond %>% filter(metro2=="Non-metro") %>% 
      mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                          " (",
                          format(lci, digits=1, nsmall=1),
                          ", ",
                          format(uci, digits=1, nsmall=1),
                          ")")) %>% # LE (95% CI) column 
      transmute(race, metro2, gender, lv_measure, lv) %>% 
      spread(lv_measure, lv)) %>% 
  transmute(race=`race...1`,
            metro2=`metro2...2`,
            gender=`gender...3`,
            e0, 
            CoV0=format(CoV0, digits=3, nsmall=3), 
            e10, 
            CoV10=format(CoV10, digits=3, nsmall=3), 
            e35, 
            CoV35=format(CoV35, digits=3, nsmall=3), 
            e65, 
            CoV65=format(CoV65, digits=3, nsmall=3)) %>% 
  arrange(match(race, sort_order1))
write.csv(BA_nonmetromen_cond , "Tables & Figures/BA_nonmetromen_cond.csv", row.names=FALSE)


# NOW WOMEN - METRO
BA_metrowomen_cond <- BA_women_cond %>% filter(metro2=="Metropolitan") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  transmute(race, metro2, gender, le_measure, le_ci) %>% 
  spread(le_measure, le_ci) %>% 
  bind_cols(
    BA_women_cond %>% filter(metro2=="Metropolitan") %>% 
      mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                          " (",
                          format(lci, digits=1, nsmall=1),
                          ", ",
                          format(uci, digits=1, nsmall=1),
                          ")")) %>% # LE (95% CI) column 
      transmute(race, metro2, gender, lv_measure, lv) %>% 
      spread(lv_measure, lv)) %>% 
  transmute(race=`race...1`,
            metro2=`metro2...2`,
            gender=`gender...3`,
            e0, 
            CoV0=format(CoV0, digits=3, nsmall=3), 
            e10, 
            CoV10=format(CoV10, digits=3, nsmall=3), 
            e35, 
            CoV35=format(CoV35, digits=3, nsmall=3), 
            e65, 
            CoV65=format(CoV65, digits=3, nsmall=3)) %>%   arrange(match(race, sort_order1))
write.csv(BA_metrowomen_cond, "Tables & Figures/BA_metrowomen_cond.csv", row.names=FALSE)

# WOMEN - NONMETRO
BA_nonmetro_women_cond <- BA_women_cond %>% filter(metro2=="Non-metro") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  transmute(race, metro2, gender, le_measure, le_ci) %>% 
  spread(le_measure, le_ci) %>% 
  bind_cols(
    BA_women_cond %>% filter(metro2=="Non-metro") %>% 
      mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                          " (",
                          format(lci, digits=1, nsmall=1),
                          ", ",
                          format(uci, digits=1, nsmall=1),
                          ")")) %>% # LE (95% CI) column 
      transmute(race, metro2, gender, lv_measure, lv) %>% 
      spread(lv_measure, lv)) %>% 
  transmute(race=`race...1`,
            metro2=`metro2...2`,
            gender=`gender...3`,
            e0, 
            CoV0=format(CoV0, digits=3, nsmall=3), 
            e10, 
            CoV10=format(CoV10, digits=3, nsmall=3), 
            e35, 
            CoV35=format(CoV35, digits=3, nsmall=3), 
            e65, 
            CoV65=format(CoV65, digits=3, nsmall=3)) %>%   arrange(match(race, sort_order1))
write.csv(BA_nonmetro_women_cond , "Tables & Figures/BA_nonmetro_women_cond.csv", row.names=FALSE)



#### BA_cond: FIGURES ####
# LIFE EXPECTANCY FIRST
# MEN COND FIGURE 
df <- BA_men_cond %>% mutate(age=case_when(le_measure=="e0" ~ 0,
                                     le_measure=="e10" ~ 10,
                                     le_measure=="e35" ~ 35,
                                     le_measure=="e65" ~ 65)) 

BA_men_cond_le <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
       aes(age, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH American Indian/Alaskan Native", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_x_continuous(breaks=c(0,10,35,65))+ 
  labs(x="Age",
       y="Life Expectancy among Men",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# WOMEN COND FIGURE 
df <- BA_women_cond %>% mutate(age=case_when(le_measure=="e0" ~ 0,
                                           le_measure=="e10" ~ 10,
                                           le_measure=="e35" ~ 35,
                                           le_measure=="e65" ~ 65)) 

BA_women_cond_le <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                         aes(age, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH American Indian/Alaskan Native", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_x_continuous(breaks=c(0,10,35,65))+ 
  labs(x="Age",
       y="Life Expectancy among Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


BA_cond_le <- ggarrange(BA_men_cond_le, BA_women_cond_le, nrow=2, 
                          common.legend = TRUE, 
                          legend="bottom")
ggsave("Tables & Figures/BA_cond_le.pdf", BA_cond_le, width=15, height=15)

# NOW LIFESPAN VARIATION
# MEN COND FIGURE 
df <- BA_men_cond %>% mutate(age=case_when(le_measure=="e0" ~ 0,
                                           le_measure=="e10" ~ 10,
                                           le_measure=="e35" ~ 35,
                                           le_measure=="e65" ~ 65)) 

BA_men_cond_lv <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                         aes(le, lv, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH American Indian/Alaskan Native", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  #scale_x_continuous(breaks=c(0,10,35,65))+ 
  geom_text_repel(aes(label=age))+
  labs(x="LE",
       y="Lifespan Variation among Men",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# WOMEN COND FIGURE 
df <- BA_women_cond %>% mutate(age=case_when(le_measure=="e0" ~ 0,
                                             le_measure=="e10" ~ 10,
                                             le_measure=="e35" ~ 35,
                                             le_measure=="e65" ~ 65)) 

BA_women_cond_lv <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                           aes(age, lv, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH American Indian/Alaskan Native", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_x_continuous(breaks=c(0,10,35,65))+ 
  labs(x="Age",
       y="Lifespan Variation among Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


BA_cond_lv <- ggarrange(BA_men_cond_lv, BA_women_cond_lv, nrow=2, 
                        common.legend = TRUE, 
                        legend="bottom")
ggsave("Tables & Figures/BA_cond_lv.pdf", BA_cond_lv, width=15, height=15)




#### A1_menle_table : TABLE #### 
df <- results_longitudinal %>% filter(gender=="Men") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(year5, race, metro2, le_ci) %>% 
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) 

A1_menle_table <-
  df %>% select(gender, year5, race, metro_le_ci) %>% 
  spread(year5, metro_le_ci) %>% 
  rename(metro_1990_1994=`1990-1994`,
         metro_1995_1999=`1995-1999`,
         metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(gender, year5, race, nonmetro_le_ci) %>% 
      spread(year5, nonmetro_le_ci) %>% 
      rename(nonmetro_1990_1994=`1990-1994`,
             nonmetro_1995_1999=`1995-1999`,
             nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`)) %>% 
  arrange(match(race, sort_order1))
write.csv(A1_menle_table, "Tables & Figures/A1_menle_table.csv", row.names=FALSE)


#### A1_womenle_table : TABLE #### 
df <- results_longitudinal %>% filter(gender=="Women") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(year5, race, metro2, le_ci) %>% 
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) 

A1_womenle_table <-
  df %>% select(gender, year5, race, metro_le_ci) %>% 
  spread(year5, metro_le_ci) %>% 
  rename(metro_1990_1994=`1990-1994`,
         metro_1995_1999=`1995-1999`,
         metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(gender, year5, race, nonmetro_le_ci) %>% 
      spread(year5, nonmetro_le_ci) %>% 
      rename(nonmetro_1990_1994=`1990-1994`,
             nonmetro_1995_1999=`1995-1999`,
             nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`)) %>% 
  arrange(match(race, sort_order1))
write.csv(A1_womenle_table, "Tables & Figures/A1_womenle_table.csv", row.names=FALSE)


#### A1_menlv_table : TABLE#### 
# Analysis 1 table: lv longitudinal by race/urbanicity
df <- results_longitudinal %>% filter(gender=="Men") %>% 
  transmute(year5, race, metro2, lv=as.numeric(format(lv, digits=3, nsmall=4))) %>% 
  spread(metro2, lv) %>% 
  rename(metro_lv=`Metropolitan`,
         nonmetro_lv=`Non-metro`) 

A1_menlv_table <-
  df %>% select(year5, race, metro_lv) %>% 
  spread(year5, metro_lv) %>% 
  rename(metro_1990_1994=`1990-1994`,
         metro_1995_1999=`1995-1999`,
         metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(year5, race, nonmetro_lv) %>% 
      spread(year5, nonmetro_lv) %>% 
      rename(nonmetro_1990_1994=`1990-1994`,
             nonmetro_1995_1999=`1995-1999`,
             nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`)) %>% 
  arrange(match(race, sort_order1))
write.csv(A1_menlv_table, "Tables & Figures/A1_menlv_table.csv", row.names=FALSE)


#### A1_womenlv_table : TABLE#### 
# Analysis 1 table: lv longitudinal by race/urbanicity
df <- results_longitudinal %>% filter(gender=="Women") %>% 
  transmute(year5, race, metro2, lv=as.numeric(format(lv, digits=3, nsmall=4))) %>% 
  spread(metro2, lv) %>% 
  rename(metro_lv=`Metropolitan`,
         nonmetro_lv=`Non-metro`) 

A1_womenlv_table <-
  df %>% select(year5, race, metro_lv) %>% 
  spread(year5, metro_lv) %>% 
  rename(metro_1990_1994=`1990-1994`,
         metro_1995_1999=`1995-1999`,
         metro_2000_2004=`2000-2004`,
         metro_2005_2009=`2005-2009`,
         metro_2010_2014=`2010-2014`,
         metro_2015_2019=`2015-2019`) %>% 
  full_join(
    df %>% select(year5, race, nonmetro_lv) %>% 
      spread(year5, nonmetro_lv) %>% 
      rename(nonmetro_1990_1994=`1990-1994`,
             nonmetro_1995_1999=`1995-1999`,
             nonmetro_2000_2004=`2000-2004`,
             nonmetro_2005_2009=`2005-2009`,
             nonmetro_2010_2014=`2010-2014`,
             nonmetro_2015_2019=`2015-2019`)) %>% 
  arrange(match(race, sort_order1))
write.csv(A1_womenlv_table, "Tables & Figures/A1_womenlv_table.csv", row.names=FALSE)


#### A1_le_longtrends : FIGURE ####
# MEN LIFE EXPECTANCY TRENDS
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4))
df <- df %>% mutate(year5=factor(year5, 
                                 levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                                 labels=c("1990- 1994",
                                          "1995- 1999",
                                          "2000- 2004", 
                                          "2005- 2009",
                                          "2010- 2014",
                                          "2015- 2019"))) %>% 
  filter(gender=="Men")

# le longitudinal figure 
menle_longtrend <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets,
                       aes(year5, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_y_continuous(limits=c(65,90),
                     oob=rescale_none)+
  facet_wrap(~metro2)+
  geom_text_repel(data=filter(df, year5%in%"2015- 2019"),
                  aes(label=paste0(le=format(le, digits=1, nsmall=1)), color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Life Expectancy for Men",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# WOMEN LIFE EXPECTANCY TRENDS
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4))
df <- df %>% mutate(year5=factor(year5, 
                                 levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                                 labels=c("1990- 1994",
                                          "1995- 1999",
                                          "2000- 2004", 
                                          "2005- 2009",
                                          "2010- 2014",
                                          "2015- 2019"))) %>% 
  filter(gender=="Women")

# le longitudinal figure 
womenle_longtrend <-ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets,
                           aes(year5, le, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_y_continuous(limits=c(65,90),
                     oob=rescale_none)+
  facet_wrap(~metro2)+
  geom_text_repel(data=filter(df, year5%in%"2015- 2019"),
                  aes(label=paste0(le=format(le, digits=1, nsmall=1)), color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Life Expectancy for Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# arranging figures above into one figure with two rows (gender) and two cols (metro)
A1_le_longtrends<- ggarrange(menle_longtrend, womenle_longtrend, nrow=2, 
                               common.legend = TRUE, 
                               legend="bottom")
ggsave("Tables & Figures/A1_le_longtrends.pdf", A1_le_longtrends, width=18, height=15)


# UB tests: longitudinal figure showing disparities vs overall
results_longitudinal %>% 
  pivot_wider(id_cols=c(year5, metro2, gender), names_from=race, values_from=le) %>% 
  pivot_longer(cols=H:NHW, names_to = "race") %>% 
  mutate(disp=value-Overall,
         race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))) %>% 
  ggplot(aes(year5, disp, group=race)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  #geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(color=race))+
  facet_grid(gender~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Life Expectancy Disparity vs Overall",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
results_longitudinal %>% 
  pivot_wider(id_cols=c(year5, metro2, gender), names_from=race, values_from=le) %>% 
  pivot_longer(cols=H:NHW, names_to = "race") %>% 
  mutate(disp=value-Overall,
         race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))) %>% 
  ggplot(aes(year5, disp, group=metro2)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=metro2), size=3) +
  #geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(color=metro2))+
  facet_grid(gender~race)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Life Expectancy Disparity vs Overall",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
results_longitudinal %>% 
  mutate(lv=lv/le) %>% 
  pivot_wider(id_cols=c(year5, metro2, gender), names_from=race, values_from=lv) %>% 
  pivot_longer(cols=H:NHW, names_to = "race") %>% 
  mutate(disp=value-Overall,
         race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))) %>% 
  ggplot(aes(year5, disp, group=metro2)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=metro2), size=3) +
  #geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_line(aes(color=metro2))+
  facet_grid(gender~race)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Lifespan Variation Disparity vs Overall",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())


#### A1_lv_longtrends : FIGURE ####
# MEN LIFESPAN VARIATION TRENDS
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4))
df <- df %>% mutate(year5=factor(year5, 
                                 levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                                 labels=c("1990- 1994",
                                          "1995- 1999",
                                          "2000- 2004", 
                                          "2005- 2009",
                                          "2010- 2014",
                                          "2015- 2019"))) %>% 
  filter(gender=="Men")
menlv_longtrend <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
                       aes(year5, lv, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_y_continuous(limits=c(0.15,0.33),
                     oob=rescale_none)+
  geom_text_repel(data=filter(df, year5%in%"2015- 2019"),
                  aes(label=format(lv, digits=3, nsmall=3), color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Lifespan Variation in Men",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# WOMEN LIFESPAN VARIATION TRENDS
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4))
df <- df %>% mutate(year5=factor(year5, 
                                 levels=c(1990, 1995, 2000, 2005, 2010, 2015),
                                 labels=c("1990- 1994",
                                          "1995- 1999",
                                          "2000- 2004", 
                                          "2005- 2009",
                                          "2010- 2014",
                                          "2015- 2019"))) %>% 
  filter(gender=="Women")
womenlv_longtrend <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
                            aes(year5, lv, group=race)) +
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols)+
  scale_y_continuous(limits=c(0.15,0.33),
                     oob=rescale_none)+
  geom_text_repel(data=filter(df, year5%in%"2015- 2019"),
                  aes(label=format(lv, digits=3, nsmall=3), color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_x_discrete(labels=label_wrap(10))+ 
  labs(x="Year",
       y="Lifespan Variation in Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())+
  guides(linetype=FALSE, color=guide_legend(nrow=1))


# arranging figures above into one figure with two rows
A1_lv_longtrends<- ggarrange(menlv_longtrend, womenlv_longtrend, nrow=2, 
                               common.legend = TRUE, 
                               legend="bottom")
ggsave("Tables & Figures/A1_lv_longtrends.pdf", A1_lv_longtrends, width=18, height=15)



#### A1_le_vs_lv : FIGURE #### 
# Analysis 1 figure: association b/w LE and LV by race/ethnicity, 1990-2019
# MEN ASSOCIATION
# df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
#   filter(gender=="Men") %>% 
#   select(year5, race, metro2, le, lci, uci, lv) %>% 
#   arrange(race, metro2, year5)
# 
# A1_menle_vs_lv <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
#                       aes(lv, le, group=race)) +
#   geom_point(aes(color=race), size=3, alpha=ifelse(df$year5%in%c("1990", "1995", "2000", "2005"), 0.5, 1)) +
#   geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
#   geom_path(aes(color=race),  alpha=ifelse(df$year5%in%c("1990", "1995", "2000", "2005"), 0.5, 1))+
#   geom_text_repel(data=filter(df, year5%in%c("2010", "2015")),
#                   aes(label=year5))+
#   facet_wrap(~metro2)+
#   scale_fill_manual(values=race_cols)+
#   scale_color_manual(labels=c("Overall", "White", "Black", "American Indian/Alaskan Native", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
#   labs(x="Lifespan Variation in Men",
#        y="Life Expectancy in Men",
#        color="", fill="")+
#   isabel_theme+
#   theme(legend.position="bottom", legend.title=element_blank())+
#   guides(linetype=FALSE, color=guide_legend(nrow=1))
# 
# 
# # WOMEN ASSOCIATION
# df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
#   filter(gender=="Women") %>% 
#   select(year5, race, metro2, le, lci, uci, lv) %>% 
#   arrange(race, metro2, year5) 
# 
# A1_womenle_vs_lv <-  ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
#                             aes(lv, le, group=race)) +
#   geom_point(aes(color=race), size=3, alpha=ifelse(df$year5%in%c("1990", "1995", "2000", "2005"), 0.5, 1)) +
#   geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
#   geom_path(aes(color=race),  alpha=ifelse(df$year5%in%c("1990", "1995", "2000", "2005"), 0.5, 1))+
#   geom_text_repel(data=filter(df, year5%in%c("2010", "2015")),
#                   aes(label=year5))+
#   facet_wrap(~metro2)+
#   scale_fill_manual(values=race_cols)+
#   scale_color_manual(labels=c("Overall", "White", "Black", "American Indian/Alaskan Native", "Asian/Pacific Islander", "Hispanic"), values=race_cols)+
#   labs(x="Lifespan Variation in Women",
#        y="Life Expectancy in Women",
#        color="", fill="")+
#   isabel_theme+
#   theme(legend.position="bottom", legend.title=element_blank())+
#   guides(linetype=FALSE, color=guide_legend(nrow=1))
# 
# A1_le_vs_lv  <- ggarrange(A1_menle_vs_lv, A1_womenle_vs_lv, nrow=2, 
#                              common.legend = TRUE, 
#                              legend="bottom")
# ggsave("Tables & Figures/A1_le_vs_lv.pdf", A1_le_vs_lv, width=15, height=15)

# A1_le_vs_lv is LE vs LV (coefficient of variation)
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  #filter(gender=="Men") %>% 
  select(year5, gender, race, metro2, le, lci, uci, lv, sd) %>% 
  mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"),
                     labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic")),
         year5b=ifelse(year5%in%c(1990, 2015), year5, "")) %>% 
  arrange(race, metro2, year5)  # manually computing CV: isabel to fix
 
A1_le_vs_lv <- ggplot(df, aes(x=lv, y=le, group=gender)) +
  geom_point(aes(fill=race, shape=gender), size=3, color="black") +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_path(aes(color=race, linetype=gender))+
  geom_text_repel(aes(label=year5b)) +
  scale_shape_manual(values=c(21, 24), 
                     guide=guide_legend(override.aes=list(fill="black", color="black")))+
  facet_grid(metro2~race, labeller = label_wrap_gen(width=26))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Lifespan Variation (coefficient of variation)",
       y="Life Expectancy (years)",
       color="", fill="", shape="Gender", linetype="Gender")+
  guides(fill="none", color="none")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
ggsave("Tables & Figures/A1_le_vs_lv.pdf", A1_le_vs_lv, width=20, height=15)


# A1_le_vs_lv2 is LE vs LV (standard deviation)
df <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  #filter(gender=="Men") %>% 
  select(year5, gender, race, metro2, le, lci, uci, sd) %>% 
  mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"),
                     labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic")),
         year5b=ifelse(year5%in%c(1990, 2015), year5, "")) %>% 
  arrange(race, metro2, year5)  

A1_le_vs_lv2 <- ggplot(df, aes(x=sd, y=le, group=gender)) +
  geom_point(aes(fill=race, shape=gender), size=3, color="black") +
  geom_linerange(aes(ymin=lci, ymax=uci, color=race))+
  geom_path(aes(color=race, linetype=gender))+
  geom_text_repel(aes(label=year5b)) +
  scale_shape_manual(values=c(21, 24), 
                     guide=guide_legend(override.aes=list(fill="black", color="black")))+
  facet_grid(metro2~race, labeller = label_wrap_gen(width=26))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Lifespan Variation (standard deviation)",
       y="Life Expectancy (years)",
       color="", fill="", shape="Gender", linetype="Gender")+
  guides(fill="none", color="none")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
ggsave("Tables & Figures/A1_le_vs_lv2.pdf", A1_le_vs_lv2, width=20, height=15)



# A1_longdisp : FIGURE #### 
# Analysis 1 figure: metro-nonmetro disparities over time in LE by race/ethnicity, 2000-2019

# MEN METRO-NONMETRO DISPARITY IN LE & LV
metro_long_menle <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  mutate(year5=factor(year5, 
                      levels=c(1990, 1995, 2000, 2005, 2010, 2015), 
                      labels=c("1990- 1994",
                               "1995- 1999",
                               "2000- 2004", 
                               "2005- 2009",
                               "2010- 2014",
                               "2015- 2019"))) %>% 
  filter(gender=="Men") %>% 
  select(year5, race, metro2, le) %>% 
  spread(metro2, le) %>% 
  rename(metro_le=`Metropolitan`,
         nonmetro_le=`Non-metro`) %>% 
  mutate(abs_le=metro_le-nonmetro_le) 

menle_long_disp <- ggplot(metro_long_menle %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
                       aes(year5, abs_le, group=race)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_x_discrete(labels=label_wrap(10))+
  labs(x="Year",
       y=str_wrap("Metropolitan to Nonmetropolitan Disparity in Life Expectancy in Men", 38),
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
  
menle_long_disp <- menle_long_disp + 
  annotate("text", x=4.5, y=6.70, size=5, label="Observations above 0-line indicate\n higher LE among men in metro areas")+
  annotate("text", x=4.5, y=-0.50, size=5, label="Observations below 0-line indicate\n higher LE among men in non-metro areas")

metro_long_menlv <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  mutate(year5=factor(year5, 
                      levels=c(1990, 1995, 2000, 2005, 2010, 2015), 
                      labels=c("1990- 1994",
                               "1995- 1999",
                               "2000- 2004", 
                               "2005- 2009",
                               "2010- 2014",
                               "2015- 2019"))) %>% 
  filter(gender=="Men") %>% 
  select(year5, race, metro2, lv) %>% 
  spread(metro2, lv) %>% 
  rename(metro_lv=`Metropolitan`,
         nonmetro_lv=`Non-metro`) %>% 
  mutate(abs_lv=metro_lv-nonmetro_lv)

menlv_long_disp <- ggplot(metro_long_menlv %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
                       aes(year5, abs_lv, group=race)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+
  labs(x="Year",
       y=str_wrap("Metropolitan to Nonmetropolitan Disparity in Lifespan Variation in Men", 38),
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

menlv_long_disp <- menlv_long_disp + 
  annotate("text", x=4.5, y=0.05, size=5, label="Observations above 0-line indicate\n higher LV among men in metro areas")+
  annotate("text", x=4.5, y=-0.06, size=5, label="Observations below 0-line indicate\n higher LV among men in non-metro areas")

A1_menlongdisp <- ggarrange(menle_long_disp, menlv_long_disp, ncol=2,
          common.legend = TRUE,
          legend="none")


# WOMEN METRO-NONMETRO DISPARITY IN LE & LV
metro_long_womenle <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  mutate(year5=factor(year5, 
                      levels=c(1990, 1995, 2000, 2005, 2010, 2015), 
                      labels=c("1990- 1994",
                               "1995- 1999",
                               "2000- 2004", 
                               "2005- 2009",
                               "2010- 2014",
                               "2015- 2019"))) %>% 
  filter(gender=="Women") %>% 
  select(year5, race, metro2, le) %>% 
  spread(metro2, le) %>% 
  rename(metro_le=`Metropolitan`,
         nonmetro_le=`Non-metro`) %>% 
  mutate(abs_le=metro_le-nonmetro_le) 

womenle_long_disp <- ggplot(metro_long_womenle %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
                          aes(year5, abs_le, group=race)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_x_discrete(labels=label_wrap(10))+
  labs(x="Year",
       y=str_wrap("Metropolitan to Nonmetropolitan Disparity in Life Expectancy in Women", 38),
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

womenle_long_disp <- womenle_long_disp + 
  annotate("text", x=4.5, y=6.70, size=5, label="Observations above 0-line indicate\n higher LE among women in metro areas")+
  annotate("text", x=4.5, y=-0.50, size=5, label="Observations below 0-line indicate\n higher LE among women in non-metro areas")

metro_long_womenlv <- results_longitudinal %>% mutate(year5=substr(year5, 1, 4)) %>% 
  mutate(year5=factor(year5, 
                      levels=c(1990, 1995, 2000, 2005, 2010, 2015), 
                      labels=c("1990- 1994",
                               "1995- 1999",
                               "2000- 2004", 
                               "2005- 2009",
                               "2010- 2014",
                               "2015- 2019"))) %>% 
  filter(gender=="Women") %>% 
  select(year5, race, metro2, lv) %>% 
  spread(metro2, lv) %>% 
  rename(metro_lv=`Metropolitan`,
         nonmetro_lv=`Non-metro`) %>% 
  mutate(abs_lv=metro_lv-nonmetro_lv)

womenlv_long_disp <- ggplot(metro_long_womenlv %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
                          aes(year5, abs_lv, group=race)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_line(aes(linetype=race, color=race))+
  scale_linetype_manual(values=c(2,1,1,1,1,1))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  scale_x_discrete(labels=label_wrap(10))+
  labs(x="Year",
       y=str_wrap("Metropolitan to Nonmetropolitan Disparity in Lifespan Variation in Women", 40),
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

womenlv_long_disp  <- womenlv_long_disp  + 
  annotate("text", x=4.5, y=0.25, size=5, label="Observations above 0-line indicate\n higher LV among women in metro areas")+
  annotate("text", x=4.5, y=-0.20, size=5, label="Observations below 0-line indicate\n higher LV among women in non-metro areas")

A1_womenlongdisp <- ggarrange(womenle_long_disp, womenlv_long_disp, ncol=2,
                            common.legend = TRUE,
                            legend="bottom")


A1_longdisp <- ggarrange(A1_menlongdisp, A1_womenlongdisp, nrow=2,
          common.legend = TRUE,
          legend="bottom")
ggsave("Tables & Figures/A1_longdisp.pdf", A1_longdisp, width=20, height=15)




#### A2_lelv_table : TABLE #### 
# Analysis 2 table: le/lv by gender
df <- results_gender %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",                      
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`)) 

A2_le_table <- df %>% select(race, gender, metro_le_ci) %>% 
  spread(gender, metro_le_ci) %>% 
  rename(metro_men_le=`Men`,
         metro_women_le=`Women`) %>% 
  full_join(
    df %>% select(race, gender, nonmetro_le_ci) %>% 
      spread(gender, nonmetro_le_ci) %>% 
      rename(nonmetro_men_le=`Men`,
             nonmetro_women_le=`Women`)) %>% 
  full_join(
    df %>% transmute(race, gender, metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1))) %>% 
      spread(gender, metro_lv) %>% 
      rename(metro_men_lv=`Men`,
             metro_women_lv=`Women`)) %>% 
  full_join(
    df %>% transmute(race, gender, nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1))) %>% 
      spread(gender, nonmetro_lv) %>% 
      rename(nonmetro_men_lv=`Men`,
             nonmetro_women_lv=`Women`)) %>% 
  # blending baseline and gender together 
  select(-metro_men_lv:-nonmetro_women_lv) %>% 
  full_join(BA_lelv_table, by="race") %>% select(-abs_le:-rel_lv) %>% 
  select(race, metro_le_ci, nonmetro_le_ci, metro_men_le, nonmetro_men_le, metro_women_le, nonmetro_women_le) %>% 
  arrange(match(race, sort_order1))
write.csv(A2_lelv_table, "Tables & Figures/A2_le_table.csv", row.names=FALSE)

# METRO-NONMETRO ABSOLUTE RELATIVE DISPARITIES TABLE FOR BLENDED BASELINE & GENDER
df <- results_gender %>% 
  select(race, metro2, le) %>%
  spread(metro2, le) %>% 
  rename(metro_le=`Metropolitan`,
         nonmetro_le=`Non-metro`) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`)) 

A2_ledisp_table <- df %>% select(race, gender, metro_le) %>% 
  spread(gender, metro_le) %>% 
  rename(metro_men_le=`Men`,
         metro_women_le=`Women`) %>% 
  full_join(
    df %>% select(race, gender, nonmetro_le) %>% 
      spread(gender, nonmetro_le) %>% 
      rename(nonmetro_men_le=`Men`,
             nonmetro_women_le=`Women`)) %>% 
  # blending baseline and gender together 
  full_join(
    results_baseline %>% select(race, metro2, le) %>% 
      spread(metro2, le) %>% 
      rename(metro_le=`Metropolitan`,
             nonmetro_le=`Non-metro`)) %>% 
  select(race, metro_le, nonmetro_le, metro_men_le, nonmetro_men_le, metro_women_le, nonmetro_women_le) %>%
  mutate(abs_disp_mw = metro_le-nonmetro_le,
         rel_disp_mw=metro_le/nonmetro_le,
         abs_disp_m=metro_men_le-nonmetro_men_le,
         rel_disp_m=metro_men_le/nonmetro_men_le, 
         abs_disp_w=metro_women_le-nonmetro_women_le, 
         rel_disp_w=metro_women_le/nonmetro_women_le) %>% 
  arrange(match(race, sort_order1)) 


A2_lvdisp_table <- df %>% select(race, gender, metro_lv) %>% 
  spread(gender, metro_lv) %>% 
  rename(metro_men_lv=`Men`,
         metro_women_lv=`Women`) %>% 
  full_join(
    df %>% select(race, gender, nonmetro_lv) %>% 
      spread(gender, nonmetro_lv) %>% 
      rename(nonmetro_men_lv=`Men`,
             nonmetro_women_lv=`Women`)) %>% 
  # blending baseline and gender together 
full_join(
  results_baseline %>% select(race, metro2, lv) %>% 
    spread(metro2, lv) %>% 
    rename(metro_lv=`Metropolitan`,
           nonmetro_lv=`Non-metro`)) %>% 
  select(race, metro_lv, nonmetro_lv, metro_men_lv, nonmetro_men_lv, metro_women_lv, nonmetro_women_lv) %>%
  mutate(abs_disp_mw = metro_lv-nonmetro_lv,
         rel_disp_mw=metro_lv/nonmetro_lv,
         abs_disp_m=metro_men_lv-nonmetro_men_lv,
         rel_disp_m=metro_men_lv/nonmetro_men_lv, 
         abs_disp_w=metro_women_lv-nonmetro_women_lv, 
         rel_disp_w=metro_women_lv/nonmetro_women_lv) %>% 
  arrange(match(race, sort_order1)) 
write.csv(A2_lvdisp_table, "Tables & Figures/A2_lvdisp_table.csv", row.names=FALSE)


#### A2_MvF_disparities : FIGURE #### 
# Analysis 2 figure: MvF scatterplot (LE and LV)  

## first, creating LE scatterplot
# spreading LE and 95% CIs by gender 
df <- results_gender %>% full_join(
  results_gender %>% select(race, metro2, gender, le) %>% 
    spread(gender,le) %>% 
    rename(men_le=`Men`,
           women_le=`Women`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, lci) %>% 
      spread(gender,lci) %>% 
      rename(men_lci=`Men`,
             women_lci=`Women`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, uci) %>% 
      spread(gender,uci) %>% 
      rename(men_uci=`Men`,
             women_uci=`Women`)) %>% 
  select(race, metro2, gender, men_le, men_lci, men_uci, women_le, women_lci, women_uci) 
# generating scatterplot
A2_le_scatter <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
                        aes(x=men_le, y=women_le)) +
  #geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(men_le=seq(70, 90, by=1),
                                   women_le=seq(70, 90, by=1)))+
  geom_linerange(aes(ymin=women_lci, ymax=women_uci))+
  geom_linerange(aes(xmin=men_lci, xmax=men_uci))+
  geom_point(aes(fill=race, shape=metro2), color="black", size=3) +
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH American Indian/Alaskan Native", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  # usama's recommendation regarding points
  scale_shape_manual(values=c(21, 24), 
                     guide=guide_legend(override.aes=list(fill="black", color="black")))+
  ggtitle("(A) Life Expectancy, Men vs Women")+
  labs(x="Life Expectancy for Men (95% CI)",
       y="Life Expectancy for Women (95% CI)",
       color="", fill="", shape="")+
  guides(linetype=FALSE, color=guide_legend(nrow=1))+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())


## then, creating LV scatterplot
# spreading LV by gender 
df <- results_gender %>% select(race, metro2, gender, lv) %>% 
    spread(gender,lv) %>% 
    rename(men_lv=`Men`,
           women_lv=`Women`)
# generating scatterplot
A2_lv_scatter <- 
  ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets, 
         aes(x=men_lv, y=women_lv)) +
  geom_abline(intercept = 0, slope=1, lty=1)+
  # geom_line(lty=1, data=data.frame(men_lv=seq(0.18, 0.30, by=0.1),
  #                                   women_lv=seq(0.18, 0.30, by=0.1)))+
  geom_point(aes(fill=race, shape=metro2), color="black", size=3) +
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  scale_color_manual(values=race_cols)+
  # usama's recommendation regarding points
  scale_shape_manual(values=c(21, 24), 
                     guide=guide_legend(override.aes=list(fill="black", color="black")))+
  ggtitle("(B) Lifespan Variation, Men vs Women")+
  labs(x="Lifespan Variation for Men",
       y="Lifespan Variation for Women",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# arranging figures above into one figure with two panels
A2_MvF_disparities <- ggarrange(A2_le_scatter, A2_lv_scatter, ncol=2, 
          common.legend = TRUE, 
          legend="bottom")
ggsave("Tables & Figures/A2_MvF_disparities.pdf", A2_MvF_disparities, width=15, height=8)


#### A2_bland_altman: FIGURE ####
# try bland-altman plot
# x-axis=mean of the two variables
# y-axis=difference
df <- results_gender %>% full_join(
  results_gender %>% select(race, metro2, gender, le) %>% 
    spread(gender,le) %>% 
    rename(men_le=`Men`,
           women_le=`Women`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, lci) %>% 
      spread(gender,lci) %>% 
      rename(men_lci=`Men`,
             women_lci=`Women`)) %>% 
  full_join(
    results_gender %>% select(race, metro2, gender, uci) %>% 
      spread(gender,uci) %>% 
      rename(men_uci=`Men`,
             women_uci=`Women`)) %>% 
  select(race, metro2, gender, men_le, men_lci, men_uci, women_le, women_lci, women_uci) 

# bland altman plot
A2_bland_altman <- ggplot(df %>% rowwise() %>% 
                            mutate(mean=mean(c(men_le, women_le)),
                                   dif=women_le-men_le,
                                   race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                          aes(x=mean, y=dif)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(fill=race, shape=metro2), color="black", size=3) +
  geom_smooth(method=lm, se=FALSE, col="black", size=0.5)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  scale_color_manual(values=race_cols)+
  # usama's recommendation regarding points
  scale_shape_manual(values=c(21, 24), 
                     guide=guide_legend(override.aes=list(fill="black", color="black")))+
  labs(x="Average Life Expectancy, Men and Women",
       y="Life Expectancy Difference (Women-Men)",
       color="", fill="", shape="")+
  #guides(color=F)+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
ggsave("Tables & Figures/A2_bland_altman.pdf", A2_bland_altman, width=10, height=8)



#### A3_le_table : TABLE #### 
# Analysis 3 table: le by race/urbanicity (6 metro categories)
# MEN FIRST
A3_menle_table <- results_urbanicity %>% filter(gender=="Men") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro6, gender, le_ci) %>% 
  spread(metro6, le_ci) %>% 
  arrange(match(race, sort_order1)) %>% 
  bind_cols(., A1_menle_table) %>% 
  transmute(race=`race...1`, 
            `All metro counties`= metro_2015_2019,
            `Large central`, 
            `Large fringe`, 
            `Medium`,
            `Small`, 
            `Micro`,
            `Noncore`,
            `All nonmetro counties`= nonmetro_2015_2019) %>% 
  arrange(match(race, sort_order1))
write.csv(A3_menle_table, "Tables & Figures/A3_menle_table.csv", row.names=FALSE)

# WOMEN NEXT
A3_womenle_table <- results_urbanicity %>% filter(gender=="Women") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, metro6, gender, le_ci) %>% 
  spread(metro6, le_ci) %>% 
  arrange(match(race, sort_order1)) %>% 
  bind_cols(., A1_womenle_table) %>% 
  transmute(race=`race...1`, 
            `All metro counties`= metro_2015_2019,
            `Large central`, 
            `Large fringe`, 
            `Medium`,
            `Small`, 
            `Micro`,
            `Noncore`,
            `All nonmetro counties`= nonmetro_2015_2019) %>% 
  arrange(match(race, sort_order1))
write.csv(A3_womenle_table, "Tables & Figures/A3_womenle_table.csv", row.names=FALSE)


#### A3_lv_table : TABLE #### 
# Analysis 3 table: lv by race/urbanicity (6 metro categories)
# MEN FIRST
A3_menlv_table <- results_urbanicity %>% filter(gender=="Men") %>% 
  select(race, metro6, gender, lv) %>% 
  spread(metro6, lv) %>% 
  arrange(match(race, sort_order1)) %>% 
  bind_cols(., A1_menlv_table) %>% 
  transmute(race=`race...1`, 
            `All metro counties`= metro_2015_2019,
            `Large central`, 
            `Large fringe`, 
            `Medium`,
            `Small`, 
            `Micro`,
            `Noncore`,
            `All nonmetro counties`= nonmetro_2015_2019) %>% 
  arrange(match(race, sort_order1))
write.csv(A3_menlv_table, "Tables & Figures/A3_menlv_table.csv", row.names=FALSE)

# WOMEN NEXT
A3_womenlv_table <- results_urbanicity %>% filter(gender=="Women") %>% 
  select(race, metro6, gender, lv) %>% 
  spread(metro6, lv) %>% 
  arrange(match(race, sort_order1)) %>% 
  bind_cols(., A1_womenlv_table) %>% 
  transmute(race=`race...1`, 
            `All metro counties`= metro_2015_2019,
            `Large central`, 
            `Large fringe`, 
            `Medium`,
            `Small`, 
            `Micro`,
            `Noncore`,
            `All nonmetro counties`= nonmetro_2015_2019) %>% 
  arrange(match(race, sort_order1))
write.csv(A3_womenlv_table, "Tables & Figures/A3_womenlv_table.csv", row.names=FALSE)




#### A3_le_fig : FIGURE ####
### trying out dot plot instead of bars ### 
A3_le_fig <- ggplot(results_urbanicity %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"),
                                                              labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"))), # reorders facets
       aes(x=metro6, y=le, group=1)) + 
  geom_point(aes(color=race, shape=gender), size=3) +
  facet_wrap(~race, nrow=1, labeller = label_wrap_gen(width=26))+
  geom_line(aes(linetype=gender, group=gender, color=race))+
  geom_linerange(aes(ymin=lci, ymax=uci, color=race), position=position_dodge(width=0.7))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.7), vjust=-0.5) +
  scale_y_continuous(limits=c(results_urbanicity %>% ungroup() %>% select(le) %>% min,
                              results_urbanicity %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_x_discrete(labels=function(results_urbanicity) str_wrap(results_urbanicity, width = 10))+
  labs(x="Levels of Urbanization",
       y="Life Expectancy",
       color="", fill="", shape="Gender", linetype="Gender")+
  guides(fill="none", color="none")+guides(fill="none", color="none")+
  isabel_theme +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())


#### A3_lv_fig : FIGURE ####
# Analysis 3 figure: one panel per race, one bar per metro category
A3_lv_fig <- ggplot(results_urbanicity %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"),
                                                              labels=c("Overall", "White", "Black", "Asian/Pacific Islander", "Hispanic"))), # reorders facets
                    aes(x=metro6, y=lv, group=1)) + 
  geom_point(aes(color=race, shape=gender), size=3) +
  facet_wrap(~race, nrow=1, labeller = label_wrap_gen(width=26))+
  geom_line(aes(linetype=gender, group=gender, color=race))+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  geom_text(aes(label=format(lv, digits=3, nsmall=3)), position=position_dodge(width=0.7), vjust=-0.5) +
  scale_y_continuous(limits=c(results_urbanicity %>% ungroup() %>% select(lv) %>% min,
                              results_urbanicity %>% ungroup() %>% select(lv) %>% max),
                     oob=rescale_none)+
  scale_x_discrete(labels=function(results_urbanicity) str_wrap(results_urbanicity, width = 10))+
  labs(x="Levels of Urbanization",
       y="Lifespan Variation",
       color="", fill="", shape="Gender", linetype="Gender")+
  guides(fill="none", color="none")+guides(fill="none", color="none")+
  isabel_theme +
  theme(legend.position="bottom", legend.title=element_blank())+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())


# arranging figures above into one figure 
A3_le_lv<- ggarrange(A3_le_fig, A3_lv_fig, nrow=2, 
                             common.legend = TRUE, 
                             legend="bottom")
ggsave("Tables & Figures/A3_le_lv.pdf", A3_le_lv, width=20, height=15)


#### A4_le_table : TABLE ####
# Analysis 4 table: le & lv by race/urbanicity and census region 

# MEN FIRST 
df <- results_geography %>% filter(gender=="Men") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, gender, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_geography %>% filter (gender=="Men") %>% 
      select(race, metro2, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`))

A4_menle_table <- df %>% select(census_region, race, metro_le_ci) %>% 
  spread(census_region, metro_le_ci) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% select(census_region, race, nonmetro_le_ci) %>% 
      spread(census_region, nonmetro_le_ci) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(race,
         NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro) %>% 
  arrange(match(race, sort_order1))
write.csv(A4_menle_table, "Tables & Figures/A4_menle_table.csv", row.names=FALSE)

A4_menlv_table <- df %>% transmute(census_region, race, metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1))) %>% 
  spread(census_region, metro_lv) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% transmute(census_region, race, nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1))) %>% 
      spread(census_region, nonmetro_lv) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(race,NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro) %>% 
  arrange(match(race, sort_order1))
write.csv(A4_menlv_table, "Tables & Figures/A4_menlv_table.csv", row.names=FALSE)


# NOW WOMEN 
df <- results_geography %>% filter(gender=="Women") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column 
  select(race, gender, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_geography %>% filter (gender=="Women") %>% 
      select(race, metro2, lv) %>% 
      spread(metro2, lv) %>% 
      rename(metro_lv=`Metropolitan`,
             nonmetro_lv=`Non-metro`))

A4_womenle_table <- df %>% select(census_region, race, metro_le_ci) %>% 
  spread(census_region, metro_le_ci) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% select(census_region, race, nonmetro_le_ci) %>% 
      spread(census_region, nonmetro_le_ci) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(race,
         NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro) %>% 
  arrange(match(race, sort_order1))
write.csv(A4_womenle_table, "Tables & Figures/A4_womenle_table.csv", row.names=FALSE)

A4_womenlv_table <- df %>% transmute(census_region, race, metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1))) %>% 
  spread(census_region, metro_lv) %>% 
  rename(NE_metro=`Northeast`,
         MW_metro=`Midwest`,
         S_metro=`South`,
         W_metro=`West`) %>% 
  full_join(
    df %>% transmute(census_region, race, nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1))) %>% 
      spread(census_region, nonmetro_lv) %>% 
      rename(NE_nonmetro=`Northeast`,
             MW_nonmetro=`Midwest`,
             S_nonmetro=`South`,
             W_nonmetro=`West`)) %>% 
  select(race,NE_metro, NE_nonmetro,
         MW_metro, MW_nonmetro,
         S_metro, S_nonmetro, 
         W_metro, W_nonmetro) %>% 
  arrange(match(race, sort_order1))
write.csv(A4_womenlv_table, "Tables & Figures/A4_womenlv_table.csv", row.names=FALSE)



#### A4_le_bars : FIGURE ####
# Analysis 4 figure: le comparison within regions and metro between races
# MEN FIRST
A4_menle_bars <- ggplot(results_geography %>% filter(gender=="Men") %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                      aes(x=census_region, y=le, fill=as.factor(race))) + 
  #geom_col(color="black", position=position_dodge(width=0.8), width=0.8) +
  geom_bar_pattern(position = position_dodge(preserve="single"),
                   aes(pattern_density=race, fill=race),
                   color = "black", 
                   stat="identity",
                   pattern_fill = "white",
                   pattern_angle = 45,
                   #pattern_density = 0,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6) + 
  scale_pattern_density_manual(values=c(0.6, 0, 0, 0, 0, 0))+
  geom_linerange(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.8))+
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Census Region",
       y="Life Expectancy among Men",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  guides(pattern_density=F)+ # only want stripes to show on "Overall" legend
  isabel_theme+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())

# NOW WOMEN
A4_womenle_bars <- ggplot(results_geography %>% filter(gender=="Women") %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                        aes(x=census_region, y=le, fill=as.factor(race))) + 
  #geom_col(color="black", position=position_dodge(width=0.8), width=0.8) +
  geom_bar_pattern(position = position_dodge(preserve="single"),
                   aes(pattern_density=race, fill=race),
                   color = "black", 
                   stat="identity",
                   pattern_fill = "white",
                   pattern_angle = 45,
                   #pattern_density = 0,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6) + 
  scale_pattern_density_manual(values=c(0.6, 0, 0, 0, 0, 0))+
  geom_linerange(aes(ymin=lci, ymax=uci), position=position_dodge(width=0.8))+
  geom_text(aes(label=format(le, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(le) %>% min,
                              results_geography %>% ungroup() %>% select(le) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Census Region",
       y="Life Expectancy among Women",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  guides(pattern_density=F)+ # only want stripes to show on "Overall" legend
  isabel_theme+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())

A4_le_bars <- ggarrange(A4_menle_bars, A4_womenle_bars, nrow=2, 
                       common.legend = TRUE, 
                       legend="bottom")
ggsave("Tables & Figures/A4_le_bars.pdf", A4_le_bars , width=17, height=13)


#### A4_lv_bars3 : FIGURE ####

# MEN FIRST 
A4_menlv_bars <- ggplot(results_geography %>% filter(gender=="Men") %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                      aes(x=census_region, y=lv, fill=as.factor(race))) + 
  #geom_col(color="black", position=position_dodge(width=0.8), width=0.8) +
  geom_bar_pattern(position = position_dodge(preserve="single"),
                   aes(pattern_density=race, fill=race),
                   color = "black", 
                   stat="identity",
                   pattern_fill = "white",
                   pattern_angle = 45,
                   #pattern_density = 0,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6) + 
  scale_pattern_density_manual(values=c(0.6, 0, 0, 0, 0, 0))+
  geom_text(aes(label=format(lv, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(lv) %>% min,
                              results_geography %>% ungroup() %>% select(lv) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Census Region",
       y="Lifespan Variation among Men",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  guides(pattern_density=F)+ # only want stripes to show on "Overall" legend
  isabel_theme+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())

# NOW WOMEN
A4_womenlv_bars <- ggplot(results_geography %>% filter(gender=="Women") %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAIAN", "NHAPI", "H"))), # reorders facets,
                        aes(x=census_region, y=lv, fill=as.factor(race))) + 
  #geom_col(color="black", position=position_dodge(width=0.8), width=0.8) +
  geom_bar_pattern(position = position_dodge(preserve="single"),
                   aes(pattern_density=race, fill=race),
                   color = "black", 
                   stat="identity",
                   pattern_fill = "white",
                   pattern_angle = 45,
                   #pattern_density = 0,
                   pattern_spacing = 0.01,
                   pattern_key_scale_factor = 0.6) + 
  scale_pattern_density_manual(values=c(0.6, 0, 0, 0, 0, 0))+
  geom_text(aes(label=format(lv, digits=1, nsmall=1)), position=position_dodge(width=0.8), vjust=-1) +
  scale_y_continuous(limits=c(results_geography %>% ungroup() %>% select(lv) %>% min,
                              results_geography %>% ungroup() %>% select(lv) %>% max),
                     oob=rescale_none)+
  scale_fill_manual(values=race_cols)+
  scale_color_manual(values=race_cols)+
  labs(x="Census Region",
       y="Lifespan Variation among Women",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  facet_wrap(~metro2)+
  guides(pattern_density=F)+ # only want stripes to show on "Overall" legend
  isabel_theme+
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())

A4_lv_bars <- ggarrange(A4_menlv_bars, A4_womenlv_bars, nrow=2, 
                        common.legend = TRUE, 
                        legend="bottom")
ggsave("Tables & Figures/A4_lv_bars .pdf", A4_lv_bars  , width=17, height=13)



#### BASIC DESCRIPTIVE TABLE ####
# SELECTED CHARACTERISTICS FOR URBAN-RURAL CATEGORIES
descrip_dta<- master_dta %>% 
  mutate(gender=as.numeric(sex),
         age_5yr_group=as.numeric(age_5yr_group),
         year5=case_when(
           year%in%c(2015:2019) ~ as.character('2015-2019'),
           year%in%c(2010:2014) ~ as.character('2010-2014'),
           year%in%c(2005:2009) ~ as.character('2005-2009'),
           year%in%c(2000:2004) ~ as.character('2000-2004'),
           year%in%c(1995:1999) ~ as.character('1995-1999'),
           year%in%c(1990:1994) ~ as.character('1990-1994'))) %>% 
  # factoring metro2, metro6, gender, census_region to help in ggplot later
  mutate(metro2=factor(metro2, levels=c(1,0),
                       labels=c("Metropolitan", "Non-metro")),
         metro6=factor(metro6, levels=c(1:6),
                       labels=c("Large central", "Large fringe", "Medium", "Small", "Micro", "Noncore")),
         gender=factor(gender, levels=c(1,0),
                    labels=c("Male", "Female")),
         census_region=factor(census_region, levels=c(1:4),
                              labels=c("Northeast", "Midwest", "South", "West")))



### OVERALL COLUMN
# finding # of counties (metro v nonmetro) then summing for total # of counties
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(fips, metro2) %>% 
  distinct() %>% 
  group_by(metro2) %>% 
  count()

# finding total # pop (2015-2019)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(pop_denom) %>% sum()/5

# finding total # deaths (2015-2019)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(count) %>% sum()/5

# finding total population by race %
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(race, pop_denom) %>% 
  group_by(race) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  mutate(pop_denom=pop_denom/5) %>% 
  mutate(pop_share=(pop_denom/321968939)*100)

# finding total death count by race %
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(race, count) %>% 
  group_by(race) %>% 
  summarise(count=sum(count)) %>% 
  mutate(count=count/5) %>% 
  mutate(death_share=(count/2770397)*100) 


### OVERALL METRO, OVERALL NONMETRO COLUMNS
# finding total population size of counties (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, pop_denom) %>% 
  group_by(metro2) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  mutate(yearly_pop=pop_denom/5)

# finding total death counts of counties (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, count) %>% 
  group_by(metro2) %>% 
  summarise(count=sum(count)) %>% 
  mutate(yearly_deaths=count/5)

# finding total population by race and share of racial % (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, race, pop_denom) %>% 
  group_by(metro2, race) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  spread(metro2, pop_denom) %>% 
  mutate(metro_share=((`Metropolitan`/5)/276992013)*100,
         nonmetro_share=((`Non-metro`/5)/44976925)*100)

# finding total deaths by race and share of racial death % (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, race, count) %>% 
  group_by(metro2, race) %>% 
  summarise(count=sum(count)) %>% 
  spread(metro2, count) %>% 
  mutate(metro_share=((`Metropolitan`/5)/2255180)*100,
         nonmetro_share=((`Non-metro`/5)/515216)*100)


### 6 URBANICITY COLUMNS
# finding # of counties (6 urbanicity categories)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(fips, metro6) %>% 
  distinct() %>% 
  group_by(metro6) %>% 
  count()

# finding total population size of counties (6 urbanicity categories)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro6, pop_denom) %>% 
  group_by(metro6) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  mutate(yearly_pop=pop_denom/5)

# finding total death counts of counties (6 urbanicity categories)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro6, count) %>% 
  group_by(metro6) %>% 
  summarise(count=sum(count)) %>% 
  mutate(yearly_deaths=count/5)

# finding total population by race and share of racial % (6 urbanicity categories)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro6, race, pop_denom) %>% 
  group_by(metro6, race) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  spread(metro6, pop_denom) %>% 
  transmute(race,
            `Large central`=((`Large central`/5)/99599953)*100,
            `Large fringe`=((`Large fringe`/5)/80781822)*100,
            `Medium`=((`Medium`/5)/67364296)*100,
            `Small`=((`Small`/5)/29245942)*100,
            `Micro`=((`Micro`/5)/26713192)*100,
            `Noncore`=((`Noncore`/5)/18263733)*100)

# finding total death count by race and share of racial death % (6 urbanicity categories)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro6, race, count) %>% 
  group_by(metro6, race) %>% 
  summarise(count=sum(count)) %>% 
  spread(metro6, count) %>% 
  transmute(race,
            `Large central`=((`Large central`/5)/719257)*100,
            `Large fringe`=((`Large fringe`/5)/644833)*100,
            `Medium`=((`Medium`/5)/604345)*100,
            `Small`=((`Small`/5)/286745)*100,
            `Micro`=((`Micro`/5)/291408)*100,
            `Noncore`=((`Noncore`/5)/223809)*100)



# finding total population size of counties (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, census_region, pop_denom) %>% 
  group_by(census_region, metro2) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  mutate(yearly_pop=pop_denom/5) 

# finding total population by race and share of racial % (metro v nonmetro)
descrip_dta %>% filter(year5%in%"2015-2019") %>% 
  select(metro2, census_region, race, pop_denom) %>% 
  group_by(census_region, metro2, race) %>% 
  summarise(pop_denom=sum(pop_denom)) %>% 
  spread(metro2, pop_denom) %>% 
  mutate(metro_share=case_when(
      census_region=="Northeast" ~ ((`Metropolitan`/5)/51469549)*100,
      census_region=="Midwest" ~ ((`Metropolitan`/5)/53067409)*100,
      census_region=="South" ~ ((`Metropolitan`/5)/103702193)*100,
      census_region=="West" ~ ((`Metropolitan`/5)/70403480)*100),
    nonmetro_share=case_when(
      census_region=="Northeast" ~ ((`Non-metro`/5)/4563587)*100,
      census_region=="Midwest" ~ ((`Non-metro`/5)/15040698)*100,
      census_region=="South" ~ ((`Non-metro`/5)/19706041)*100,
      census_region=="West" ~ ((`Non-metro`/5)/6744851)*100))



#### EXPLORATORY STUFF ####
tmp <- dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=(count/pop_denom), 
         age_5yr_group=factor(as.numeric(age_5yr_group))) %>% 
  arrange(metro2, race, age_5yr_group) 

# TREND LINE FIGURE OF AGE-SPECIFIC MORTALITY RATES
ggplot(tmp %>% mutate(race=factor(race, levels=c("NHW", "NHB", "NHAIAN", "NHAPI", "H"))), 
       aes(x=age_5yr_group, y=mx)) + 
  geom_point(aes(color=race), size=3) +
  geom_line(aes(color=race))+
  facet_wrap(~metro2)+
  scale_fill_manual(values=race_cols[2:6])+
  scale_color_manual(values=race_cols[2:6])+
  labs(x="age-specific mortality rate",
       y="5-year age group",
       color="", fill="")+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())

# BAR GRAPH FIGURE OF AGE-SPECIFIC MORTALITY RATES -- NONMETRO
ggplot(tmp %>% filter(metro2%in%"Non-metro") %>% mutate(race=factor(race, levels=c("NHW", "NHB", "NHAIAN", "NHAPI", "H"))), 
       aes(x=age_5yr_group, y=mx, fill=as.factor(race))) + 
  geom_col(color="black", position=position_dodge(width=0.7), width=0.7) +
  scale_y_continuous(limits=c(0,0.05))+
  scale_fill_manual(values=race_cols[2:6])+
  scale_color_manual(values=race_cols[2:6])+
  labs(x="age-specific mortality rate",
       y="5-year age group",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  isabel_theme + 
  # emphasizes only y=lines (horizontal) 
  theme(panel.grid.major.y=element_line(color="gray60")) +
  theme(panel.grid.major.x=element_blank())



# ATTEMPT AT USAMA'S FIGURE FOR AGE SPECIFIC MORTALITY RATES
agegrp_order <- c("0-4 years", "5-14 years", "15-24 years", "25-34 years", "35-44 years",
                "45-54 years", "55-64 years","65-74 years", "75-84 years", "85+ years")

tmp <- dta %>% filter(year5=='2015-2019') %>% 
  mutate(age_5yr_group=case_when(
    age_5yr_group%in%c(0,1) ~ "0-4 years", 
    age_5yr_group%in%c(5,10) ~ "5-14 years", 
    age_5yr_group%in%c(15,20) ~ "15-24 years", 
    age_5yr_group%in%c(25,30) ~ "25-34 years", 
    age_5yr_group%in%c(35,40) ~ "35-44 years", 
    age_5yr_group%in%c(45,50) ~ "45-54 years", 
    age_5yr_group%in%c(55,60) ~ "55-64 years", 
    age_5yr_group%in%c(65,70) ~ "65-74 years", 
    age_5yr_group%in%c(75,80) ~ "75-84 years", 
    age_5yr_group%in%85 ~ "85+ years", 
  )) %>% 
  group_by(age_5yr_group, metro2, race) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=(count/pop_denom)*100000) %>% 
  arrange(metro2, race, age_5yr_group) %>% 
  arrange(match(age_5yr_group, agegrp_order))


ggplot(tmp %>% filter(metro2%in%"Metropolitan") %>% mutate(race=factor(race, levels=c("NHW", "NHB", "NHAIAN", "NHAPI", "H"))), 
       aes(x=mx, y=factor(age_5yr_group, level=agegrp_order), group=age_5yr_group, fill=as.factor(race))) + 
  geom_point(aes(color=race), size=3) +
  #geom_path(aes(color=age_5yr_group)) +
  #scale_x_continuous(breaks=seq(0, 0.1, 10, 100, 1000))+
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 100, 1000))+
  scale_fill_manual(values=race_cols[2:6])+
  scale_color_manual(values=race_cols[2:6])+
  labs(x="age-specific mortality rate per 100,000",
       y="5-year age group",
       color="", fill="")+
  theme(legend.position="bottom", legend.title=element_blank())+
  isabel_theme 





results_gender %>% filter(gender%in%"Men") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_gender %>% filter(gender%in%"Men") %>% 
      select(race, metro2, le) %>% 
      spread(metro2, le) %>% 
      rename(metro_le=`Metropolitan`,
             nonmetro_le=`Non-metro`) %>% 
      full_join(
        results_gender %>% filter(gender%in%"Men") %>% 
          select(race, metro2, lv) %>% 
          spread(metro2, lv) %>% 
          rename(metro_lv=`Metropolitan`,
                 nonmetro_lv=`Non-metro`)) %>% 
      mutate(abs_le=metro_le-nonmetro_le,
             rel_le=metro_le/nonmetro_le,
             abs_lv=metro_lv-nonmetro_lv,
             rel_lv=metro_lv/nonmetro_lv)) %>%
  ungroup() %>% 
  transmute(race, 
            metro_le_ci, 
            nonmetro_le_ci,
            abs_le=as.numeric(format(abs_le, digits=2, nsmall=2)), 
            rel_le=as.numeric(format(rel_le, digits=2, nsmall=2)), 
            metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1)), 
            nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1)), 
            abs_lv=as.numeric(format(abs_lv, digits=2, nsmall=2)), 
            rel_lv=as.numeric(format(rel_lv, digits=2, nsmall=2))) %>% 
  arrange(match(race, sort_order1))


results_gender %>% filter(gender%in%"Women") %>% 
  mutate(le_ci=paste0(le=format(le, digits=1, nsmall=1), 
                      " (",
                      format(lci, digits=1, nsmall=1),
                      ", ",
                      format(uci, digits=1, nsmall=1),
                      ")")) %>% # LE (95% CI) column
  select(race, metro2, le_ci) %>%
  spread(metro2, le_ci) %>% 
  rename(metro_le_ci=`Metropolitan`,
         nonmetro_le_ci=`Non-metro`) %>% 
  full_join(
    results_gender %>% filter(gender%in%"Women") %>% 
      select(race, metro2, le) %>% 
      spread(metro2, le) %>% 
      rename(metro_le=`Metropolitan`,
             nonmetro_le=`Non-metro`) %>% 
      full_join(
        results_gender %>% filter(gender%in%"Women") %>% 
          select(race, metro2, lv) %>% 
          spread(metro2, lv) %>% 
          rename(metro_lv=`Metropolitan`,
                 nonmetro_lv=`Non-metro`)) %>% 
      mutate(abs_le=metro_le-nonmetro_le,
             rel_le=metro_le/nonmetro_le,
             abs_lv=metro_lv-nonmetro_lv,
             rel_lv=metro_lv/nonmetro_lv)) %>%
  ungroup() %>% 
  transmute(race, 
            metro_le_ci, 
            nonmetro_le_ci,
            abs_le=as.numeric(format(abs_le, digits=2, nsmall=2)), 
            rel_le=as.numeric(format(rel_le, digits=2, nsmall=2)), 
            metro_lv=as.numeric(format(metro_lv, digits=1, nsmall=1)), 
            nonmetro_lv=as.numeric(format(nonmetro_lv, digits=1, nsmall=1)), 
            abs_lv=as.numeric(format(abs_lv, digits=2, nsmall=2)), 
            rel_lv=as.numeric(format(rel_lv, digits=2, nsmall=2))) %>% 
  arrange(match(race, sort_order1))



#### 2013 vs 1990 NCHS Urban-Rural Classifications ####
load('1990_2019_nchs_90class.rdata')

dta <- master_dta_90 %>% select(-fips) %>% 
  mutate(gender=as.numeric(sex),
         age_5yr_group=as.numeric(age_5yr_group),
         year5=case_when(
           year%in%c(2015:2019) ~ as.character('2015-2019'),
           year%in%c(2010:2014) ~ as.character('2010-2014'),
           year%in%c(2005:2009) ~ as.character('2005-2009'),
           year%in%c(2000:2004) ~ as.character('2000-2004'),
           year%in%c(1995:1999) ~ as.character('1995-1999'),
           year%in%c(1990:1994) ~ as.character('1990-1994'))) %>% 
  # factoring metro2, metro6, gender, census_region to help in ggplot later
  mutate(metro2=factor(metro2, levels=c(1,0),
                       labels=c("Metropolitan", "Non-metro")),
         metro6=factor(metro6, levels=c(1:6),
                       labels=c("Large central", "Large fringe", "Medium", "Small", "Micro", "Noncore")),
         gender=factor(gender, levels=c(1,0),
                       labels=c("Men", "Women")),
         census_region=factor(census_region, levels=c(1:4),
                              labels=c("Northeast", "Midwest", "South", "West"))) %>% select(-sex)

#### results_gender_90#### 
results_gender_90 <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% filter(year5=='2015-2019') %>% 
  group_by(age_5yr_group, metro2, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro2, race, gender, age_5yr_group) %>% 
  group_by(race, metro2, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>% 
      group_by(age_5yr_group, metro2, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, race, gender, age_5yr_group) %>% 
      group_by(race, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%  
      group_by(age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%  
      group_by(age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro2, gender, age_5yr_group) %>% 
      group_by(metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )

#### results_longitudinal_90 #### 
results_longitudinal_90 <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% group_by(year5, age_5yr_group, metro2, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(year5, metro2, race, gender, age_5yr_group) %>% 
  group_by(year5, race, metro2, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, race, gender, age_5yr_group) %>% 
      group_by(year5, race, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, gender, age_5yr_group) %>% 
      group_by(year5, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% group_by(year5, age_5yr_group, metro2, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(year5, metro2, gender, age_5yr_group) %>% 
      group_by(year5, metro2, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )


#### results_urbanicity_90 #### 
results_urbanicity_90 <- 
  # MEN BY RACE/URBANICITY LE/LV
  dta %>% filter(year5=='2015-2019') %>%
  group_by(age_5yr_group, metro6, race, gender) %>% 
  summarise(count=sum(count, na.rm=T),
            pop_denom=sum(pop_denom, na.rm=T)) %>% 
  ungroup() %>% 
  filter(gender=="Men") %>% 
  # creating age-specific death rates variable, Mx 
  mutate(mx=count/pop_denom, 
         age_5yr_group=as.numeric(age_5yr_group)) %>% 
  arrange(metro6, race, gender, age_5yr_group) %>% 
  group_by(race, metro6, gender) %>% 
  group_modify(~le_lv(., age_num=0, sex="M")) %>% 
  # WOMEN BY RACE/URBANICITY LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, race, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>% 
      # creating age-specific death rates variable, Mx 
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, race, gender, age_5yr_group) %>% 
      group_by(race, metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W"))
  ) %>% 
  # MEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Men") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, gender, age_5yr_group) %>% 
      group_by(metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="B")) %>% 
      mutate(race=rep("Overall"))
  ) %>% 
  # WOMEN BY URBANICITY OVERALL LE/LV
  bind_rows(
    dta %>% filter(year5=='2015-2019') %>%
      group_by(age_5yr_group, metro6, gender) %>% 
      summarise(count=sum(count, na.rm=T),
                pop_denom=sum(pop_denom, na.rm=T)) %>% 
      ungroup() %>% 
      filter(gender=="Women") %>%
      mutate(mx=count/pop_denom, 
             age_5yr_group=as.numeric(age_5yr_group)) %>% 
      arrange(metro6, gender, age_5yr_group) %>% 
      group_by(metro6, gender) %>% 
      group_modify(~le_lv(., age_num=0, sex="W")) %>% 
      mutate(race=rep("Overall"))
  )


#### 2013 vs 1990 scatterplots #### 
# LE - METRO2  
df1 <- results_gender %>% transmute(race, 
                                    metro2, 
                                    gender, 
                                    le_2013=le)
df2 <- results_gender_90 %>% transmute(race, 
                                       metro2, 
                                       gender, 
                                       le_1990=le)
df <- df1 %>% full_join(df2, by=c("race", "metro2", "gender"))


# generating scatterplot
le_13v90 <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
                   aes(x=le_2013, y=le_1990)) +
  #geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(le_2013=seq(70, 90, by=1),
                                   le_1990=seq(70, 90, by=1)))+
  geom_point(aes(color=race), size=3) +
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  facet_grid(gender~metro2)+
  # usama's recommendation regarding points
  labs(x="Life Expectancy (1990 NCHS Classification)",
       y="Life Expectancy (2013 NCHS Classification)",
       color="", fill="", shape="")+
  guides(linetype=FALSE, color=guide_legend(nrow=1))+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
ggsave("Tables & Figures/le_13v90.pdf", le_13v90, width=18, height=15)


# generating bland altman plot
ggplot(df %>% rowwise() %>% 
         mutate(mean=mean(c(le_2013, le_1990)),
                dif=le_2013-le_1990,
                race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets,
       aes(x=mean, y=dif)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_smooth(method=lm, se=FALSE, col="black", size=0.5)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  scale_color_manual(values=race_cols)+
  facet_grid(gender~metro2)+
  labs(x="Average LE",
       y="LE Difference (2013-1990)",
       color="", fill="", shape="")+
  #guides(color=F)+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())



## LE -  METRO6
# spreading LE and 95% CIs by gender 
df1 <- results_urbanicity %>% transmute(race, 
                                        metro6, 
                                        gender, 
                                        le_2013=le)
df2 <- results_urbanicity_90 %>% transmute(race, 
                                           metro6, 
                                           gender, 
                                           le_1990=le)
df <- df1 %>% full_join(df2, by=c("race", "metro6", "gender"))


# generating scatterplot
ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
       aes(x=le_2013, y=le_1990)) +
  #geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(le_2013=seq(70, 90, by=1),
                                   le_1990=seq(70, 90, by=1)))+
  geom_point(aes(color=race), size=3) +
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  facet_grid(gender~metro6)+
  # usama's recommendation regarding points
  labs(x="2015-2019 LE (1990 NCHS Classification)",
       y="2015-2019 LE (2013 NCHS Classification)",
       color="", fill="", shape="")+
  guides(linetype=FALSE, color=guide_legend(nrow=1))+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())


# generating bland altman plot
ggplot(df %>% rowwise() %>% 
         mutate(mean=mean(c(le_2013, le_1990)),
                dif=le_2013-le_1990,
                race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets,
       aes(x=mean, y=dif)) +
  geom_hline(yintercept = 0, lty=2)+
  geom_point(aes(color=race), size=3) +
  geom_smooth(method=lm, se=FALSE, col="black", size=0.5)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  scale_color_manual(values=race_cols)+
  facet_grid(gender~metro6)+
  labs(x="Average LE",
       y="LE Difference (2013-1990)",
       color="", fill="", shape="")+
  #guides(color=F)+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())


## LV - METRO2
df1 <- results_gender %>% transmute(race, 
                                    metro2, 
                                    gender, 
                                    lv_2013=lv)
df2 <- results_gender_90 %>% transmute(race, 
                                       metro2, 
                                       gender, 
                                       lv_1990=lv)
df <- df1 %>% full_join(df2, by=c("race", "metro2", "gender"))


# generating scatterplot
lv_13v90 <- ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
                   aes(x=lv_2013, y=lv_1990, group=1)) +
  geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(lv_2013=seq(0.15, 0.35, by=0.5),
                                   lv_1990=seq(0.15, 0.35, by=0.5)))+
  geom_point(aes(color=race), size=3) +
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  facet_grid(gender~metro2)+
  # usama's recommendation regarding points
  labs(x="Lifespan Variation (1990 NCHS Classification)",
       y="Lifespan Variation (2013 NCHS Classification)",
       color="", fill="", shape="")+
  guides(linetype=FALSE, color=guide_legend(nrow=1))+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())
ggsave("Tables & Figures/lv_13v90.pdf", lv_13v90, width=18, height=15)


## LV - METRO6
df1 <- results_urbanicity %>% transmute(race, 
                                        metro6, 
                                        gender, 
                                        lv_2013=lv)
df2 <- results_urbanicity_90 %>% transmute(race, 
                                           metro6, 
                                           gender, 
                                           lv_1990=lv)
df <- df1 %>% full_join(df2, by=c("race", "metro6", "gender"))


# generating scatterplot
ggplot(df %>% mutate(race=factor(race, levels=c("Overall", "NHW", "NHB", "NHAPI", "H"))), # reorders facets, 
       aes(x=lv_2013, y=lv_1990, group=1)) +
  geom_abline(intercept = 0, slope=1, lty=1)+
  geom_line(lty=1, data=data.frame(lv_2013=seq(0.15, 0.35, by=0.5),
                                   lv_1990=seq(0.15, 0.35, by=0.5)))+
  geom_point(aes(color=race), size=3) +
  scale_color_manual(labels=c("Overall", "NH White", "NH Black", "NH Asian/Pacific Islander", "Hispanic"), values=race_cols)+
  scale_fill_manual(values=race_cols,
                    guide=guide_legend(override.aes=list(color=race_cols)))+
  facet_grid(gender~metro6)+
  # usama's recommendation regarding points
  labs(x="2015-2019 LV (1990 NCHS Classification)",
       y="2015-2019 LV (2013 NCHS Classification)",
       color="", fill="", shape="")+
  guides(linetype=FALSE, color=guide_legend(nrow=1))+
  isabel_theme+
  theme(legend.position="bottom", legend.title=element_blank())




