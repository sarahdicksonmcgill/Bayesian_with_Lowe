
# Package loading ---------------------------------------------------------
library(dplyr)
library(purrr)
library(stringr)

# Read in data ------------------------------------------------------------
  #For Sarah
setwd("C:/Users/sarah/OneDrive - McGill University/Bayesian With Lowe/Data")
  
  #For Sammy
#Sammy put your file path to One Drive here

peru_2000 <- readRDS("Peru/PE2000DHS_clean.rds")
peru_2004 <- readRDS("Peru/PE2004DHS_clean.rds")
peru_2007 <- readRDS("Peru/PE2007DHS_clean.rds")
peru_2009 <- readRDS("Peru/PE2009DHS_clean.rds")
peru_2010 <- readRDS("Peru/PE2010DHS_clean.rds")
peru_2011 <- readRDS("Peru/PE2011DHS_clean.rds")
peru_2012 <- readRDS("Peru/PE2012DHS_clean.rds")


# Clean and collate data for each country ---------------------------------


# Peru --------------------------------------------------------------------

keep_vars <- c("surveyID", "age", "sex", "residenceType", "education",
               "careForFamilyHIVAIDS_binary", "teacherHIVAIDS_any_binary",
               "workHIVAIDS_binary")
vars2 <- c("family_care", "teacher", "work")

peru_list <- mget(ls(pattern = "^peru_\\d{4}$"))

peru <- peru_list %>%
  map(~ .x %>%
        filter(heardHIVAIDS == "yes") %>%
        select(any_of(keep_vars))) %>%
  bind_rows(.id = "dataset") %>% 
  mutate(country = stringr::word(dataset, 1, sep = "_"),
         year = stringr::word(dataset, -1, sep = "_"),
         family_care = careForFamilyHIVAIDS_binary,
         teacher = teacherHIVAIDS_any_binary,
         work = workHIVAIDS_binary) %>% 
  dplyr::select(country, year, age, sex, residenceType, education,
                family_care, teacher, work)

miss_by_survey <- peru %>%
  group_by(country, year) %>%
  summarise(n = n(),
            across(all_of(vars2), ~ round(mean(is.na(.x)) * 100, 1)),
            .groups = "drop")

#saveRDS(peru, file = "peru.rds")


