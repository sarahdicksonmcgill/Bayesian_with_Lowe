
# Package loading ---------------------------------------------------------
#Testing 123
library(dplyr)
library(purrr)
library(stringr)
library(gtsummary)
library(ggplot2)

# Read in data ------------------------------------------------------------
  #For Sarah
setwd("C:/Users/sarah/OneDrive - McGill University/Bayesian With Lowe/Data")
  
  #For Sammy
#Sammy put your file path to One Drive here

peru_2000 <- readRDS("Peru/PE2000DHS_clean.rds")
peru_2007 <- readRDS("Peru/PE2007DHS_clean.rds")
peru_2009 <- readRDS("Peru/PE2009DHS_clean.rds")
#peru_2010 <- readRDS("Peru/PE2010DHS_clean.rds")
peru_2011 <- readRDS("Peru/PE2011DHS_clean.rds")
peru_2012 <- readRDS("Peru/PE2012DHS_clean.rds")


# Clean and collate data for each country ---------------------------------


# Peru --------------------------------------------------------------------

keep_vars <- c("surveyID", "age", "sex", "residenceType", "education",
               "careForFamilyHIVAIDS_binary", "teacherHIVAIDS_any_binary",
               "workHIVAIDS_binary", "surveyWgt_ind")
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
         work = workHIVAIDS_binary,
         any = ifelse(family_care == 1 | work == 1 | teacher == 1, 1, 0),
         all = ifelse(family_care == 1 & work == 1 & teacher == 1, 1, 0)) %>% 
  dplyr::select(country, year, age, sex, residenceType, education,
                family_care, teacher, work, any, all, surveyWgt_ind)

# miss_by_survey <- peru %>%
#   group_by(country, year) %>%
#   summarise(n = n(),
#             across(all_of(vars2), ~ round(mean(is.na(.x)) * 100, 1)),
#             .groups = "drop")

saveRDS(peru, file = "peru.rds")


# EDA/Table 1 -------------------------------------------------------------

peru %>%
  dplyr::select(-surveyWgt_ind) %>% 
  tbl_summary(
    label = list(
      age                          ~ "Age (years)",
      sex                          ~ "Sex",
      residenceType                ~ "Rurality",
      education                    ~ "Education level",
      family_care  ~ "Would not care for relative with HIV",
      teacher    ~ "Believes a teacher with HIV should stop teaching",
      work           ~ "Believes an employer could/should dismiss worker with HIV",
      any ~ "Responded yes to any stigmatizing question",
      all ~ "Responded yes to all stigmatizing questions"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    type = list(c(family_care, teacher, work) ~ "dichotomous"),
    value = list(c(family_care, teacher, work) ~ 1),
    missing = "ifany",
    missing_text = "Missing",
    missing_stat = "{N_miss} ({p_miss}%)"
  ) %>%
  modify_header(label = "**Characteristic**",
                stat_0 = "**N = {style_number(N)}**") %>%
  bold_labels()

stigma_vars <- c("family_care", "teacher", "work", "any", "all")

prev <- peru %>%
  pivot_longer(all_of(stigma_vars), names_to = "question", values_to = "response") %>%
  filter(!is.na(response)) %>%
  group_by(question, year) %>%
  summarise(n = n(),
            n_yes = sum(response == 1),
            .groups = "drop") %>%
  mutate(prev  = n_yes / n,
         se    = sqrt(prev * (1 - prev) / n),
         lower = pmax(0, prev - 1.96 * se),
         upper = pmin(1, prev + 1.96 * se),
         across(c(prev, lower, upper), ~ round(.x * 100, 1))) %>%
  select(-se)

prev %>% 
  filter(question %in% c("work", "teacher", "family_care")) %>% 
  ggplot(aes(x = year, y = prev, colour = question, group = question)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  labs(title = "Prevelance of stigmatizing attitudes in Peru, by year",
       x = "Year",
       y = "Prevelance of stigma (%)") +
  scale_colour_manual(
    name   = "Question",
    values = c(family_care = "#0072B2", teacher = "#E69F00", work = "#009E73"),
    labels = c(family_care = "Would not care for relative",
               teacher     = "Teacher should not teach",
               work        = "Employer could dismiss")
  )


