
# looking at defining DEGs at higher cutoffs for reviewers

library(limma)
library(Glimma)
library(edgeR)

library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 <- grcm38
library(tidyverse)


## all - originally at Log2FC of 0.2, or 15%
TREATMENT_limma_results1 <- readRDS("results/TREATMENT_limma_results.RDS")

#originally at Log2FC of 0.2, or 15%
## 13 up, 43 down
TREATMENT_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC>0.2),
            Down = sum(logFC< - 0.2)) %>% 
  mutate(.,Total = Up + Down) 

# at Log2FC of 0.33, or just above 25%
## 3 up, 24 down
TREATMENT_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.33),
            Down = sum(logFC< -0.33)) %>% 
  mutate(.,Total = Up + Down) 

treatment_csv_25pct <- TREATMENT_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  filter(logFC> 0.33 | logFC< -0.33)


# at Log2FC of 0.59, or just above 50%
## 0 up, 2 down
TREATMENT_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.59),
            Down = sum(logFC< -0.59)) %>% 
  mutate(.,Total = Up + Down) 




### males
Elevel_M_limma_results1 <- readRDS("results/Elevel_M_limma_results.RDS")

# initial - 68 up, 70 down
Elevel_M_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.2),
            Down = sum(logFC< -0.2)) %>% 
  mutate(.,Total = Up + Down) 

# .33 - 18 up, 51 down
Elevel_M_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.33),
            Down = sum(logFC< -0.33)) %>% 
  mutate(.,Total = Up + Down) 

male_csv_25pct <- Elevel_M_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  filter(logFC> 0.33 | logFC< -0.33)




### e-level low
Elevel_low_limma_results1 <- readRDS("results/Elevel_low_limma_results.RDS")

# initial - 123 up, 75 down
Elevel_low_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.2),
            Down = sum(logFC< -0.2)) %>% 
  mutate(.,Total = Up + Down) 

# .33 - 24 up, 36 down
Elevel_low_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.33),
            Down = sum(logFC< -0.33)) %>% 
  mutate(.,Total = Up + Down) 

low_csv_25pct <- Elevel_low_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  filter(logFC> 0.33 | logFC< -0.33)




### e-level high
Elevel_high_limma_results1 <- readRDS("results/Elevel_high_limma_results.RDS")

# initial - 106 up, 83 down
Elevel_high_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.2),
            Down = sum(logFC< -0.2)) %>% 
  mutate(.,Total = Up + Down) 

# .33 - 34 up, 22 down
Elevel_high_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  summarise(.,Up = sum(logFC> 0.33),
            Down = sum(logFC< -0.33)) %>% 
  mutate(.,Total = Up + Down) 

high_csv_25pct <- Elevel_high_limma_results1 %>% 
  filter(., P.Value<0.05) %>% 
  filter(., P.Value != 0) %>%
  filter(logFC> 0.33 | logFC< -0.33)
