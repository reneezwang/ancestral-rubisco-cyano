##### Wang et al. (2023) PNAS
##### Code for Supplemental Section 4c
##### Fitting our data with other models

library(ggplot2)
library(ggpubr)

##### ==========================================================================
##### Load data from Wang et al. (2023) 

## WILL NOT RUN
## EDIT ME WITH YOUR FILE PATH NAME
allData <- readxl::read_xlsx("/Volumes/ReneeWang/PHD_ALL/Rubisco_ALL/ANCRubisco_forGithub/data.xlsx",
                             sheet = "data")

##### ==========================================================================
##### Define functions for each model:
##### 1) Sharkey & Berry (1985)
##### 2) Erez et al. (1998)
##### 3) Eichner et al. (2015)

##### Sharkey & Berry (1985) Model
sharkey_model_orig <- function(x,
                               e_rub,
                               e_db) {
  e_db+e_rub*x
}
##### Inverse Sharkey & Berry (1985) Model
sharkey_model_orig_inv <- function(e_p,
                                   e_rub,
                                   e_db) {
  (e_p-e_db)/e_rub
}

##### Erez et al. (1998) Model
erez_model <- function(x,
                       e_rub,
                       e_db,
                       a) {
  x*e_rub + a*e_db
}
##### Inverse Erez et al. (1998) Model
erez_model_inv <- function(e_p,
                       e_rub,
                       e_db,
                       a) {
  (e_p-a*e_db)/e_rub
}

##### Modified Sharkey & Berry (1985) Model from Eichner et al. (2015)
sharkey_model <- function(x,
                          e_rub,
                          a_cyt,
                          e_db) {
  x*e_rub + a_cyt*e_db
}
##### Inverse Modified Sharkey & Berry (1985) Model from Eichner et al. (2015)
sharkey_model_inv <- function(e_p,
                          e_rub,
                          a_cyt,
                          e_db) {
  (e_p-a_cyt*e_db)/e_rub
}

##### Eichner et al 2015
eichner_model <- function(x,
                          a_carb,
                          e_cyt,
                          L_carb,
                          e_rub,
                          a_cyt,
                          e_db) {
  x*(a_carb*e_cyt + L_carb*e_rub) + a_cyt*e_db
}
##### Inverse Eichner et al 2015
eichner_model_inv <- function(e_p,
                          a_carb,
                          e_cyt,
                          L_carb,
                          e_rub,
                          a_cyt,
                          e_db) {
  (e_p-a_cyt*e_db)/(a_carb*e_cyt+L_carb*e_rub)
}

##### C isotope record model; Main Text Eqn. 1
trad_model <- function(e_rub,
                       b,
                       x) {
  e_rub - b/x
}

##### Traditional model; Main text Eqn. 2
wang_trad_model <- function(e_diff,
                            e_rub,
                            x) {
  (1-x)*e_diff + x*e_rub
}


##### Calculate L_cyt from E_p based on each model =============================
## Make small dataframe for plotting
library(tidyverse)
allData.WT.allModelOutputs <- subset(allData,Strain=="WT") %>%
  group_by(Sample_ID) %>%
  mutate(
    L_cyt_sharkey_model_orig = sharkey_model_orig_inv(
      e_p = Ep_CO2_bio,
      e_rub = 25.18,
      e_db = -7.9),
    L_cyt_sharkey_acyt0 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 25.18,
      a_cyt = 0,
      e_db = -9),
    L_cyt_sharkey_acyt05 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 25.18,
      a_cyt = 0.5,
      e_db = -9),
    L_cyt_sharkey_acyt1 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 25.18,
      a_cyt = 1,
      e_db = -9),
    L_cyt_eichner_Lcarb0 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 0,
      e_rub = 25.18,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_eichner_Lcarb05 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 0.5,
      e_rub = 25.18,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_eichner_Lcarb1 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 1,
      e_rub = 25.18,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_erez_X1 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 1,
      e_db = 8,
      e_rub = 25.18),
    L_cyt_erez_X05 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 0.5,
      e_db = 8,
      e_rub = 25.18),
    L_cyt_erez_X0 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 0,
      e_db = 8,
      e_rub = 25.18)
    )

allData.ANC.allModelOutputs <- subset(allData,Strain=="ANC") %>%
  group_by(Sample_ID) %>%
  mutate(L_cyt_sharkey_model_orig = sharkey_model_orig_inv(
    e_p = Ep_CO2_bio,
    e_rub = 17.23,
    e_db = -7.9),
    L_cyt_sharkey_acyt0 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 17.23,
      a_cyt = 0,
      e_db = -9),
    L_cyt_sharkey_acyt05 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 17.23,
      a_cyt = 0.5,
      e_db = -9),
    L_cyt_sharkey_acyt1 = sharkey_model_inv(
      e_p = Ep_CO2_bio,
      e_rub = 17.23,
      a_cyt = 1,
      e_db = -9),
    L_cyt_eichner_Lcarb0 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 0.2,
      e_rub = 17.23,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_eichner_Lcarb05 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 0.5,
      e_rub = 17.23,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_eichner_Lcarb1 = eichner_model_inv(
      e_p = Ep_CO2_bio,
      a_carb = 1,
      e_cyt = 30,
      L_carb = 1,
      e_rub = 17.23,
      a_cyt = 0.8,
      e_db = -9),
    L_cyt_erez_X1 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 1,
      e_db = 8,
      e_rub = 17.23),
    L_cyt_erez_X05 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 0.5,
      e_db = 8,
      e_rub = 17.23),
    L_cyt_erez_X0 = erez_model_inv(
      e_p = Ep_CO2_bio,
      a = 0,
      e_db = 8,
      e_rub = 17.23)
    )


##### Supplemental Figure S7 ===================================================
##### Compare our 'traditional' model to the Sharkey Model & Eichner 'base' model
pPlant <- ggplot() +
  xlim(0,1) + 
  ylim(-10,25) +
  ## Original Sharkey Model
  geom_function(fun = sharkey_model_orig,
                args = list(e_rub=25,
                            e_db=-7.9),
                color = "blue",
                linetype = "solid") +
  ## Our traditional model
  geom_function(fun = wang_trad_model,
                args = list(e_rub=25,
                            e_diff=1),
                color = "green",
                linetype = "solid") +
  ## Reference lines
  geom_hline(yintercept = 25,
             color = "red") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "S7A. Plant-based model vs. Sharkey & Berry (1985) model",
       x = "Ci Leakage out of cell",
       y = "Epsilon_p (permil)") 

##### Modified Sharkey & Berry (1985) model by Eichner et al. (2015)
pEichner <- ggplot() +
  xlim(0,1) + 
  ylim(-10,25) +
  ## Modified Sharkey Model from Eichner
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=0,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=0.5,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=1,
                            e_db=-9),
                linetype = "solid") +
  ## Reference lines
  geom_hline(yintercept = 25,
             color = "red") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "S7B. Modified Sharkey & Berry model by Eichner et al. (2015)",
       x = "Ci Leakage out of cell",
       y = "Epsilon_p (permil)") 

##### Main text Eqn. 1; C Isotope Record Model
pCIsotope <- ggplot() +
  xlim(0,100) + 
  ylim(-10,25) +
  geom_function(fun = trad_model,
                args = list(e_rub=25,b=50),
                linetype = "dotted") +
  geom_function(fun = trad_model,
                args = list(e_rub=25,b=100),
                linetype = "dashed") +
  geom_function(fun = trad_model,
                args = list(e_rub=25,b=200),
                linetype = "solid") +
  geom_hline(yintercept = 25,
             color = "#3953A4") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  labs(title = "S7D. C Isotope Record Model",
       x = "[CO2(aq)] (umol/kg)",
       y = "Epsilon_p (permil)") 

##### Main text Eqn. 2; Traditional box model
pTraditional <- ggplot() +
  xlim(0,1) + 
  ylim(-10,25) +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=0,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=0.5,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25,
                            a_cyt=1,
                            e_db=-9),
                linetype = "solid") +
  geom_hline(yintercept = 25,
             color = "#3953A4") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  labs(title = "S7C. Traditional Box Model",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

ggarrange(pPlant,pEichner,
  pTraditional,pCIsotope)

##### ===========================================================================
##### Figure S8: Fitting our data using the traditional model
##### (Modified Sharkey & Berry Model from Eichner et al.)
## WT
p1 <- ggplot() +
  xlim(0,2) + 
  ylim(-10,27) +
  ## Model
  geom_function(fun = sharkey_model,
                args = list(e_rub=25.18,
                            a_cyt=0,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25.18,
                            a_cyt=0.5,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=25.18,
                            a_cyt=1,
                            e_db=-9),
                linetype = "solid") +
  ## Our WT data for a_cyt=0
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for a_cyt=0.5
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for a_cyt=1
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_vline(xintercept = 1,
             linetype = "solid") +
  geom_hline(yintercept = 25.18,
             linetype = "solid",
             color="red") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Sharkey Model WT",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

## ANC
p2 <- ggplot() +
  xlim(0,2) + 
  ylim(-10,27) +
  ## Model
  geom_function(fun = sharkey_model,
                args = list(e_rub=17.23,
                            a_cyt=0,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=17.23,
                            a_cyt=0.5,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = sharkey_model,
                args = list(e_rub=17.23,
                            a_cyt=1,
                            e_db=-9),
                linetype = "solid") +
  ## Our ANC data for a_cyt=0
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for a_cyt=0.5
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for a_cyt=1
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_sharkey_acyt1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_vline(xintercept = 1,
             linetype = "solid") +
  geom_hline(yintercept = 17.23,
             linetype = "solid",
             color="red") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Sharkey Model ANC",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

ggarrange(p1,p2,common.legend = TRUE)


##### ===========================================================================
##### Supplemental S9A: Fitting our data using Sharkey & Berry (1985)
p1 <- ggplot() +
  xlim(0,2) + 
  ylim(-10,27) +
  ## Model for WT
  geom_function(fun = sharkey_model_orig,
                args = list(e_rub=25.18,
                            e_db=-7.9),
                color = "black",
                linetype = "solid") +
  ## Our WT data
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_sharkey_model_orig,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 25.18,
             linetype = "solid",
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Sharkey Model: WT",
       x = "F3/F1",
       y = "Epsilon_p (permil)") 

p2 <- ggplot() +
  xlim(0,2) + 
  ylim(-10,27) +
  ## Model for ANC
  geom_function(fun = sharkey_model_orig,
                args = list(e_rub=17.23,
                            e_db=-7.9),
                color = "black",
                linetype = "dotted") +
  ## Our ANC data
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_sharkey_model_orig,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 17.23,
             linetype = "solid",
             color="red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Sharkey Model: ANC",
       x = "F3/F1",
       y = "Epsilon_p (permil)") 

ggarrange(p1,p2,common.legend = TRUE)

##### ===========================================================================
##### Supplemental S9B: Fitting our data with the Erez et al. model
p1 <- ggplot() +
  xlim(-0.1,1.5) + 
  ylim(0,40) +
  ## Plot function
  geom_function(fun = erez_model,
                args = list(e_rub=25.18,
                            a=1, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "dotted") +
  geom_function(fun = erez_model,
                args = list(e_rub=25.18,
                            a=0.5, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "dashed") +
  geom_function(fun = erez_model,
                args = list(e_rub=25.18,
                            a=0, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "solid") +
  ## Plot our data
  ## Our WT data for X=0
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_erez_X0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for X=0.5
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_erez_X05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for X=1
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_erez_X1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 25.18,
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Erez Model: WT",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

## ANC
p2 <- ggplot() +
  xlim(-0.1,1.5) + 
  ylim(-0,40) +
  ## Plot functions
  geom_function(fun = erez_model,
                args = list(e_rub=17.23,
                            a=1, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "dotted") +
  geom_function(fun = erez_model,
                args = list(e_rub=17.23,
                            a=0.5, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "dashed") +
  geom_function(fun = erez_model,
                args = list(e_rub=17.23,
                            a=0, ##Erez instead defines X=1 as all CO2 uptake
                            e_db=8),
                linetype = "solid") +
  ## Plot our data
  ## Our ANC data for L_carb = 0
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_erez_X0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for L_carb = 0.5
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_erez_X05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for L_carb = 1
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_erez_X1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 17.23,
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  geom_vline(xintercept = 0,
             linetype = "dotted") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Erez Model: ANC",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

ggarrange(p1,p2,common.legend=TRUE)

##### ===========================================================================
##### Supplemental S9C: Fitting our data using the Eichner et al. model
## WT
p1 <- ggplot() +
  xlim(0,1.2) + 
  ylim(-10,30) +
  ## Plot function
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=1,
                            e_rub=25.18,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=0.5,
                            e_rub=25.18,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=0,
                            e_rub=25.18,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "solid") +
  ## Plot our data
  ## Our WT data for L_carb = 0
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for L_carb = 0.5
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our WT data for L_carb = 1
  geom_point(data=allData.WT.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 25.18,
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Eichner Model: WT",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

## ANC
p2 <- ggplot() +
  xlim(0,1.2) + 
  ylim(-10,30) +
  ## Plot functions
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=1,
                            e_rub=17.23,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "dotted") +
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=0.5,
                            e_rub=17.23,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "dashed") +
  geom_function(fun = eichner_model,
                args = list(a_carb=1,
                            e_cyt=30,
                            L_carb=0,
                            e_rub=17.23,
                            a_cyt=0.8,
                            e_db=-9),
                linetype = "solid") +
  ## Plot our data
  ## Our ANC data for L_carb = 0
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb0,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for L_carb = 0.5
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb05,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Our ANC data for L_carb = 1
  geom_point(data=allData.ANC.allModelOutputs,
             aes(x=L_cyt_eichner_Lcarb1,
                 y=Ep_CO2_bio,
                 color=Condition,
                 fill=Condition),
             size = 4,
             shape = 21,
             alpha=0.7) +
  ## Reference lines
  geom_hline(yintercept = 17.23,
             color = "red") +
  geom_vline(xintercept = 1,
             linetype = "solid") +
  ## Formatting
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black", fill=NA),
        aspect.ratio = 1) +
  scale_color_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  scale_fill_manual(values=c("#4E813B", "#000000", "#3953A4")) +
  labs(title = "Eichner Model: ANC",
       x = "L_cyt",
       y = "Epsilon_p (permil)") 

ggarrange(p1,p2,common.legend=TRUE)



