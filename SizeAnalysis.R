#### Flow Cytometry Data Size Variable Exploration and Analysis

####Step 1: To begin with we will load the packages necessary for data manipulation and analysis:

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)
library(lmerTest)

####Step 2: Next I will load in the the required dataset (here named "BDGM_Full", this contains all variables created in RScript FlowDataVariableCreation):
#note, I call this BDFull in this script for simplicity

BDFull <- read.csv("data/BDGM_Full.csv")

BDFull$Strain <- as.factor(BDFull$Strain)



####Step 3: Exploring Data Variables SSC and FSC

##Now, I have SSC, let's take a look at that

boxplot(BDFull$SSC_Med ~ BDFull$treatment_combo)

#Ok, but let me colour the data by strain, get a better idea of SSC distribution

# Boxplot with treatment on x-axis, strain as color (grouped)
ggplot(BDFull, aes(x = treatment_combo, y = SSC_Med, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(title = "SSC by Treatment Group",
       x = "Treatment",
       y = "SSC") +
  theme_minimal()

####

##ok, Interestingly it seems that SSC is higher in the dark treatment indicating higher granularity/internal complexity


## Now how about I compare fsc (just called size) and ssc
#a 

ggplot(BDFull, aes(x = Size_Median, y = SSC_Med, color = Strain)) +
  geom_point(alpha = 0.7) +
  labs(title = "SSC vs FSC by Strain",
       x = "FSC (Size)",
       y = "SSC") +
  theme_minimal()

##ok, maybe it would help to separate the plots by treatment

ggplot(BDFull, aes(x = Size_Median, y = SSC_Med, color = Strain)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ treatment_combo) +
  labs(title = "SSC vs FSC by Strain & Treatment",
       x = "FSC",
       y = "SSC") +
  theme_minimal()


BDFull$Strain <- as.factor(BDFull$Strain)


ggplot(BDFull, aes(x = treatment_combo, y = Size_Median, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(x = "Treatment",
       y = "Median Cell Size (FSC)",
       legend.title = "Strain") +
  theme_minimal()




