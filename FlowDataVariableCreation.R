## Flow Cytometry Data Variable Creation


## With this R Script I create the necessary and appropriate variables
# that I will subseqently use for the flow data analysis.


####Step 1: To begin with we will load the packages necessary for data manipulation and analysis:
library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)


####Step 2: Load in the data called "FlowData - sheet 1"

FlowData <- read_csv("data/FlowData - Sheet1.csv")

head(FlowData)



####Step 3: Creating a variable of Lipid content Fold Change

## So now, I want to create a variable that represents the fold change in lipid content between 
#light treatment and dark treatment cells

## So... how do I do this? I don't technically have paired data because my cultures are different, not
#not the same cultures measured at different time points?

##Therefore maybe I can 

#1) Calculate mean Bodipy for Light samples per group
light_means <- FlowData %>%
  filter(Light_treatment == "Light") %>%
  group_by(Strain, Temp) %>%
  summarise(mean_lightBD = mean(Bodipy_Median, na.rm = TRUE), .groups = "drop")


#2) Join the Light group means to the original dataset
BDlight_ref <- FlowData %>%
  left_join(light_means, by = c("Strain", "Temp"))


#3) Calculate fold change (only for Dark samples)
BDlight_ref <- BDlight_ref %>%
  mutate(fold_change = (Bodipy_Median / mean_lightBD))

# check if this worked
head(BDlight_ref)


### ok... Now I really need to visualise this

library(ggplot2)

ggplot(BDlight_ref, aes(x = interaction(Light_treatment, Temp), 
                        y = fold_change,
                        fill = Light_treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, aes(color = Light_treatment)) +
  facet_wrap(~ Strain, labeller = label_both) +
  labs(
    x = "Treatment (Light × Temp)",
    y = "Fluorescence Relative to Light Mean (Fold Change)",
    title = "BODIPY Fluorescence Relative to Light Mean per Strain",
    fill = "Light Treatment",
    color = "Light Treatment"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




####Step 4: Creating Lipid Concentration Variable

## So what I really need to do is standardise lipid content by size...

## Since the instrument calibration was done on the same beads for all treatment and measurements
##maybe I can start with standardising to arbitrary units, at least until I figure out this size
##calibration.

## Ok so let me start with a unit-less comparative ratio of unit lipid per unit size:

BDlight_ref <- BDlight_ref %>%
  mutate(
    lipid_per_size = Bodipy_Median / Size_Median
  )

head(BDlight_ref)




####Step 5: standardising Lipid content to a max of 1:


# Well first I need to separate my data to keep only the samples that were actually stained:

FlowDataBD <- BDlight_ref[BDlight_ref$BD == "bd",]

#cool, now I will min-max normalise the data:

FlowDataBD <- FlowDataBD %>%
  mutate( Bodipy_scaled = Bodipy_Median / max(Bodipy_Median))


# I also want to make a single tretment column, to make plotting easier
FlowDataBD <- FlowDataBD %>%
  mutate(treatment = paste(Light_treatment, Temp, sep = ""))

FlowDataBD$Strain <- as.factor(FlowDataBD$Strain)


## And make a plot to visualise the variable I just created
ggplot(FlowDataBD, aes(x = treatment, y = Bodipy_scaled, fill = Strain)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Scaled Lipid content (in BODIPY Fluorescence units, Max = 1)")


ggplot(FlowDataBD, aes(x = treatment, y = Bodipy_scaled)) +
  geom_boxplot() +
  theme_minimal() +
  ylab("Scaled BODIPY Fluorescence (Max = 1)")

## Ok, so based off this it seems like there's a clear effect of treatment,
# but also substantial and consistent variation between Strains.



####Step 6: Now I want to add in my growth rate and prop. dead variables from datasets EBData5 and GRData 5:

## So first I need to load in both my datasets

EBData5 <- read_csv("data/EBData5.csv")

GRData5 <- read.csv("data/GRData5.csv")

BDSizeData <- read.csv("data/BDSizeData.csv")


## Ok, now I isolate the 5th timepoint of the EB data and GR data

EBtp5 <- EBData5 %>% 
  filter(timepoint == 5)


GRtp5 <- GRData5 %>%
  filter(timepoint == 5)

#beautiful

### Now I need to join the grwoth rate and mortality variables onto the BDlight_ref dataset by sample ID

BDFull <- BDSizeData %>%
  left_join(EBtp5 %>% select(ID, prop_dead), by = "ID")

#and

BDFull <- BDFull %>%
  left_join(GRtp5 %>% select(ID, growthrateB), by = "ID")


##I should check if the left-join worked properly

head(BDFull)




####Step 7: Finally I need to save the BDlight_ref dataset as a csv, to make loading it and using it for analyses simpler

write.csv(BDFull, "E:/BDGM_Full.csv", row.names = FALSE)





