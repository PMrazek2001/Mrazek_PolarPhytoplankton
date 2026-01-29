#### Lipid Growth rate Aalysis

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

#this dataset contains all the necessary columns and data for the analysis in this R code file.

head(BDFull)

BDFull$Strain <- as.factor(BDFull$Strain)


##And I need to subset to only Dark treatment data
BDDark <- subset(BDFull, Light_treatment == "Dark")



#I already have lipid per size
#fold change (=prop retained)
#mean Bodipy in the light
#Bodipy scaled



####Step 3: Testing my Hypotheses H0 - H3



###### H0 - would be that only size perfectly explains variation in the growthrate data, 
#so let's do that:

GH0 <- lm(growthrateB ~ Size_Median, data = BDDark)
summary(GH0)
##size is significant but model explains relatively little variation (37%)
plot(GH0)
#though the residual plots are alright.

##Let me try the full interaction model

GH0a <- lm(growthrateB ~ Size_Median*Strain*Temp, data = BDDark)
summary(GH0a)
##This explains a whooping 79% of variation
#the model does not seem to be "broken" either
#Size, some strains, and temperature are significant
plot(GH0a)
#generally ok, but some points have higher levarage than ideal
AIC(GH0a)

anova(GH0a)


##
####### Calculating the Variation attributable to each component of the model

anova(GH0a)


## Prop. of total variation = Sum of squares for the variable/total sum of squares

GRLanova_table <- anova(GH0a)

GRLSS_total <- sum(GRLanova_table$'Sum Sq')

# now calculate the proportion for each component:

GRLanova_table$Proportion <- GRLanova_table$"Sum Sq"/ GRLSS_total

GRLanova_table




###### H1 - is that only initial lipid content explains all the variation in growthrate
#the variable I need to use hear is "mean_lightBD"

GH1 <- lm(growthrateB ~ mean_lightBD*Strain*Temp, data = BDDark)
summary(GH1)
#this model breaks (as the mortality one)

GH1a <- lm(growthrateB ~ mean_lightBD*Strain + Temp, data = BDDark)
summary(GH1a)
#this one breaks as well

GH1b <- lm(growthrateB ~ mean_lightBD*Temp + Strain, data = BDDark)
summary(GH1b)
#This one doesn't break butthe R2 is low
#Nothing is significant here
plot(GH1b)
AIC(GH1)

##It would seem initial lipid content is not the best explanatory variable
#i.e. it is not important in determining growth rate





##### H2 - final lipid content best describes variation in growthrate in the dark
#The variable I want to use here is just Median Bodipy,as I am using the Dark dataset, 
#this is simply the final bodipy value


GH2 <- lm(growthrateB ~ Bodipy_Median*Strain*Temp, data = BDDark)
summary(GH2)
#no interactions significant, nothing significant in fact
plot(GH2)
#residual leverage isnt good either
#let me try remove the 3-way interaction
AIC(GH2)


GH2a <- lm(growthrateB ~ Bodipy_Median + Strain*Temp, data = BDDark)
summary(GH2a)
#here one strain716:temperature interaction is significant and the intercept is,
#but r2 is lower than the H0
AIC(GH2a)

GH2b <- lm(growthrateB ~ Bodipy_Median*Strain + Temp, data = BDDark)
summary(GH2b)
#nope
plot(GH2b)
#even though the residual plots look decent
AIC(GH2b)


GH2c <- lm(growthrateB ~ Bodipy_Median*Temp + Strain, data = BDDark)
summary(GH2c)
#same as before, low R2 compared to H0
plot(GH2c)
#I've had better residual plots too
AIC(GH2c)






###### H3 - change in lipid content (or proxy of this) is the best lipid explanatory variable
#for growth rate
#The variable to approximate change in lipid content is called "fold_change" and is calculated
#as the proportion of lipid retained between the average of the light cultures and each 
#individual dark culture (grouped by Strain and treatment)

GH3 <- lm(growthrateB ~ fold_change*Strain*Temp, data = BDDark)
summary(GH3)
#r2 is better, near 60%, but still lower than H0, nothing significant
AIC(GH3)


GH3a <- lm(growthrateB ~ fold_change*Strain + Temp, data = BDDark)
summary(GH3a)
#r2 at 58% several strains and fold:change interactions significant
plot(GH3a)
#residual plots are ok


GH3b <- lm(growthrateB ~ fold_change + Strain*Temp, data = BDDark)
summary(GH3b)
#this is not the one


GH3c <- lm(growthrateB ~ fold_change*Temp + Strain, data = BDDark)
summary(GH3c)
plot(GH3c)
#residuals are ok, but not this one either.





##### H4 <- Finally, I will investigate my lipid_per_size variable.
#Here I want to see if I take size into account if the model with lipids is better


GH4 <- lm(growthrateB ~ lipid_per_size*Strain*Temp, data = BDDark)
summary(GH4)
#whop that's terrible
plot(GH4)
#yeah really bad


GH4a <- lm(growthrateB ~ lipid_per_size*Temp + Strain, data = BDDark)
summary(GH4a)
#this one is not too bad
plot(GH4a)
#residual plots are decent too



####### GH5 : Let me just try add size into the H3 to see if it breaks or not

GH5 <- lm(growthrateB ~ fold_change*Strain*Temp*Size_Median, data = BDDark)
summary(GH5)
#yeah... that breaks


GH5a <- lm(growthrateB ~ fold_change*Strain*Temp + Size_Median, data = BDDark)
summary(GH5a)


GH5b <- lm(growthrateB ~ fold_change*Temp*Size_Median + Strain, data = BDDark)
summary(GH5b)
#a lot of significance here, including the 3-way interaction, but I am cautious of this model
plot(GH5b)
#well residual plots look ok
anova(GH5b)
#temp:size might not be important 
AIC(GH5b)

GH5c <- lm(growthrateB ~ fold_change*Size_Median + Strain + Temp, data = BDDark)
summary(GH5c)


#### Based on this I am inclined to go with GH0a <- lm(growthrateB ~ Size_Median*Strain*Temp, data = BDDark), as it has a pretty high
#R2 = 79%, which is at least 10% higher than any other model here.



####Step 4: Creating a graph to represent my data and model

## Now the plot I want to make for this is actually one of growth rate against size, grouping the data by strain and adding trend-lines based on size


library(ggplot2)

colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(BDDark, aes(x = Size_Median, 
                   y = growthrateB, 
                   color = Strain)) +
  geom_point(alpha = 0.8, size = 2, stroke = 1) + 
  
  # Trend-lines grouped by temperature
  geom_smooth(aes(group = Temp, linetype = as.factor(Temp)),
              method = "lm", se = FALSE, color = "darkgrey") +
  
  labs(x = "Cell Size (fluorescence units)",
       y = "Growth Rate",
       color = "Strain",
       linetype = "Temperature (°C)") +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal(base_size = 13)

