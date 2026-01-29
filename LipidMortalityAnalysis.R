### Flow cytometry Lipid Data for Mortlity Hypotheses:

####Step 1: To begin with we will load the packages necessary for data manipulation and analysis:

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)
library(lmerTest)
install.packages("lmerTest")
library(lmerTest)


####Step 2: Next I will load in the the required dataset (here named "BDGM_Full", this contains all variables created in RScript FlowDataVariableCreation):
#note, I call this BDFull in this script for simplicity

BDFull <- read.csv("data/BDGM_Full.csv")

#this dataset contains all the necessary columns and data for the analysis in this R code file.

head(BDFull)

BDFull$Strain <- as.factor(BDFull$Strain)


##And I need to subset to only Dark treatment data
BDDark <- subset(BDFull, Light_treatment == "Dark")

##OK, so now let me  create the variables that I will need:

#I already have lipid per size
#fold change (=prop retained)
#mean Bodipy in the light
#Bodipy scaled

###### H0 - would be that only size perfectly explains variation in the mortality data, 
#so let's do that:

H0 <- lm(prop_dead ~ Size_Median, data = BDDark)
summary(H0)

H0a <- lm(prop_dead ~ Size_Median*Strain*Temp, data = BDDark)
summary(H0a)
AIC(H0a)
#r2 = 61%
#-53.61882

##Well Size is not *significant* but the model R2 is quite good
plot(H0a)
#not very good residual plots


###### H1 - is that only initial lipid content explains all the variation in mortality
#the variable I need to use hear is "mean_lightBD"

H1 <- lm(prop_dead ~ mean_lightBD*Strain*Temp, data = BDDark)
summary(H1)
##this one is broken

H1a <- lm(prop_dead ~ mean_lightBD*Temp + Strain, data = BDDark)
summary(H1a)
# only Temp is important, interaction isnt

H1b <- lm(prop_dead ~ mean_lightBD*Strain + Temp, data = BDDark)
summary(H1b)
# nothing is important and I am getting errors (probably running out of df?)

H1c <- lm(prop_dead ~ mean_lightBD + Strain + Temp, data = BDDark)
summary(H1c)
AIC(H1c)
#again, only temperature is significant.
#approximately 45% of variation is explained by this model



##### H2 - final lipid content best describes variation in mortality in the dark
#The variable I want to use here is just Median Bodipy,as I am using the Dark dataset, 
#this is simply the final bodipy value

H2 <- lm(prop_dead ~ Bodipy_Median*Strain*Temp, data = BDDark)
summary(H2)
#uff... well nothing is "significant" but strain 595 interaction with Temp and also with Temp
#and bodipy is near 0.05. 
plot(H2)
#residuals are pretty bad though, some points have very high leverage

H2a <- lm(prop_dead ~ Bodipy_Median*Temp + Strain, data = BDDark)
summary(H2a)
#here strain 633, 716, temperature are significant, and Bodipy is near significance.
#R2 = 48% 
plot(H2a)
#residual plots are much better too
## This looks like a pretty good model as it is reasonable, explains a decent amount of variation
#given how limited the design was, and the model fit plots are acceptable too

##It also fits with my previous results, with a couple of strains differing, and then
#within the dark, temperature has a strong effect on mortality


H2b <- lm(prop_dead ~ Bodipy_Median*Strain + Temp, data = BDDark)
summary(H2b)
#temperature is highly important, nothing else
#r2 = 47%
plot(H2b)
#residual plots aren't as good as the previous one


##Since none of the actual interactions were significant above, let me try a model
#without interactions:

######################################################################################################
H2c <- lm(prop_dead ~ Bodipy_Median + Strain + Temp, data = BDDark)
summary(H2c)
#Ok, this looks kind of promising. Strains 633, 716, and a little bit 654 are important
#temperature is very important and bodipy a little bit but weirdly it seems like weak evidence for
#higher bodipy = higher prop dead... which doesn't make much biological sense...
#R2 = 50%
plot(H2c)
#Residual plots look good too

### Ok, model H2c is promising, let's see what H3 looks like



###### H3 - change in lipid content (or proxy of this) is the best lipid explanatory variable
#for mortality. 
#The variable to approximate change in lipid content is called "fold_change" and is calculated
#as the proportion of lipid retained between the average of the light cultures and each 
#individual dark culture (grouped by Strain and treatment)

H3 <- lm(prop_dead ~ fold_change*Strain*Temp, data = BDDark)
summary(H3)
#ok, the two things that approach significance are Strain595:Temp and fold_change:Strain595:Temp
#R2 = 47%
plot(H3)
#residuals are not great, Q-Q plot is wonky and leverage plot shows some overly significant datapoints


H3a <- lm(prop_dead ~ fold_change*Temp + Strain, data = BDDark)
summary(H3a)
#strain 716 seems to be driving the whole effect, which is possible... but seems somewhat
#unlikely as most other models indicate temperature is important. 
#R2 = 46%
plot(H3a)
#residuals are ok


H3b <- lm(prop_dead ~ fold_change*Strain + Temp, data = BDDark)
summary(H3b)
#Here Temperature is significant. The rest not
#R2 = 40%
plot(H3b)
#Residual plots are ok but I've had better. the Q-Q plot is wonky


H3c <- lm(prop_dead ~ fold_change + Strain + Temp, data = BDDark)
summary(H3c)
#Here Temperature is important and strain 716.
#R2 = 48%
AIC(H3c)
plot(H3c)
#Q-Q plots are ok, but I've seen better



###### H4 - let's see if model fit seems better if I remove lipids completely:


########################################################################################
H4 <- lm(prop_dead ~ Strain*Temp, data = BDDark)
summary(H4)
#ironically this might actually be the best fit?
#Strain 595 is significant, and the effect of temperature
#R2 = 53%
plot(H4)
#decent residual plots too

H4a <- lm(prop_dead ~ Strain+Temp, data = BDDark)
summary(H4a)
plot(H4a)




##### H5 <- Finally, I will investigate my lipid_per_size variable

H5 <- lm(prop_dead ~ lipid_per_size*Strain*Temp, data = BDDark)
summary(H5)
#R2 = 44%
plot(H5)
#not very good residuals


H5a <- lm(prop_dead ~ lipid_per_size*Temp + Strain, data = BDDark)
summary(H5a)
#R2 = 42%
plot(H5a)
#much better residuals


##########################################################################################
H5b <- lm(prop_dead ~ lipid_per_size*Strain + Temp, data = BDDark)
summary(H5b)
#R2 = 54% (so far the best)
plot(H5b)
#residuals are acceptable too
##Temperature is still the most important explanatory variable, which is reassuring
AIC(H5b)


plot(BDDark$prop_dead ~ BDDark$lipid_per_size)

#anova
aov(H5b)
anova(H5b)

### Plot to illustrate this model: 
##Problem is I can't remember what the plot that Sinead suggested was...
#I have lipid content... treatment... 
boxplot(BDDark$lipid_per_size~BDDark$treatment)



##### Making the plot to represent these models

library(ggplot2)

colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(BDDark, aes(x = lipid_per_size, 
                   y = prop_dead, 
                   color = Strain)) +
  geom_point(alpha = 0.8, size = 2, stroke = 1) + 
  
  # Trend-lines grouped by temperature
  geom_smooth(aes(group = Temp, linetype = as.factor(Temp)),
              method = "lm", se = FALSE, color = "darkgrey") +
  
  labs(x = "Lipid concentration per cell (fluorescence units)",
       y = "Proportion of dead cells",
       color = "Strain",
       linetype = "Temperature (°C)") +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal(base_size = 13)

