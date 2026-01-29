## Growth Rate after re-illumination Analysis

####Step 1: To begin with we will load the packages necessary for data manipulation and analysis:

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)
library(lmerTest)
install.packages("lmerTest")
library(lmerTest)


####Step 2: Next I will load in the the required dataset (here named "GRData5"):
GRData5 <- read_csv("data/GRData5.csv")
#this dataset contains all the necessary columns and data for the analysis in this R code file.

#and adjust the format of the dataset if necessary using the following code:

names(GRData5)[1] <- "ID"

GRData5$strain <- as.factor(GRData5$strain)



####Step 3: I want to quickly visualise the data:

plot(growthrateA~day, data =GRData5)

plot(growthrateB~day, data = GRData5)

#Note that growthrateA is NOT corrected for mortliaty, growthrateB IS corrected for mortality

## Now I will properly visualise the data:

## Raw Growth Rates

ggplot(GRData5, aes(x = day, y = growthrateA, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Growth rate after re-inoculation into Light",
       x = "Timepoint (days)",
       y = "Growth Rate") +
  theme_minimal()

## Corrected Growth Rates

ggplot(GRData5, aes(x = day, y = growthrateB, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Growth rate after re-inoculation into Light",
       x = "Timepoint (days)",
       y = "Growth Rate") +
  theme_minimal()



####Step 4: Statistical analysis:

#I will start with a simple linear model of my hypothesis, that growth rate declines with time spent in the dark

#An initial simple model:
GRAMod1 <- lm(growthrateA ~ strain + temperature + treatment + day, data = GRData5)
summary(GRAMod1)

#let's take a look at model diagnostics
plot(GRAMod1)
#not too bad, but could be better


## I need to check the distribution of the response variable

hist(GRData5$growthrateA)

hist(GRData5$growthrateB)
# this approaches normality, which means I can proceed with linear models


#So now I will proceed to a model that camptures the extent of my hypothesis

GRBMod4 <- lm(growthrateB ~ temperature*treatment*day + strain, data = GRData5)
summary(GRBMod4)

#I want to check the model diagnostics
plot(GRBMod4)

plot(fitted(GRBMod4), residuals(GRBMod4))
abline(h = 0, col = "red")

qqnorm(residuals(GRBMod4))
qqline(residuals(GRBMod4), col = "red")
#which are acceptable

## Deviance and AIC
deviance(GRBMod4)
#1.398032
AIC(GRBMod4)
#-733.2519


# And check the overdispersion
overdispersion <- deviance(GRBMod4) / df.residual(GRBMod4)
overdispersion
#0.004854279
#practically no overdispersion, great



####Step 5: Interpreting coefficients from the model:

# Calculating the Variation attributable to each component of the model

anova(GRBMod4)


## Prop. of total variation = Sum of squares for the variable/total sum of squares

GRanova_table <- anova(GRBMod4)

GRSS_total <- sum(GRanova_table$'Sum Sq')

# now calculate the proportion for each component:

GRanova_table$Proportion <- GRanova_table$"Sum Sq"/ GRSS_total

GRanova_table

#growth rate for strain 716 after 60 days of darkness at 1C:
GR716 = 0.2937110 + (1*0.0323467) + (-0.0204134) + (60*0.0001329) + (1*0.0027184) + (1*60*-0.0003536) + (60*-0.0024276) + (60*1*-0.0000767) + (-0.0623062)
GR716
# = 0.0825565

#growth rate for strain 590 after 90 days of darkness at 1C:
GR590 = 0.2937110 + (1*0.0323467) + (-0.0204134) + (60*0.0001329) + (1*0.0027184) + (1*60*-0.0003536) + (60*-0.0024276) + (60*1*-0.0000767)
GR590
# = 0.1448627

## Concentration calculation: Nt = N0 * e^(mu*t)

Nt716 = 100 * exp(0.0825565*7)
Nt716
# = 178.2284

Nt590 = 100 * exp(0.1448627*7)
Nt590
# = 275.6713

#fold difference:
Nt590/Nt716
# = 1.546731


##### At 4C:
#####growth rate for strain 716 after 60 days of darkness at 4C:
GR716_4C = 0.2937110 + (4*0.0323467) + (-0.0204134) + (60*0.0001329) + (4*0.0027184) + (4*60*-0.0003536) + (60*-0.0024276) + (60*4*-0.0000767) + (-0.0623062)
GR716_4C
# = 0.1102978

#growth rate for strain 590 after 90 days of darkness at 1C:
GR590_4C = 0.2937110 + (4*0.0323467) + (-0.0204134) + (60*0.0001329) + (4*0.0027184) + (4*60*-0.0003536) + (60*-0.0024276) + (60*4*-0.0000767)
GR590_4C
# = 0.172604

## Concentration calculation: Nt = N0 * e^(mu*t)

Nt716 = 100 * exp(0.0825565*7)
Nt716
# = 178.2284

Nt590 = 100 * exp(0.1448627*7)
Nt590
# = 275.6713

#fold difference:
Nt590/Nt716
# = 1.546731



####Step 6: Making a graph to represent model results:

## Load necessary libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)


# Define colorblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Calculate growth rates (based on my data)
GRData1 <- GRData5 %>%
  group_by(strain, treatment, temperature, replicate) %>%
  arrange(day) %>%
  mutate(growth_rate = growthrateB) %>%
  filter(!is.na(growth_rate))

# Fit the linear models to the growth rates
regression_models <- GRData1 %>%
  group_by(strain, treatment, temperature) %>%
  do(model = lm(growth_rate ~ day, data = .))

# Create predictions for plotting trend lines
predictions <- regression_models %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pred_data = list(
    data.frame(day = seq(min(GRData1$day), max(GRData1$day), length.out = 100),
               growth_rate = predict(model, newdata = data.frame(day = seq(min(GRData1$day), max(GRData1$day), length.out = 100)))
    )
  )) %>%
  unnest(pred_data)

# Calculate mean and standard error for plotting points
summary_data <- GRData1 %>%
  group_by(strain, treatment, temperature, day) %>%
  summarize(
    mean_growth_rate = mean(growth_rate, na.rm = TRUE),
    sd_growth_rate = sd(growth_rate, na.rm = TRUE),
    n = n(),
    se_growth_rate = sd_growth_rate / sqrt(n),
    .groups = 'drop'
  )

# Function to create plots by temperature
create_plot <- function(summary_data, predictions, temperature, title) {
  summary_filtered <- filter(summary_data, temperature == !!temperature)
  predictions_filtered <- filter(predictions, temperature == !!temperature)
  
  ggplot() +
    geom_line(data = predictions_filtered, aes(x = day, y = growth_rate, color = strain, linetype = treatment), size = 1.2) +
    geom_point(data = summary_filtered, aes(x = day, y = mean_growth_rate, color = strain, shape = treatment), size = 4, stroke = 1.5) +
    geom_errorbar(data = summary_filtered, aes(x = day, ymin = mean_growth_rate - se_growth_rate, ymax = mean_growth_rate + se_growth_rate), width = 0.3) +
    labs(
      x = "Day",
      y = "Growth Rate",
      color = "Strain",
      shape = "Light Treatment",
      linetype = "Light Treatment"
    )+ scale_color_manual(values = colorblind_palette) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 15),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    )
}

# Create plots for each temperature
plot_temp_1 <- create_plot(summary_data, predictions, 0.5, "Growth Rate at 1C Temperature")
plot_temp_4 <- create_plot(summary_data, predictions, 4, "Growth Rate at 4C Temperature")

# Combine plots into one display
grid.arrange(plot_temp_1, plot_temp_4, ncol = 2, nrow = 1)


#Yup


