## Mortality Analysis

####Step 1: To begin with we will load the packages necessary for data manipulation and analysis:

library(lme4)
library(ggplot2)
library(glmmTMB)
library(dplyr)
library(tidyverse)
install.packages("lmerTest")
library(lmerTest)


####Step 2: Next I will load in the the required dataset (here named "EBData5"):
EBData5 <- read_csv("data/EBData5.csv")
#this dataset contains all the necessary columns and data for the analysis in this R code file.

#and adjust the format of the dataset if necessary using the following code:

names(EBData5)[1] <- "ID" #to correct naming of column 1

EBData5$strain = as.factor(EBData5$strain) #to make sure Strain column is treated as a factor (Not numeric!)



####Step 3: I want to quickly visualise the data:

#the first most simple graph I can make to give me an idea what the data looks like
plot(perc_dead~day, data = EBData5)

#then I want to make a more detailed plot separating the data by my 4 treatments, I will use ggplot for this

ggplot(EBData5, aes(x = day, y = prop_dead, color = strain, group = interaction(strain, replicate))) +
  geom_line() +
  facet_grid(temperature ~ treatment) +
  labs(title = "Proportion of Dead Cells Over Time",
       x = "Time (days)",
       y = "Proportion of Dead Cells") +
  theme_minimal()

#this suggests there might be a pattern of mortality increasing over time, perhaps exacerbated by temperature



####Step 4: Statistical analysis:

#I will start with a simple linear model of my hypothesis, that proportion of dead cells is affected by time spent in
# the dark and at different temperatures, for now without any interactions

EBMod1 <- lm(perc_dead ~ strain + temperature + treatment + day, data = EBData5)
summary(EBMod1)
#and get model diagnostics:
plot(EBMod1)

#this is not bad for a preliminary model, but it does not capture the entirety of my question (i.e. the interactions),
#and I don't think that a percentage/proportion can be modeled with a linear model, as it violates assumptions of normal distribution

#I will check the distribution of the response variable
hist(EBData5$perc_dead)
# that looks like it might be a poisson, or logistic distribution

#for a poisson distribution the mean should be equal to the variance, i can check this:
mean(EBData5$perc_dead)
#23.23685

var(EBData5$perc_dead)
#435.4093

#this should not be treated as a poisson distribution then.

#another option is to treat it as a binomial regression with a logit link:

## I can analyse proportional data with a logistic regression, i.e. binomial distribution and a logit link in R
#this will bound my analysis at 0 and 1.

EBMod5 <- glm(cbind(dead,live) ~ temperature*treatment*day + strain, family = binomial(link="logit"), data = EBData5)
summary(EBMod5)

#and I want to view the model diagnostic plots:
plot(EBMod5)
#which are acceptable, except for one point with high leverage
plot(fitted(EBMod5), residuals(EBMod5))
abline(h = 0, col = "red")

qqnorm(residuals(EBMod5))
qqline(residuals(EBMod5), col = "red")



####Step 5: Interpreting coefficients from the model:

#as this is a logistic regression, I need to back transform the "estimate" to get the effect size, as an odds ratio

#Odds ratio calculation:
#Dark treatment odds ratio:
exp(0.6809679)
#=1.975789 = 197%
# in the Dark, cells are almost twice as likely to be dead than in the light 
#(by the end of the experiment)


##I will now calculate the odds ratios and confidence intervals for all coefficients:
exp(cbind(Odds_Ratio = coef(EBMod5), confint(EBMod5)))

#then I can calculate percentages sing the following equation: Percentage = (OR/(1+OR))∗100
Pdark = (1.97578915/(1 + 1.97578915))*100
Pdark
#66.39547



####Step 6: Making a graph to represent model results:

library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)


# first I will fit regression models for each combination of strain, treatment, and temperature. 
#This does not need to be the same as the statistical model, as I need this to visualise trend lines
regression_models <- EBData5 %>%
  group_by(strain, treatment, temperature) %>%
  do(model = lm(prop_dead ~ day, data = .))

# Create predictions for plotting trend lines
predictions <- regression_models %>%
  ungroup() %>%
  rowwise() %>%
  mutate(pred_data = list(
    data.frame(day = seq(min(EBData5$day), max(EBData5$day), length.out = 100),
               prop_dead = predict(model, newdata = data.frame(day = seq(min(EBData5$day), max(EBData5$day), length.out = 100)))
    )
  )) %>%
  unnest(pred_data)

# Calculate mean and standard error for plotting points
summary_data <- EBData5 %>%
  group_by(strain, treatment, temperature, day) %>%
  summarize(
    mean_prop_dead = mean(prop_dead, na.rm = TRUE),
    sd_prop_dead = sd(prop_dead, na.rm = TRUE),
    n = n(),
    se_prop_dead = sd_prop_dead / sqrt(n),
    .groups = 'drop'
  )

# Filter and plot for temperature 1C
summary_temp_1 <- filter(summary_data, temperature == 1)
predictions_temp_1 <- filter(predictions, temperature == 1)

plot_temp_05 <- ggplot() +
  geom_line(data = predictions_temp_1, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_1, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_1, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )




# Filter and plot for temperature 4
summary_temp_4 <- filter(summary_data, temperature == 4)
predictions_temp_4 <- filter(predictions, temperature == 4)

plot_temp_4 <- ggplot() +
  geom_line(data = predictions_temp_4, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_4, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_4, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Day",
    y = "Mean Proportion Dead",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  )

# Combine plots into one display
grid.arrange(plot_temp_05, plot_temp_4, ncol = 2)



### Finally I want to create a colourblind-friendly version of these plots

# First I define a colourblind-friendly palette
colorblind_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


#I also want to set the y-axis limits, to have this comparable for both graphs:
ymin <- min(
  summary_temp_1$mean_prop_dead - summary_temp_1$se_prop_dead,
  summary_temp_4$mean_prop_dead - summary_temp_4$se_prop_dead
)

ymax <- max(
  summary_temp_1$mean_prop_dead + summary_temp_1$se_prop_dead,
  summary_temp_4$mean_prop_dead + summary_temp_4$se_prop_dead
)

# Filter and plot for temperature 1C
summary_temp_1 <- filter(summary_data, temperature == 1)
predictions_temp_05 <- filter(predictions, temperature == 1)

plot_temp_1 <- ggplot() +
  geom_line(data = predictions_temp_1, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_1, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_1, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Days spent in the Dark",
    y = "Proportion Dead cells",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95")
  ) + 
  coord_cartesian(ylim = c(ymin, ymax))


# Filter and plot for temperature 4
summary_temp_4 <- filter(summary_data, temperature == 4)
predictions_temp_4 <- filter(predictions, temperature == 4)

plot_temp_4 <- ggplot() +
  geom_line(data = predictions_temp_4, aes(x = day, y = prop_dead, color = strain, linetype = treatment), size = 1.2) +
  geom_point(data = summary_temp_4, aes(x = day, y = mean_prop_dead, color = strain, shape = treatment), size = 4, stroke = 1.5) +
  geom_errorbar(data = summary_temp_4, aes(x = day, ymin = mean_prop_dead - se_prop_dead, ymax = mean_prop_dead + se_prop_dead), width = 0.3) +
  labs(
    x = "Days spent in the Dark",
    y = "Proportion Dead cells",
    color = "Strain",
    shape = "Light Treatment",
    linetype = "Light Treatment"
  ) +
  scale_color_manual(values = colorblind_palette) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95") 
  ) + 
  coord_cartesian(ylim = c(ymin, ymax))

# Combine plots into one display
grid.arrange(plot_temp_1, plot_temp_4, ncol = 2)


#1C plot: 
plot_temp_1

#4C plot:
plot_temp_4

#Lovely




