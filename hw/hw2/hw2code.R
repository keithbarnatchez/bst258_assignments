library(tidyverse)
nhefs <- read.csv('../data/nhefs.csv')
#-------------------------------------------------------------------------------
# Data cleaning

# Get set of covariates for treatment + outcome regressions
covs <- c('sex','age','race','education','smokeyrs','smokeintensity',
          'active','exercise','wt71') 

# Set outcome and treatment
outcome <- 'wt82_71' ; treatment <- 'qsmk'

# Set factor variables 
factor_vars <- c('sex','race','education',
                 'active','exercise')

# Filter out the vars we want and only keep complete cases
analysis_data <- nhefs %>%
  select(covs, outcome,treatment) %>%
  mutate_at(all_of(factor_vars), as.factor) %>%
  filter(!is.na(wt82_71))
#-------------------------------------------------------------------------------
# Weight construction

# Make formula for logistic regression (include every covariate and quadratic terms)
logistr <- as.formula(paste(treatment, '~', paste(covs, collapse = '+'),
                            '+ I(age^2) + I(wt71^2) + I(smokeintensity^2) + 
                            I(smokeyrs^2)'))

# Get IPW weights
ipw_mod <- glm(logistr, data = analysis_data, family = binomial)
ghat <- predict(ipw_mod, type = 'response')

analysis_data$ipw <- ifelse(analysis_data$qsmk == 1, 1/ghat, 1/(1-ghat))
analysis_data$stbl_ipw <- with(analysis_data,
                               ipw / mean(qsmk/ghat))

# Make df for plotting weights
plotdf <- with(analysis_data,
               data.frame(weight = c(ipw,stbl_ipw),
                                        type = rep(c('IPW','Stabilized'),
                                                   each = nrow(analysis_data))))
# Plot the weights
plotdf %>%
  ggplot(aes(x = weight)) + facet_wrap(~type, scales = 'free') +
  geom_histogram(fill='steelblue',color='black') +
  labs(title = 'Histogram of IPW Weights',
       x = 'IPW Weight',
       y = 'Frequency') +
  theme_minimal()
#-------------------------------------------------------------------------------
# ATE estimation

# Get the IPW estimate
get_ipw <- function(data, weight) {
  n <- nrow(data)
  est <- with(data, mean(qsmk * wt82_71 * weight) - mean((1 - qsmk) * wt82_71 * weight))
  se <- sqrt(with(data, var(qsmk * wt82_71 * weight - (1 - qsmk) * wt82_71 * weight) / n))
  ci <- est + c(-1, 1) * 1.96 * se
  list(est = est, se = se, ci = ci)
}

# Applying the function to both weights
ipw_est <- get_ipw(analysis_data, analysis_data$ipw)
stipw_est <- get_ipw(analysis_data, analysis_data$stbl_ipw)

comparison_df <- data.frame(
  Method = c('IPW', 'Stabilized IPW'),
  Estimate = c(ipw_est$est, stipw_est$est),
  `Standard Error` = c(ipw_est$se, stipw_est$se),
  `Lower 95% CI` = c(ipw_est$ci[1], stipw_est$ci[1]),
  `Upper 95% CI` = c(ipw_est$ci[2], stipw_est$ci[2])
)

# Using knitr to create a table
knitr::kable(comparison_df, caption = "Comparison of IPW and Stabilized IPW Estimates")
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Question 2: Doubly robust estimation

outstr <- as.formula(paste(outcome, '~', paste(covs, collapse = '+'),
                            '+ qsmk + I(age^2) + I(wt71^2) + I(smokeintensity^2) + 
                            I(smokeyrs^2) + I(qsmk*smokeintensity)'))

# fit the outcome model
outcome_mod <- lm(outstr, data = analysis_data)

# Extract predictions under each treatment
m1 <- predict(outcome_mod, newdata = analysis_data %>% mutate(qsmk = 1))
m0 <- predict(outcome_mod, newdata = analysis_data %>% mutate(qsmk = 0))

# Form the DR estimator
dr_est <- list() 
dr_est$est <- with(analysis_data,
                   ipw_est$est + 
                   mean(
                     (1-qsmk/ghat)*m1 - (1 - (1-qsmk)/(1-ghat))*m0
                   )
                   )
dr_est$se <- sqrt(with(analysis_data,
                       var(
                         (1-qsmk/ghat)*m1 - (1 - (1-qsmk)/(1-ghat))*m0
                       )/n))
dr_est$ci <- dr_est$est + c(-1,1)*1.96*dr_est$se
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Part 2: G-formula / plug-in

outstr <- as.formula(paste(outcome, '~', paste(covs, collapse = '+'),
      '+ qsmk + I(age^2) + I(wt71^2) + I(smokeintensity^2) + I(smokeyrs^2) + I(qsmk*smokeintensity)'))

# 

# Estimate outcome model
outcome_mod <- lm(outstr, data = analysis_data)

# Get the plug-in ATE estimate
m1 <- predict(outcome_mod, newdata = analysis_data %>% mutate(qsmk = 1))
m0 <- predict(outcome_mod, newdata = analysis_data %>% mutate(qsmk = 0))

meff <- marginaleffects::avg_comparisons(outcome_mod, 
                                 variables=list(qsmk= c(0,1) ))

plugin_est <- list(est = mean(m1 - m0),
                   sd = meff$std.error)
plugin_est$ci <- plugin_est$est + c(-1,1)*1.96*plugin_est$sd

