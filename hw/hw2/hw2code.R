library(tidyverse)
nhefs <- read.csv('../data/nhefs.csv')
#-------------------------------------------------------------------------------
# Data cleaning

# Get set of covariates for treatment + outcome regressions
covs <- c('sex','age','race','education','smokeyrs','smokeintensity',
          'active','exercise','wt71') 

# Set outcome and treatment
outcome <- 'wt82_71'
treatment <- 'qsmk'

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
n <- nrow(analysis_data)
ipw_est <- list(mean = with(analysis_data,
                mean(qsmk*wt82_71*ipw) - 
                  mean((1-qsmk)*wt82_71*ipw)),
                se = sqrt(with(analysis_data,
                               sd(qsmk*wt82_71*ipw - (1-qsmk)*wt82_71*ipw)/sqrt(n))))

# Get a 95% CI
ipw_est$ci <- ipw_est$mean + c(-1,1)*1.96*ipw_est$se
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Question 2: Doubly robust estimation

outstr <- as.formula(paste(outcome, '~', paste(covs, collapse = '+'),
                            '+ I(age^2) + I(wt71^2) + I(smokeintensity^2) + 
                            I(smokeyrs^2)'))

# fit the outcome model
outcome_mod <- lm(outstr, data = analysis_data)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Part 2: G-formula / plug-in

outstr <- as.formula(paste(outcome, '~', paste(covs, collapse = '+'),
                            '+ I(age^2) + I(wt71^2) + I(smokeintensity^2) + 
                            I(smokeyrs^2) + qsmk:smokeintensity'))

# Estimate outcome model
outcome_mod <- lm(outstr, data = analysis_data)
