set.seed(20211215)
# Simulate data -----------------------------------------------------------
df = data.frame( # one of each combination of species, condition and individual
  expand.grid(Species = LETTERS[1:2],
              Condition = paste(1:2),
              ID = paste0(sample(letters, 10), 1:10)
  )
)
#ensure different individuals for each species, like in real life #TODO avoid doubling
df = within(df,
            {
              ID = sapply(1:length(Species),
                          function(i)
                            {
                            switch(as.character(Species[i]), 
                                  A = paste0(ID[i], 'A'),
                                  B = paste0(ID[i], 'B')
                            )
                          }
                          )
            }
            )
#simulate data in (0,1) by averaging a binomial distribution
simulation_means <- c(A1 = 0.8,
                      A2 = 0.81,
                      B1 = 0.8,
                      B2 = 0.75)
df = within(df,
            {
              mean_vector = sapply(X = 1:length(Species), #for each species
                                   function(i) # apply the data function
                                   {
                                     switch(as.character(Species[i]), #depending on species
                                   A = switch(as.character(Condition[i]), # small condition difference
                                              '1' = rbinom(1, 1e2, simulation_means['A1'])/1e2,
                                              '2' = rbinom(1, 1e2, simulation_means['A2'])/1e2),
                                   B = switch(as.character(Condition[i]), # large condition difference
                                              '1' = rbinom(1, 1e2, simulation_means['B1'])/1e2,
                                              '2' = rbinom(1, 1e2, simulation_means['B2'])/1e2)
                                   )
                                   }
              )
              SpeciesCondition = paste0(Species, Condition) # make a compound variable
            }
            )

# Inspect data ------------------------------------------------------------
#looks about right
boxplot(mean_vector~Condition*Species, data = df)

# Perform nonparametric test ----------------------------------------------


#load the package for the Prentice test
require(muStat)
with(df, prentice.test(y = mean_vector, groups = SpeciesCondition, blocks = ID) )
# statistic: chi-square = 4.5, df = 2, p-value = 0.1054
#Prentice test is not necessarily sensitive enough


# Try parametric GLMM -----------------------------------------------------
#needs careful installing
# install.packages('TMB', type = 'source')#perform once only; don't restart
require(glmmTMB)    
#use the beta distribution for fractions
#fixed effects of condition, species, and interaction on mean
#random effects of individual on mean precision and change by condition
#N.B. each individual is only one species, so no random effects on species
#include "dispformula", allowing variance to differ by condition, species and interaction
m.beta <- glmmTMB(formula = mean_vector ~ Condition*Species + (1+Condition|ID), 
              dispformula = ~Condition*Species,
              data = df, 
              family= beta_family(link = "logit")
              )
m.beta0 <- glmmTMB(formula = mean_vector ~ (1|ID), 
              dispformula = ~1,
              data = df, 
              family= beta_family(link = "logit")
              )
anova(m.beta,m.beta0) # the model is a good fit
# Df     AIC     BIC logLik deviance  Chisq Chi Df Pr(>Chisq)   
# m.beta0  3 -125.51 -120.45 65.757  -131.51                            
# m.beta  11 -129.14 -110.56 75.569  -151.14 19.625      8    0.01185 **
summary(m.beta)
#describes random effects on mean, and fixed effects on mean and 1/variance
#N.B. these are qlogis(mean_vector) scaled values; plogis(x) recovers them
# Random effects:
#   Conditional model:
# Groups Name        Variance Std.Dev. Corr  
# ID     (Intercept) 0.002521 0.05021        
# Condition2  0.038456 0.19610  -1.00 
# Number of obs: 40, groups:  ID, 20
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          1.32494    0.08099  16.360  < 2e-16 ***
#   Condition2           0.17232    0.12794   1.347  0.17802    
# SpeciesB             0.01887    0.09760   0.193  0.84668    
# Condition2:SpeciesB -0.45991    0.15545  -2.958  0.00309 ** 
# 
# Dispersion model:
#   Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           4.5507     0.4808   9.464   <2e-16 ***
#   Condition2            0.1990     1.4131   0.141    0.888    
# SpeciesB              0.8558     0.6770   1.264    0.206    
# Condition2:SpeciesB   0.1236     3.3655   0.037    0.971    

# Perform posthocs --------------------------------------------------------
require(emmeans)
#use the package for calculating "estimated marginal means":EM means
#these are the general means across a condition, interaction or population
#which can be used to compare conditions
emm.interaction <- emmeans(object = m.beta,
                          specs = pairwise ~ Condition*Species)
emm.condition <- emmeans(object = m.beta,
                        specs = pairwise~Condition)
emm.species <- emmeans(object = m.beta,
                      specs = pairwise~Species)
summary(emm.interaction)$contrasts
# contrast  estimate     SE df t.ratio p.value
# 1 A - 2 A  -0.1723 0.1279 29 -1.347  0.5416 
# 1 A - 1 B  -0.0189 0.0976 29 -0.193  0.9974 
# 1 A - 2 B   0.2687 0.1035 29  2.596  0.0662 
# 2 A - 1 B   0.1535 0.1064 29  1.442  0.4846 
# 2 A - 2 B   0.4410 0.1082 29  4.075  0.0018 *condition 2 differs between A and B
# 1 B - 2 B   0.2876 0.0928 29  3.100  0.0210 *condition 1 differs from 2 in B
summary(emm.condition)$contrasts
# contrast estimate     SE df t.ratio p.value
#  1 - 2      0.0576 0.0803 29 0.718   0.4786  no general difference between conditions
summary(emm.species)$contrasts
# contrast estimate    SE df t.ratio p.value
#  A - B       0.211 0.0677 29 3.120   0.0041  *species do differ, excluding effects of condition

emm.interaction.phi <- emmeans(object = m.beta,
                           specs = pairwise ~ Condition*Species,
                           component = 'disp')
# contrast  estimate    SE df t.ratio p.value
# 1 A - 2 A   -0.199 1.413 29 -0.141  0.9990 
# 1 A - 1 B   -0.856 0.677 29 -1.264  0.5924 
# 1 A - 2 B   -1.178 4.248 29 -0.277  0.9924 
# 2 A - 1 B   -0.657 1.569 29 -0.419  0.9748 
# 2 A - 2 B   -0.979 3.128 29 -0.313  0.9891 In this case 2A and 2B did not differ in variance
# 1 B - 2 B   -0.323 4.427 29 -0.073  0.9999 
emm.condition.phi <- emmeans(object = m.beta,
                           specs = pairwise ~ Condition,
                           component = 'disp')
# contrast estimate   SE df t.ratio p.value
#  1 - 2      -0.261 2.82 29 -0.092  0.9270  #The conditions did not differ in variance
emm.species.phi <- emmeans(object = m.beta,
                           specs = pairwise ~ Species,
                           component = 'disp')
# contrast estimate   SE df t.ratio p.value
# A - B      -0.918 1.51 29 -0.607  0.5488  #the species did not differ in variance

# Plot the fitted model ---------------------------------------------------
#calculate predictions using a full dataset, in case we want to define individual effects
newdt = with(df,
             expand.grid(Condition = Condition,
                         Species = Species,
                         ID = ID)
              )
#model precitions for each combination
prd = predict(object = m.beta,
              newdata = newdt)
#combine with dataset of combinations for reference
pred_data = within(newdt, {mean_vector = prd})
#calculate quantiles across individual predictions
agg = aggregate(formula = mean_vector~Condition*Species,
                data = pred_data,
                FUN = quantile,
                p = c(0.025,0.5,0.975) # calculate 95%CI and median of estimates
                )

# . Open plot -------------------------------------------------------------

#Plot raw data
stripchart(mean_vector~Condition*Species,
           data = df,
           vertical = TRUE,
           method = 'jitter',
           pch = 20)
#pairs for species A
for(ind in subset(df, Species %in% 'A')$ID)
{
  with(subset(df, ID %in% ind),
       {
         lines(x = 1:2,
               y = mean_vector,
               col = adjustcolor('darkblue',
                                 alpha.f = 30/255),
               lwd = 2
               )
       }
       )
}
#pairs for species B
for(ind in subset(df, Species %in% 'B')$ID)
{
  with(subset(df, ID %in% ind),
       {
         lines(x = 2+1:2,
               y = mean_vector,
               col = adjustcolor('darkblue',
                                 alpha.f = 30/255),
               lwd = 2
               )
       }
       )
}
#add each prediction, extracted from estimated marginal means
for(i in 1:dim(summary(emm.interaction$emmeans))[1])
{
  with(summary(emm.interaction$emmeans)[i,], 
       {#95% confidence interval
        arrows(x0 = i,
               x1 = i,
               y0 = plogis(lower.CL),#fitted as logit (qlogis), convert to logistic (plogis)
               y1 = plogis(upper.CL),#fitted as logit, convert to logistic
               col = adjustcolor(col = 'orange',
                                 alpha.f = 100/255),
               code = 3,
               length = 0,
               angle = 90,
               lwd = 0.75
               )
         #standard error
        arrows(x0 = i,
               x1 = i,
               y0 = plogis(emmean - 1 *SE),
               y1 = plogis(emmean + 1 *SE),
               col = adjustcolor(col = 'darkred',
                                 alpha.f = 100/255),
               code = 3,
               length = 0.03,
               angle = 90,
               lwd = 1.5
               )
        #mean
       points(x = i,
              y = plogis(emmean),
              col = 'darkred',
              pch = 3,
              lwd = 3)
       }
      )
}
#add legend
legend(x = 'topright',
       legend = c('mean',
                  '+- standard error',
                  '95% confidence interval'),
       col = c('darkred', 'darkred', 'orange'),
       pch = c(3, NA, NA),
       lty = c(NA, 1, 1),
       cex = 0.5
)

# . Compare with original -------------------------------------------------
#for reference, compare with original means used to generate the data
for(j in 1:length(simulation_means))
{
  lines(x = j+c(-1,1)*0.2,
        y = rep(x = simulation_means[j], times = 2),
        lty = 2,
        col = 'blue')
}
#All within confidence interval (in fact within 1 standard error)
legend(x = 'top', 
       legend = 'input data',
       lty = 2,
       col = 'darkblue',
       cex = 0.5)
