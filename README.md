# mvGLMM
Statistical comparisons for mean vectors generated in orientation experiments. 
In particular Generalised Linear Mixed-effects Models (GLMM), accounting for individual differences in precision.
Mean vectors are often calculated as a summary of the precision of a circular ([von Mises](https://en.wikipedia.org/wiki/Von_Mises_distribution)) distribution.
Since they are normalised to (0,1) they cannot, themselves, be normally distributed.
## Approaches
Two approaches are developed here.

 - **Log kappa**

    Mean vectors are converted to estimates of the *kappa* parameter of a von Mises distribution. 
    This falls in the range (0, infinity). 
    These values are then converted to their natural logarithm (range from -6 to 6 for mean vector lengths between 0.001 and 0.999; full range all real numbers).
    These are approximately normally distributed, allowing for statistics that rely on normality. 
    Variances may nonetheless differ, check for homogeneity of variance before applying these methods.
    
 - **Beta GLMM**
    Mean vectors are treated as samples from a [beta distribution](https://en.wikipedia.org/wiki/Beta_distribution) (common for fractions).
    Modelling is performed via [glmmTMB](https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf).
    These are [logit](https://en.wikipedia.org/wiki/Logit) transformed and effects on both the mean and dispersion
    ([concentration parameter](https://en.wikipedia.org/wiki/Beta_distribution#Mode_and_concentration) of beta distribution) of the distribution are modelled.
    Estimates and post-hoc comparisons are extracted using [emmeans](https://cran.r-project.org/web/packages/emmeans/vignettes/basics.html).
    
 For an (under-development) approach based on the von Mises distribution of the original data set, see [Tubular](https://github.com/Foztarz/Tubular).
