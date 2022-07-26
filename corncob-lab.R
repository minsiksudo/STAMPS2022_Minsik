# Corncob lab 
## prepared by Bryan Martin, Sarah Teichman & Amy Willis 

# --------------------- Vignette Information ----------------

# We thank Dr. Thea Whitman for kindly providing us with the example data set we use 
# for this vignette. You can read more about this data in Whitman, et al. 
# "Dynamics of microbial community composition and soil organic carbon mineralization
# in soil following addition of pyrogenic and fresh organic matter." 
# The ISME Journal 10.12 (2016): 2918.

# We also use IBD microbiome data from Papa, Eliseo, et al. "Non-Invasive Mapping of 
# the Gastrointestinal Microbiota Identifies Children with Inflammatory Bowel 
# Disease." PLoS One 7(6), e39242. The data are made available by Duvallet, 
# Claire, et al. (2017). "MicrobiomeHD: the human gut microbiome in health and disease"
# [Data set]. Zenodo. We thank the authors for making their data open source and 
# easily accessible.

# --------------------- Introduction -------------------------

# Modeling microbial relative abundance is challenging for reasons such as
  
#  - different sequencing depth across samples,
#  - lots of taxa unobserved in many samples,
#  - high variability in empirical relative abundances (overdispersion),
#  - within-taxon correlation (community dynamics),
#  - multiple categorical and continuous covariates to account for in a model.

# Here, we introduce `corncob`, a microbiome relative abundance model. 
# `corncob` is able to model differential 
# abundance and differential variability, and addresses each of the challenges 
# presented above.

# Note that in order to follow along with this tutorial (but not to use `corncob`!) 
# you will need to have `phyloseq` installed. (This is installed on your cloud
# accounts).

# On the cloud we already have `corncob` installed. 
# If you're working on your own version of RStudio, install `corncob` using:
# install.packages("corncob") for the stable version
#  or devtools::install_github("bryandmartin/corncob") for the development version.

# To begin, we load our example data set as a `phyloseq` object.
#install.packages("corncob")

library(corncob)
library(phyloseq)
library(dplyr)
        
        data(soil_phylo)

# If you are unfamiliar with `phyloseq`, we can view a description of the data using:

        soil_phylo

# We now see that we have an OTU abundance table with 7770 OTUs and 119 samples. We 
# can extract using `otu_table()`. Let's examine a small subset of our data in more 
# detail.

        otu_table(soil_phylo)[1:3, 1:3]

# We can also see that we have 5 sample variables. We can extract this using 
# `sample_data()`. Let's again examine a small subset in more detail.

        sample_data(soil_phylo)[1:3, ]

# Our covariates are as follows:
  
  
#  - `Plants`: Indicator of whether plants are in the soil for this sample.
#  - `Amdmt`: Categorical variable representing one of three soil additives; none (0), 
# biochar (1), and fresh biomass (2), respectively.
#  - `ID`: Categorical variable representing different plots of soil.
#  - `Day`: Categorical variable representing one of three days of measurement; 
# day 1, day 21, and day 81, respectively.
#  - `DayAmdmt`: Categorical variable combining the `Day` and `Amdmt` variables into 
# a single variable.

# Finally, we have a taxonomy table with 7 taxonomic ranks. 

        tax_table(soil_phylo)[1:3, ]

# ---------------------------- Fitting a Model -------------------------------

# Now, let's set up our model. 

# First, let's subset our samples to only include those with the `DayAmdmt` covariate
# equal to 11 or 21 and then collapse the samples to the phylum level.
        
        soil <- soil_phylo %>% 
          phyloseq::subset_samples(DayAmdmt %in% c(11,21)) %>%
          tax_glom("Phylum") 

# Let's examine the data and the taxonomy table again.
        
        soil
        tax_table(soil)[1:5, ]

# Note that collapsing the samples is not necessary, and this model can work at any 
# taxonomic rank. However, we will later be fitting a model to every taxa. We can 
# see that by agglomerating taxa to the phylum level, we have gone from from 7770 to
# 40 taxa. Here we collapse *only* in order to decrease the runtime for the purposes of this 
# tutorial.

# Now we fit our model. We will demonstrate with Proteobacteria, or OTU.1. 

# For now, we will not include any covariates, so we use `~ 1` as our model formula 
# responses. This builds a model with only an intercept. 

corncob1 <- bbdml(formula = OTU.1 ~ 1,
                 phi.formula = ~ 1,
                 data = soil)

# Note that the first two arguments are both formulas. The first formula specifies
# the variables for the abundance part of the model and the second formula, 
# `phi.formula`, specifies the variables for the dispersion part of the model. 
# Like lm, the last argument is the data (a phyloseq object). 

# ------------------- Interpreting a Model ---------------------------------

# First, let's plot the data with our model fit on the relative abundance scale. 
# To do this, we simply type:

        plot(corncob1)

# The points represent the relative abundances. The bars represent the 95% 
# prediction intervals for the observed relative abundance by sample. 
# Computing these intervals might take ~10 seconds. 

# Next, let's color the plot by the `DayAmdmt` covariate. To do this, we add the 
# option `color = "DayAmdmt"` to our plotting code. 

        plot(corncob1, color = "DayAmdmt")

# Notice that this plot also reorders our samples so that groups appear together so 
# that they are easier to compare.

# We can observe on this plot that it might be of interest to distinguish between 
# the two groups with covariates. The average empirical relative abundance for the 
# samples with `DayAmdmt = 21` tends to be lower and less variable than the samples 
# with `DayAmdmt = 11`.

# ------------------------ Adding Covariates -----------------------------------

# Let's try modeling the expected relative abundance and the variability of the 
# absolute abundance with `DayAmdmt` as a covariate. We do this by modifying 
# `formula` and `phi.formula` as:
        
        corncob1_da <- bbdml(formula = OTU.1 ~ DayAmdmt,
                            phi.formula = ~ DayAmdmt,
                            data = soil)

# Let's also plot this data on the relative abundance scale. 

        plot(corncob1_da, color = "DayAmdmt")

# Visually, the model with covariates seems to provide a much better fit to the 
# data. 

# ------------------ Parameter Interpretation ----------------------------------

# Now that we have discussed how to fit our model, let's interpret our model output. To see a 
# summary of the model, type:

summary(corncob1_da)

# This output will look familiar from the `lm`. Covariates associated with the 
# expected relative abundance are presented separately from covariates associated 
# with the variance of the absolute abundances are preceded by.

# From this model summary, we can see that the `DayAmdmt21` abundance coefficient is
# negative and statistically significant. We estimate that this taxon has 
# a lower relative abundance for `DayAmdmt`=21 compared to `DayAmdmt`=11 -- which isn't
# surprising given what we saw in the observed relative abundances.

# We can also see that the `DayAmdmt21` dispersion coefficient is negative and 
# statistically significant. This suggests that this taxon has a different dispersion
# across `DayAmdmt`, and that samples with `DayAmdmt = 21` are expected to
# have a lower variability. This also matches what we saw from the observed abundances.

# ------------------- Analysis for Multiple Taxa --------------------------------

# What if we want to test all the taxa in our data to see if they are differentially-
# relative abundant or differentially-dispersed? We use the `differentialTest` function. It 
# will perform the above tests on all taxa, and it will control the false discovery
# rate to account for multiple comparisons.

# We specify the covariates of our model using `formula` and `phi.formula` as before, 
# except we no longer include the response term because we are testing multiple taxa. 
# We also specify which covariates we want to test for by removing them in the 
# `formula_null` and `phi.formula_null` arguments. 

# The difference between the formulas and the null version of the formulas will be 
# the variables that are tested. In this case, as when we examined the single taxon,
# we will be testing the coefficients of `DayAmdmt` for both the expected relative 
# abundance and the overdispersion.

# We set `fdr_cutoff` to be our controlled false discovery rate.

set.seed(1)

# `set.seed` is a function that will make sure that if you do the same analysis 
# involving randomness multiple times, you will get the same output. With a different
# seed you could get results that are slightly different numerically. 

da_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ 1,
                                phi.formula_null = ~ DayAmdmt,
                                test = "Wald", boot = FALSE,
                                data = soil,
                                fdr_cutoff = 0.05)

# We can see the output of the function by calling it:

da_analysis

# We can see a list of differentially-abundant taxa using: 

da_analysis$significant_taxa

# In this case, we identified 14 taxa that are differentially-abundant at a 5% FDR level 
# across `DayAmdmt` (out of the 39 taxa tested).

# We can see a list of differentially-variable taxa using: 

set.seed(1)
dv_analysis <- differentialTest(formula = ~ DayAmdmt,
                                phi.formula = ~ DayAmdmt,
                                formula_null = ~ DayAmdmt,
                                phi.formula_null = ~ 1,
                                data = soil,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
dv_analysis$significant_taxa

# We can switch the OTU labels to taxonomic labels using `otu_to_taxonomy`. We 
# supply our OTU labels as strings for the `OTU` argument. We supply the `phyloseq` 
# object for the `data` argument.

otu_to_taxonomy(OTU = da_analysis$significant_taxa, data = soil)

otu_to_taxonomy(OTU = dv_analysis$significant_taxa, data = soil)

# In this case, we identified 7 taxa that are differentially-variable across 
# `DayAmdmt` (out of the 40 taxa tested).

# We can examine a subset of the p-values of our tests using:

da_analysis$p[1:5]

# We can examine a subset of the p-values after controlling for the false discovery
# rate using:

da_analysis$p_fdr[1:5]

# where the values are now adjusted to control the false discovery rate at 0.05.
# These are q-values! (You can call them adjusted p-values, but we prefer to
# call them q-values since it's less ambiguous.)

# We can also plot the model coefficients of our results:

plot(da_analysis)

# Here, we can see that the relative abundance of Bacteria_Armatimonadetes
# is estimated to be greater for DayAmdmt=21 compared to baseline (in this case, DayAmdmt11).

# Finally, we can see a list of any taxa for which we were not able to fit a model
# using:

which(is.na(da_analysis$p)) %>% names

# In this case, we weren't able to fit `OTU.4206` automatically. It's worthwhile to
# investigate the OTU individually if this is the case. First let's check what 
# phylum this represents.

otu_to_taxonomy(OTU = "OTU.4206", data = soil)

# It may be that the model is overparameterized because there aren't enough 
# observations, or it may just be that the initializations were invalid for that 
# taxa and it needs to be re-evaluated with new initializations.

# Let's first try examining the data. 

View(otu_table(soil)["OTU.4206"])

# We see that the observed counts of OTU is zero in all samples except for `S102`, 
# where we observed a single count. Let's try fitting the model individually by 
# letting the model select the initializations automatically.

check_GN04 <- bbdml(formula = OTU.4206 ~ DayAmdmt,
                    phi.formula = ~ DayAmdmt,
                    data = soil)

# While the model fits, we should be skeptical of **any** statistical model fit on 
# a single observed count!

# If you notice any issues with `corncob`, please log them on Github at 
# https://github.com/bryandmartin/corncob/issues to help us help you!

# --------- Example of answering scientific question --------------------------

# We will now walk through several scientific questions of interest and show how 
# they can be answered using hypothesis testing with `corncob`. Note that Day and 
# Amdmt are both factor covariates with levels 0, 1, and 2. (A factor is a special
# way categories can be stored in R).

# Note that some of these are rather strange tests, and shown for demonstration of 
# the flexibility of the model only. Normally, when testing for differential 
# variability across a covariate, we recommend always controlling for the effect of 
# that covariate on the abundance. We first demonstrate examples with the soil 
# dataset.

soil_full <- soil_phylo %>% 
  tax_glom("Phylum") 

# Testing for differential abundance across Day, without controlling for anything 
# else:

ex1 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ 1,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)

plot(ex1)

# Testing for differential abundance across Day, controlling for the effect of Day 
# on dispersion:

ex2 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ Day,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Day,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex2)

# Jointly testing for differential abundance and differential variability across 
# Day:

ex3 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ Day,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex3)

# Jointly testing for differential abundance and differential variability across 
# Day, controlling for the effect of Amdmt on abundance only:

ex4 <- differentialTest(formula = ~ Day + Amdmt,
                        phi.formula = ~ Day,
                        formula_null = ~ Amdmt,
                        phi.formula_null = ~ 1,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex4)

# Jointly testing for differential abundance across Day and differential abundance 
# across Amdmt, controlling for the effect of Day and Amdmt on dispersion:

ex5 <- differentialTest(formula = ~ Day + Amdmt,
                        phi.formula = ~ Day + Amdmt,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Day + Amdmt,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex5)

# Jointly testing for differential abundance across Day, and differential dispersion
# across Amdmt, controlling for the effect of Day on Dispersion:

ex6 <- differentialTest(formula = ~ Day,
                        phi.formula = ~ Day + Amdmt,
                        formula_null = ~ 1,
                        phi.formula_null = ~ Day,
                        data = soil_full,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex6)

# We now demonstrate examples with the IBD data set. We again begin by agglomerating
# for purposes of demonstration. We agglomerate to the genus level.

data(ibd_phylo)
ibd <- ibd_phylo %>% 
  tax_glom("Genus") 

# Testing for differential abundance across IBD status, without controlling for 
# anything else:

ex7 <- differentialTest(formula = ~ ibd,
                        phi.formula = ~ 1,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = ibd,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex7)

# Because the taxonomic information is so long, we can't see our plot! 

# We can make the plot cleaner using the `level` parameter. Here we will display 
# both the family and genus information.

plot(ex7, level = c("Family", "Genus"))

# Jointly testing for differential abundance and differential variability across IBD
# status, without controlling for anything else:

ex8 <- differentialTest(formula = ~ ibd,
                        phi.formula = ~ ibd,
                        formula_null = ~ 1,
                        phi.formula_null = ~ 1,
                        data = ibd,
                        test = "Wald", boot = FALSE,
                        fdr_cutoff = 0.05)
plot(ex8, level = "Genus")




