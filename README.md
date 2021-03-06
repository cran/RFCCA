[![Build Status](https://travis-ci.com/calakus/RFCCA.svg?branch=master)](https://travis-ci.com/calakus/RFCCA)

# RFCCA
R package which implements **R**andom **F**orest with **C**anonical **C**orrelation **A**nalysis (**RFCCA**).

**RFCCA** is a random forest method for estimating the canonical correlations between two sets of variables, *X* and *Y*, depending on 
the subject-related covariates, *Z*. The trees are built with a splitting rule specifically designed to partition the data to maximize 
the canonical correlation heterogeneity between child nodes.

For theoretical details and example data analysis, you can look at the vignette from within `R` by using the following command:

```R
vignette("RFCCA")
```

## Installation
This package is available on [CRAN](https://CRAN.R-project.org/package=RFCCA). Alternatively, you can install **RFCCA** from GitHub using the `devtools` package. Run the following code in `R` to install:

```R
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
devtools::install_github('calakus/RFCCA', build_vignettes = TRUE)
```   
## References

- Alakus, C., Larocque, D., Jacquemont, S., Barlaam, F., Martin, C.-O., Agbogba, K., Lippe, S., and Labbe, A. (2020). Conditional canonical correlation estimation based on covariates with random forests. arXiv preprint [arXiv:2011.11555](https://arxiv.org/abs/2011.11555).
