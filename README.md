# BayesSSD

This is the repository for the R package `BayesSSD`, which was created using R version 4.5.1. The function `SSD_longit` performs simulation-based Bayesian sample size determination for longitudinal experiments. Various patterns of attrition (drop out) can be taken into account by means of a number of parametric and non-parametric survival functions. 
You can download the package from GitHub via `devtools::install_github("ulrichlosener/BayesSSD")`.

The accompanying paper for the `BayesSSD` package authored by Ulrich Lösener and Mirjam Moerbeek is called "Bayesian Sample Size Determination for Longitudinal Trials with Attrition - the BayesSSD Package" and has been submitted to *Advances in Methods and Practices in Psychological Science*.

## Installation
To install the latest release version of `BayesSSD` from GitHub follow these steps:

```
install.packages("devtools")                                         # install devtools package

devtools::install_github("ulrichlosener/BayesSSD", upgrade="never")  # install BayesSSD from GitHub

library(BayesSSD)                                                    # load BayesSSD package
```

## Citing BayesSSD

You can cite this R-package with the following citation:

 Lösener, U. (2025). BayesSSD: Bayesian Sample Size Determination for Longitudinal Experiments with Attrition. (Version 0.1.0), R package.


## Contributing and Contact Information

If you have ideas, please get involved. You can contribute by opening an
issue on GitHub, or sending a pull request with proposed features.
Contributions in code must adhere to the [tidyverse style
guide](https://style.tidyverse.org/).

-   File a GitHub issue [here](https://github.com/ulrichlosener/BayesSSD)
-   Make a pull request [here](https://github.com/ulrichlosener/BayesSSD/pulls)

By participating in this project, you agree to abide by the [Contributor
Code of Conduct v2.0](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html).

If you have questions or comments, feel free to get in contact via [e-mail](mailto:u.c.losener1@uu.nl).
