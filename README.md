# BayesSSD

This is the repository for the R package `BayesSSD`, which was created using R version 4.5.1. The `BayesSSD` package performs simulation-based Bayesian sample size determination for multilevel models. The authors are Ulrich Lösener and Camila Barragán Ibáñez at Utrecht University.

The function `SSD_longit` performs simulation-based Bayesian sample size determination for longitudinal experiments. Various patterns of attrition (drop-out) can be taken into account by means of a number of parametric and non-parametric survival functions. 

You can download the package from GitHub via `devtools::install_github("ulrichlosener/BayesSSD")`.

## Installation
To install the latest release version of `BayesSSD` from GitHub follow these steps:

```
install.packages("devtools")                                         # install devtools package
# For older devtools versions
devtools::install_github("ulrichlosener/BayesSSD", upgrade="never")  # install BayesSSD from GitHub
# For newer devtools versions
pak::pak("urlichlosener/BayesSSD", upgrade=FALSE)
library(BayesSSD)                                                    # load BayesSSD package
```

## Citing the package

You can cite this R-package with the following citation:

 Lösener, U. , Barragan I., C. N. (2025). BayesSSD: Bayesian Sample Size Determination for Multilevel Models. (Version 0.2.0), R package.

## Literature on BayesSSD

A number of papers have been published that refer to the method employed in this package or the package itself. Articles include the following:

-   Lösener, U., & Moerbeek, M. (2025). Bayesian sample size determination for longitudinal intervention studies with linear and log-linear growth. *Behavior Research Methods, 57*(9), 239.

-   Barragan Ibañez, C. N., Lösener, U., Moeteke, N., & Moerbeek, M. (2026). Design and analysis of behavioral intervention studies: A Bayesian approach. *Plos one, 21*(2), e0342163.

-   Ibanez, C. N. B., & Moerbeek, M. (2026). Method for sample size determination for cluster-randomized trials using the Bayes factor. *Journal of Educational and Behavioral Statistics, 51*(4), 875-907.

-   Lösener, U., & Moerbeek, M. (2026). Bayesian Sample-Size Determination for Longitudinal Trials With Attrition: The `BayesSSD` Package. *Advances in Methods and Practices in Psychological Science, 9*(3), 25152459261463494.


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
