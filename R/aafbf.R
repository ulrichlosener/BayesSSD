####################### CALCULATE AAFBF #####################################

# type: String that indicates the type of comparison of hypotheses. "equality" 
#       test for only equality vs. only inequality, "inequalities" test for only
#       inequality vs. only inequality.
# estimates: R object with estimates.
# n: numeric. Effective sample size.
# sigma: list with covariances.
# b: Numerical. Fraction of information to specify the prior.
# n_eff: Effective sample size.

calc_aafbf <- function(type, estimates, sigma, b, n_eff) {
    if (type == "Inequalities") {
        # Complexities
        comp1 <- .5
        comp2 <- .5
        
        # Fit
        fit2 <- pnorm(0, mean = estimates[1], sd = sqrt(sigma[[1]]))
        fit1 <- 1 - fit2
        
        # Calculation BFs
        AAFBF1u <- (fit1 / comp1)
        AAFBF2u <- (fit2 / comp2)
        AAFBF12 <- AAFBF1u / AAFBF2u
        AAFBF21 <- 1 / AAFBF12
        
        #Calculation of PMPs
        pmp1 <- AAFBF1u / (AAFBF1u + AAFBF2u)
        pmp2 <- 1 - pmp1
        output <- list(bf.12 = AAFBF12, bf.21 = AAFBF21, pmp1 = pmp1, pmp2 = pmp2)
        
    } else if (type == "Equality") {
        b_calc <- b * 1 / n_eff                   # Calculate b
        
        # Complexities
        comp0 <- dnorm(0, mean = 0, sd = sqrt(sigma[[1]] / b_calc))     # overlap of parameter under H0 and unconstrained prior -> density of the prior under Hu at the focal point 0
        comp1 <- 1 - pnorm(0, mean = 0, sd = sqrt(sigma[[2]] / b_calc))
        
        # Fit
        fit0 <- dnorm(0, mean = estimates[[1]], sd = sqrt(sigma[[1]])) # overlap of parameter under H0 and posterior -> density of the posterior at focal point 0
        fit1 <- 1 - pnorm(0, mean = estimates[[2]], sd = sqrt(sigma[[2]])) # the fit is equal to 1 - the fit of the complement
        
        # Calculation of BFs
        AAFBF0u <- fit0 / comp0                    # AAFBF of H0 vs Hu
        AAFBF1u <- fit1 / comp1                    # AAFBF of H1 vs Hu
        AAFBF01 <- AAFBF0u / AAFBF1u
        AAFBF10 <- 1 / AAFBF01
        
        # Calculation of PMPs
        pmp0 <- AAFBF0u / (AAFBF0u + AAFBF1u)
        pmp1 <- 1 - pmp0
        output <- list(bf.10 = AAFBF10, bf.01 = AAFBF01, pmp0 = pmp0, pmp1 = pmp1)
    }
    #Output
    return(output)
}
