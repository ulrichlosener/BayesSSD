#===============================================================================
# Function for power calculation with Missing Data Patterns for three groups
#===============================================================================


#' Bayesian power calculation with attrition for three treatment conditions
#'
#' @param attrition The attrition pattern (FALSE for no attrition, otherwise "weibull", "modified_weibull", "linear_exponential", "log_logistic", "gompertz" or "non-parametric")
#' @param params The parameters passed to the survival function specified in "attrition". First parameter is omega, second is gamma.
#' @param m The number of datasets simulated. The higher m, the more accurate the power level but the higher the computation cost
#' @param N The number of subjects
#' @param t.points The points in time of measurement. Can be non-equidistant.
#' @param var.u0 The intercept variance.
#' @param var.u1 The slope variance.
#' @param cov The covariance between intercept variance and slope variance.
#' @param var.e The residual variance.
#' @param eff.sizes The effect sizes defined as the differences between the regression coefficients of interaction between time and condition.
#' @param fraction The fraction of information used to construct the prior for the Bayes Factor.
#' @param log.grow Use log-linear growth?#' @param var.u0
#' @param BFthres The Threshold a Bayes Factor needs to exceed in order to be considered convincing evidence.
#' @param seed Set a seed for reproducibility
#' @param hypothesis The hypothesis to be evaluated. Treatment groups are coded as "a", "b", "c", etc.
#' @param PMPthres The Threshold a Posterior Model Probability needs to exceed in order to be considered convincing evidence.
#'
#' @return Returns the power when using a Bayes Factor and the power when using Posterior Model Probabilities.
#'
#' @examples
#' getpower_mis_mv(attrition="weibull", params=c(.5,1),
#' m=100, N=100, t.points=c(0,1,2,3,4), var.u0=0.03,
#' var.u1=.1, var.e=.02, cov=0, eff.sizes=c(0, .5, .5),
#' BFthres=5, fraction=1, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9)


getpower_mis_mv <- function(attrition="weibull", params=c(.5,1),
                            m=100, N=100, t.points=c(0,1,2,3,4), var.u0=0.03,
                            var.u1=.1, var.e=.02, cov=0, eff.sizes=c(0, .5, .5),
                            BFthres=5, fraction=1, log.grow=F, seed=NULL,
                            hypothesis="a<b<c", PMPthres=.9){

  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility

  plan(multisession, workers = availableCores() - 1)  # Use all but one core

  Ns <- rep(N, m)  # object to use lapply on with first argument for the function (N)
  suppressMessages({
    bfs <- future_lapply(Ns, getbf_mis_mv, attrition=attrition,
                         params=params, hypothesis=hypothesis, t.points=t.points,
                         var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e,
                         eff.sizes=eff.sizes, fraction=fraction, log.grow=log.grow,
                         future.seed = TRUE)
  })

  plan(sequential)  # Reset plan to avoid unexpected parallel behavior later

  bf <- sapply(bfs, function(x) x[1])
  pmp <- sapply(bfs, function(x) x[2])

  power_bf <- mean(bf > BFthres)
  power_pmp <- mean(pmp > PMPthres)

  return(list(power_bf=power_bf,
              power_pmp=power_pmp))
}

# END OF FUNCTION --------------------------------------------------------------
