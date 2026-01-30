#' Bayesian power calculation with attrition for three treatment conditions
#'
#' @param attrition The attrition pattern (FALSE for no attrition, otherwise "weibull", "modified_weibull", "linear_exponential", "log_logistic", "gompertz" or "non-parametric")
#' @param params The parameters passed to the survival function specified in "attrition". First parameter is omega, second is gamma.
#' @param m The number of datasets simulated. The higher m, the more accurate the power level but the higher the computation cost
#' @param N The total number of subjects or a vector with group sizes for each condition
#' @param t.points The points in time of measurement. Can be non-equidistant.
#' @param var.u0 The intercept variance.
#' @param var.u1 The slope variance.
#' @param cov The covariance between intercept variance and slope variance.
#' @param var.e The residual variance.
#' @param eff.sizes The effect sizes defined as the differences between the regression coefficients of interaction between time and condition.
#' @param fraction The fraction of information used to construct the prior for the Bayes Factor.
#' @param log.grow Use log-linear growth?
#' @param BFthres The Threshold a Bayes Factor needs to exceed in order to be considered convincing evidence.
#' @param seed Set a seed for reproducibility
#' @param hypothesis The hypothesis to be evaluated. Treatment groups are coded as "a", "b", "c", etc.
#' @param PMPthres The Threshold a Posterior Model Probability needs to exceed in order to be considered convincing evidence.
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @export
#' @return Returns the power when using a Bayes Factor and the power when using Posterior Model Probabilities.
#'
#' @examples
#' get_power(attrition="weibull", params=c(.5,1),
#' m=1000, N=100, t.points=c(0,1,2,3,4), var.u0=0.03,
#' var.u1=.1, var.e=.02, cov=0, eff.sizes=c(0, .5, .5),
#' BFthres=5, fraction=1, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9)


get_power <- function(N=100,
                      hypothesis="a<b<c",
                      eff.sizes=c(0, .5, .8),
                      t.points=c(0,1,2,3,4),
                      m=1000,
                      BFthres=5,
                      PMPthres=.9,
                      attrition="weibull",
                      params=c(.5,1),
                      var.u0=0.03,
                      var.u1=.1,
                      var.e=.02,
                      cov=0,
                      fraction=1,
                      log.grow=F,
                      seed=NULL
                      ){

  if(!is.null(seed)) {set.seed(seed)} # set user-specified seed for reproducibility


  suppressWarnings({ # suppress warning "package 'future' was built under R version 4.4.3"
  suppressMessages({ # supress messages about singular fit in MLMs

  future::plan(future::sequential)  # Reset plan to avoid unexpected leftover parallel behavior

  future::plan(future::multisession, workers = future::availableCores() - 1)  # Use all but one core

  Ns <- replicate(m, N, simplify = FALSE)  # object to use lapply on with first argument for the function (N)

  # Run simulation m times
    bfs <- future.apply::future_lapply(
      Ns,
      function(ss){
        get_bf(ss,
               attrition=attrition,
               params=params, hypothesis=hypothesis, t.points=t.points,
               var.u0=var.u0, var.u1=var.u1, cov=cov, var.e=var.e,
               eff.sizes=eff.sizes, fraction=fraction, log.grow=log.grow)
      },
      future.seed = TRUE
    )

  future::plan(future::sequential) # Reset plan to avoid unexpected parallel behavior later
  })
  })
  # extract number of simplified models due to identification issues
  prop_simplified <- mean(unlist(sapply(bfs, function(x) {x[4]})))

  # extract BFs and PMPs, handling NULL cases without error
  bfc <- vapply(bfs, function(x) {
    if (length(x[[1]]) == 0) {
      NA_real_
    } else {
      as.numeric(x[[1]])
    }
  }, numeric(1))
  pmp <- vapply(bfs, function(x) {
    if (length(x[[2]]) == 0) {
      NA_real_
    } else {
      as.numeric(x[[1]])
    }
  }, numeric(1))
  bf <- vapply(bfs, function(x) {
    if (length(x[[3]]) == 0) {
      NA_real_
    } else {
      as.numeric(x[[1]])
    }
  }, numeric(1))

  # bfc <- vapply(bfs, function(x) {as.numeric(x[[1]])}, numeric(1))
  # pmp <- vapply(bfs, function(x) {as.numeric(x[[2]])}, numeric(1))
  # bf <- vapply(bfs, function(x) {as.numeric(x[[3]])}, numeric(1))

  power_bfc <- mean(bfc>BFthres, na.rm = T)
  power_pmp <- mean(pmp>PMPthres, na.rm = T)
  power_bf <- mean(bf>BFthres, na.rm = T)

  return(list(power_bfc=power_bfc,
              power_pmp=power_pmp,
              power_bf=power_bf,
              prop_simplified=prop_simplified))
}

# END OF FUNCTION --------------------------------------------------------------
