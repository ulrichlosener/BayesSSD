#' Perform simulation-based Bayesian sample size determination
#'
#' @param eta The desired power level.
#' @param attrition The attrition pattern (FALSE for no attrition, otherwise "weibull", "modified_weibull", "linear_exponential", "log_logistic", "gompertz" or "non-parametric")
#' @param params The parameters passed to the survival function specified in "attrition". First parameter is omega, second is gamma.
#' @param m The number of datasets simulated in each iteration. The higher m, the more accurate the power level but the higher the computation cost
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
#' @param sensitivity Logical. Conduct a sensitivity analysis for the parameter fraction?
#' @param tol Tolerance for the deviation of the final result from eta. Higher values may speed up performance.
#' @param N.max The maximum sample size to be considered. Lower values may speed up performance.
#' @param N.min The minimum sample size to be considered. Higher values may speed up performance.
#' @param method The method used for hypothesis evaluation. If "bfc"/"BFc", then the hypothesis is compared against its complement via the Bayes Factor. If "bf"/"BF", then the first hypothesis is compared to the second one via the Bayes Factor. If "pmp"/"PMP", then the first hypothesis is compared to the whole set of hypotheses including the complement via posterior model probabilities.
#'
#' @return Returns the sample size (number of subjects) necessary to achieve the desired power level (eta).
#' @export
#' @examples
#' BayeSSD(eta=.8, attrition="weibull", params=c(.5,1),
#' m=100, t.points=c(0,1,2,3,4), var.u0=0.01,
#' var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .8, .8),
#' BFthres=5, fraction=1, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
#' N_max=1000)

BayeSSD <- function(eta=.8, attrition="weibull", params=c(.5,1),
                    m=100, t.points=c(0,1,2,3,4), var.u0=0.01,
                    var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .5, .8),
                    BFthres=5, fraction=1, log.grow=F, seed=NULL,
                    hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
                    N.max=1000, N.min=30, method="bfc") {

  # if function is interrupted mid-run, reset parallel processing behavior and print message
  on.exit({
    future::plan(future::sequential)
    message("The Bayesian SSD was interrupted by the user.")
  })

  # error and warning messages in case of incorrect input
  if(eta<0 | eta>1) {stop("'eta' (the desired power level) must be between 0 and 1.")}
  if(m%%1!=0 | m<1) {stop("'m' must be a positive integer.")}
  if(!is.logical(log.grow)) {stop("'log.grow' must be either TRUE or FALSE.")}
  if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE.")}
  if(any(t.points<0)) {stop("all time points must be positive.")}
  if(var.u0<0 | var.u1<0 | var.e<0) {stop("all variance components must be positive.")}
  if(BFthres<0) {stop("'BFthres' must be positive.")}
  if(fraction%%1!=0 | fraction<1) {stop("'fraction' must be a positive integer, b=fraction/N.")}
  if(m<1000) {message("Results with less than 1000 generated datasets per iteration can be unreliable.")}
  if((method=="bf" | method=="BF") & (length(hypothesis)!=2)) {stop("Method 'bf' requires exactly two hypotheses.")}

  start_time <- Sys.time()

  if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility

  N <- list()

  Nmin <- N.min         # (initial) minimal sample size
  Nmax <- N.max         # (initial) maximum sample size
  condition <- FALSE    # condition initially FALSE until power criterion is reached
  j <- 1                # iteration counter

  while(condition == F){

    N[j] <- round((Nmin + Nmax)/2 - .1, digits = 0)  # current N is the mid point between Nmin and Nmax, rounded to the lower number
    # generate data and store BFs
    results <- get_power(attrition=attrition, params=params, m=m, N=unlist(N[j]),
                               log.grow=log.grow, fraction=fraction,
                               t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                               cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                               BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

    # check if condition is met
    if(method=="bfc" | method=="BFc" | method=="bf_c" | method=="BF_c"){
      ifelse(results$power_bfc>=eta,
             Nmax <- unlist(N[j]) - 1,
             Nmin <- unlist(N[j]) + 1
      )
      pow <- results$power_bfc
    } else if(method=="pmp" | method=="PMP"){
      ifelse(results$power_pmp>=eta,
             Nmax <- unlist(N[j]) - 1,
             Nmin <- unlist(N[j]) + 1
      )
      pow <- results$power_pmp
    } else if(method=="bf" | method=="BF"){
      ifelse(results$power_bf>=eta,
             Nmax <- unlist(N[j]) - 1,
             Nmin <- unlist(N[j]) + 1
      )
      pow <- results$power_bf
    }

    # Calculate time metrics
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    avg_time_per_iter <- elapsed / j
    remaining_time <- avg_time_per_iter * (16 - j)  # max_iter = 16

    # Print progress
    cat(
      sprintf("Iteration %d: N = %d | Power = %.3f | Elapsed: %.1f minutes | Remaining: ~ %.1f minutes \n",
              j, unlist(N[[j]]), pow, elapsed, remaining_time)
    )

    # Warn about simplified models due to too little observations
    if(results$prop_simplified >= .001) {
      cat(
        sprintf("Iteration %d: %.1f%% of models required simplification (independent random effects) due to high attrition (too little observations) \n",
                j, results$prop_simplified * 100)
      )
    } else if(results$prop_simplified < .001 & results$prop_simplified > 0) {
        cat(
          "< 0.1% of models required simplification (independent random effects) due to high attrition (too little observations) \n"
        )
      }

    # if N increases by only 1 or f power level is very close to desired power level, condition is met and the algorithm stops
    if ((N[j] == Nmin+1 | Nmax == Nmin) | round(abs(pow - eta), 8) <= tol) {
      condition <- TRUE
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      cat(sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f) \n",
                  j, total_time, unlist(N[[j]]), pow))
    }

    # increase iteration by 1
    j <- j+1
  }


}


# END OF FUNCTION --------------------------------------------------------------
