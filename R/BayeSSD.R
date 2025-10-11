#' Perform simulation-based Bayesian sample size determination
#'
#' @param eta The desired power level.
#' @param attrition The attrition pattern (`FALSE` for no attrition, otherwise `weibull`, `modified_weibull`, `linear_exponential`, `log_logistic`, `gompertz` or `non-parametric`)
#' @param params The parameters passed to the survival function specified in `attrition`. First parameter is omega, second gamma and third (if applicable) is kappa.
#' @param m The number of datasets simulated in each iteration. The higher `m`, the more accurate the power level but the higher the computation time
#' @param t.points The points in time of measurement. Can be non-equidistant.
#' @param var.u0 The intercept variance.
#' @param var.u1 The slope variance.
#' @param cov The covariance between intercept variance and slope variance.
#' @param var.e The residual variance.
#' @param eff.sizes The effect sizes defined as the differences between the regression coefficients of interaction between time and condition.
#' @param log.grow Logical. Use log-linear growth?
#' @param BFthres The Threshold a Bayes Factor needs to exceed in order to be considered convincing evidence.
#' @param seed Set a seed for reproducibility?
#' @param hypothesis List of the hypotheses to be evaluated. Treatment groups are coded as `a`, `b`, `c`, etc.
#' @param PMPthres The Threshold a Posterior Model Probability needs to exceed in order to be considered convincing evidence.
#' @param sensitivity Logical. Conduct a sensitivity analysis for the parameter `fraction`?
#' @param tol Tolerance for the deviation of the final result from `eta`. Higher values may speed up performance.
#' @param N.max The maximum sample size to be considered. Lower values may speed up performance.
#' @param N.min The minimum sample size to be considered. Higher values may speed up performance.
#' @param method The method used for hypothesis evaluation. If `bfc`/`BFc`, then the hypothesis is compared against its complement via the Bayes Factor. If `bf`/`BF`, then the first hypothesis is compared to the second one via the Bayes Factor. If `pmp`/`PMP`, then the first hypothesis is compared to the whole set of hypotheses including the complement via posterior model probabilities.
#'
#' @return Returns the sample size (number of subjects) necessary to achieve the desired power level `eta`.
#' @export
#' @examples
#' BayeSSD(eta=.8, attrition="weibull", params=c(.5,1),
#' m=100, t.points=c(0,1,2,3,4), var.u0=0.01,
#' var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .8, .8),
#' BFthres=5, fraction=1, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
#' N_max=1000)

BayeSSD <- function(eta=.8, attrition="weibull", params=c(.5,1),
                    m=10000, t.points=c(0,1,2,3,4), var.u0=0.01,
                    var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .5, .8),
                    BFthres=5, log.grow=F, seed=NULL,
                    hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.01,
                    N.max=1000, N.min=30, method="bfc") {

  # Use calling handlers to catch interrupts
  withCallingHandlers({

    # error and warning messages in case of incorrect input
    if(eta<0 | eta>1) {stop("'eta' (the desired power level) must be between 0 and 1.")}
    if(m%%1!=0 | m<1) {stop("'m' must be a positive integer.")}
    if(!is.logical(log.grow)) {stop("'log.grow' must be either TRUE or FALSE.")}
    if(is.logical(sensitivity)==F) {stop("'sensitivity' must be either TRUE or FALSE.")}
    if(any(t.points<0)) {stop("all time points must be positive.")}
    if(var.u0<0 | var.u1<0 | var.e<0) {stop("all variance components must be positive.")}
    if(BFthres<0) {stop("'BFthres' must be positive.")}
    if(m<1000) {message("Results with less than 1000 generated datasets per iteration can be unreliable.")}
    if((method=="bf" | method=="BF") & (length(hypothesis)!=2)) {stop("Method 'bf' requires exactly two hypotheses.")}

    start_time <- Sys.time()

    if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility

    N <- list()
    n_cond <- length(eff.sizes) # extract number of conditions
    candidate_N <- seq(from=N.min, to=N.max, by=n_cond) # set of candidate N (divisible by number of conditions)

    condition <- FALSE                               # condition initially FALSE until power criterion is reached
    j <- 1                                           # iteration counter
    pow <- 0                                         # initialize power
    av_it <- round(log((N.max - N.min + 1), base=2)) # approximation of average numbers of iterations
    N_min <- N.min
    N_max <- N.max

    if(sensitivity==F){
      ################### without sensitivity analysis #########################
      while(condition == F){

        N_mid <- round((N_min + N_max)/2 - .1, digits = 0)  # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
        N[[j]] <- candidate_N[which.min(abs(candidate_N - N_mid))]   # find the nearest candidate value for N to N_mid

        # set m according to iteration/difference between actual (pow) and desired power (eta)
        if(m>=5000){
          if(j==1 | abs(pow-eta) > .1){ # in the first iteration or if the difference between pow and eta is at least .1, set m to 1000
            current_m <- 1000
          } else {
            current_m <- m
          }
        }

        # generate data and store BFs
        results <- get_power(attrition=attrition, params=params, m=current_m, N=unlist(N[[j]]),
                             log.grow=log.grow, fraction=1,
                             t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                             cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                             BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

        # check if condition is met
        if(method=="bfc" | method=="BFc" | method=="bf_c" | method=="BF_c"){
          if(results$power_bfc>=eta){
            N_max <- unlist(N[[j]]) - n_cond
          } else {
            N_min <- unlist(N[[j]]) + n_cond
          }
          pow <- results$power_bfc
        } else if(method=="pmp" | method=="PMP"){
          if(results$power_pmp>=eta){
            N_max <- unlist(N[[j]]) - n_cond
          } else {
            N_min <- unlist(N[[j]]) + n_cond
          }
          pow <- results$power_pmp
        } else if(method=="bf" | method=="BF"){
          if(results$power_bf>=eta){
            N_max <- unlist(N[[j]]) - n_cond
          } else{
            N_min <- unlist(N[[j]]) + n_cond
          }
          pow <- results$power_bf
        }

        # Calculate time metrics
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        avg_time_per_iter <- elapsed / j
        remaining_time <- avg_time_per_iter * (av_it - j)

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

        # if power level is close enough to desired power level, condition is met and the algorithm stops
        if (round(abs(pow - eta), 8) <= tol) {
          condition <- TRUE
          total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          cat(sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f) \n",
                      j, total_time, unlist(N[[j]]), pow))
        }
        # increase iteration by 1
        j <- j+1
      }

    } else {
      ##################### with sensitivity analysis ##########################
      for (i in 1:3){
        # Re-initialize
        N <- list()
        condition <- FALSE    # condition initially FALSE until power criterion is reached
        j <- 1                # iteration counter
        pow <- 0              # initialize power
        av_it <- round(log((N.max - N.min + 1), base=2)) # approximation of average numbers of iterations
        N_min <- N.min
        N_max <- N.max
        # print info on sensitivity analysis
        cat("\n", "\n", "Sensitivity analysis for fraction = ", i, "\n")
        while(condition == F){

          N_mid <- round((N_min + N_max)/2 - .1, digits = 0)  # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
          N[[j]] <- candidate_N[which.min(abs(candidate_N - N_mid))]   # find the nearest candidate value for N to N_mid

          # set m according to iteration/difference between actual (pow) and desired power (eta)
          if(m>=5000){
            if(j==1 | abs(pow-eta) > .1){ # in the first iteration or if the difference between pow and eta is at least .1, set m to 1000
              current_m <- 1000
            } else {
              current_m <- m
            }
          }

          # generate data and store BFs
          results <- get_power(attrition=attrition, params=params, m=current_m, N=unlist(N[[j]]),
                               log.grow=log.grow, fraction=i,
                               t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                               cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                               BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

          # check if condition is met
          if(method=="bfc" | method=="BFc" | method=="bf_c" | method=="BF_c"){
            if(results$power_bfc>=eta){
              N_max <- unlist(N[[j]]) - n_cond
            } else {
              N_min <- unlist(N[[j]]) + n_cond
            }
            pow <- results$power_bfc
          } else if(method=="pmp" | method=="PMP"){
            if(results$power_pmp>=eta){
              N_max <- unlist(N[[j]]) - n_cond
            } else {
              N_min <- unlist(N[[j]]) + n_cond
            }
            pow <- results$power_pmp
          } else if(method=="bf" | method=="BF"){
            if(results$power_bf>=eta){
              N_max <- unlist(N[[j]]) - n_cond
            } else{
              N_min <- unlist(N[[j]]) + n_cond
            }
            pow <- results$power_bf
          }

          # Calculate time metrics
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          avg_time_per_iter <- elapsed / j
          remaining_time <- avg_time_per_iter * (av_it - j)

          # Print progress
          cat(
            sprintf("Iteration %d: N = %d | Power = %.3f | Elapsed: %.1f minutes | Total remaining: ~ %.1f minutes \n",
                    j, unlist(N[[j]]), pow, elapsed, remaining_time)
          )

          # Warn about simplified models due to too little observations
          if(results$prop_simplified >= .001) {
            cat(
              sprintf("%.1f%% of models required simplification (independent random effects) due to high attrition (too little observations) \n",
                      results$prop_simplified * 100)
            )
          } else if(results$prop_simplified < .001 & results$prop_simplified > 0) {
            cat(
              "< 0.1% of models required simplification (independent random effects) due to high attrition (too little observations) \n"
            )
          }

          # if power level is close enough to desired power level, condition is met and the algorithm stops
          if (round(abs(pow - eta), 8) <= tol) {
            condition <- TRUE
            total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
            cat(sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f) \n",
                        j, total_time, unlist(N[[j]]), pow))
          }

          # increase iteration by 1
          j <- j+1
        }
      }
    }

    # in case of interruption or error, reset parallel behavior
  }, interrupt = function(e) {
    message("Interrupt detected - resetting parallel processing plan")
    future::plan(future::sequential)
  }, error = function(e) {
    message("Error detected - resetting parallel processing plan")
    future::plan(future::sequential)
    stop(e)  # Re-throw the error
  })
}

# END OF FUNCTION --------------------------------------------------------------
