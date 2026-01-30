#' Perform simulation-based Bayesian sample size determination
#'
#' @param eta The desired power level.
#' @param hypothesis List of the hypotheses to be evaluated. Treatment groups are coded as `a`, `b`, `c`, etc.
#' @param eff.sizes The effect sizes defined as the differences between the regression coefficients of interaction between time and condition.
#' @param BFthres The Threshold a Bayes Factor needs to exceed in order to be considered convincing evidence.
#' @param PMPthres The Threshold a Posterior Model Probability needs to exceed in order to be considered convincing evidence.
#' @param method The method used for hypothesis evaluation. If `bfc`/`BFc`, then the hypothesis is compared against its complement via the Bayes Factor. If `bf`/`BF`, then the first hypothesis is compared to the second one via the Bayes Factor. If `pmp`/`PMP`, then the first hypothesis is compared to the whole set of hypotheses including the complement via posterior model probabilities.
#' @param attrition The attrition pattern (`FALSE` for no attrition, otherwise `weibull`, `modified_weibull`, `linear_exponential`, `log_logistic`, `gompertz` or `non-parametric`)
#' @param params The parameters passed to the survival function specified in `attrition`. First parameter is omega, second gamma and third (if applicable) is kappa.
#' @param t.points The points in time of measurement. Can be non-equidistant.
#' @param var.u0 The intercept variance.
#' @param var.u1 The slope variance.
#' @param var.e The residual variance.
#' @param cov The covariance between intercept variance and slope variance.
#' @param m The number of datasets simulated in each iteration. The higher `m`, the more accurate the power level but the higher the computation time
#' @param log.grow Logical. Use log-linear growth?
#' @param seed Set a seed for reproducibility?
#' @param sensitivity Logical. Conduct a sensitivity analysis for the b fraction?
#' @param tol Tolerance for the deviation of the final result from `eta`. Higher values may speed up performance.
#' @param N.max The maximum sample size to be considered. Lower values may speed up performance.
#' @param N.min The minimum sample size to be considered. Higher values may speed up performance.
#'
#' @return Returns the sample size (number of subjects) necessary to achieve the desired power level `eta`.
#' @export
#' @examples
#' SSD_longit(eta=.8, attrition="weibull", params=c(.5,1),
#' m=100, t.points=c(0,1,2,3,4), var.u0=0.01,
#' var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .8, .8),
#' BFthres=5, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
#' N.max=1000)

SSD_longit <- function(eta=.8, hypothesis="a<b<c", eff.sizes=c(0, .5, .8),
                        BFthres=5, PMPthres=.9, method="bfc",
                        attrition="weibull", params=c(.5,1),
                        t.points=c(0,1,2,3,4),
                        var.u0=0.01, var.u1=.1, var.e=.01, cov=0,
                        m=10000, log.grow=F, seed=NULL,
                        sensitivity=F, tol=.01,
                        N.max=1000, N.min=30, group_sizes = NULL) {

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
    if(m<1000) {message("Results with less than m=1000 generated datasets per iteration can be unreliable.")}
    if((method=="bf" | method=="BF") & (length(hypothesis)!=2)) {stop("Method 'bf' requires exactly two hypotheses.")}
    if (!is.null(group_sizes)) {
      if (length(group_sizes) != length(eff.sizes)) {stop("Length of 'group_sizes' must match number of conditions.")}
      if (any(group_sizes <= 0)) {stop("'group_sizes' must contain only positive values.")}
    }

    # extract variable names in order of appearance
    cond <- unique(unlist(strsplit(gsub("[^[:alnum:]_]", " ", hypothesis), "\\s+")))
    # check that the number of variables in hypothesis matches length of effect sizes
    if (length(cond) != length(eff.sizes)) {
      stop("Number of effect sizes must match number of unique variables in hypothesis.")
    }

    if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility

    N <- list()
    group_sizes <- group_sizes / sum(group_sizes)            # set group_sizes to correct scale if necessary
    n_cond <- length(eff.sizes)                              # extract number of conditions

    candidate_N <- seq(from = N.min,
                       to   = N.max,
                       by   = 1)

    condition <- FALSE                                        # condition initially FALSE until power criterion is reached
    stuck <- FALSE                                            # indicator if algorithm is stuck on repeated sample size
    j <- 1                                                    # iteration counter
    pow <- 0                                                  # initialize power
    av_it <- round(log(((N.max - N.min + 1)/n_cond), base=2)) # approximation of average numbers of iterations
    N_min <- N.min
    N_max <- N.max

    start_time <- Sys.time()

    if(sensitivity==F){
      ################### without sensitivity analysis #########################
      while(condition == F){

        N_mid <- round((N_min + N_max)/2 - .1, digits = 0)         # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
        N_tot <- candidate_N[which.min(abs(candidate_N - N_mid))]

        if (is.null(group_sizes)) { # balanced design
          N[[j]] <- N_tot
        } else { # unbalanced design via ratios
          N_raw <- N_tot * group_sizes
          N[[j]] <- floor(N_raw)
        }

        # set m according to iteration/difference between actual (pow) and desired power (eta)
        if(m>=5000){
          if(j==1 | abs(pow-eta) > .1){ # in the first iteration or if the difference between pow and eta is at least .1, set m to 1000
            current_m <- 1000
          } else {
            current_m <- m
          }
        } else {
          current_m <- m
        }

        # generate data and store BFs
        results <- get_power(attrition=attrition, params=params, m=current_m, N=N[[j]],
                             log.grow=log.grow, fraction=1,
                             t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                             cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                             BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

        # check if power is larger or smaller than desired power
        if(method == "bfc" | method == "BFc" | method == "bf_c" | method == "BF_c"){
          if(results$power_bfc >= eta){
            N_max <- sum(N[[j]]) - 1
          } else {
            N_min <- sum(N[[j]]) + 1
          }
          pow <- results$power_bfc
        } else if(method == "pmp" | method == "PMP"){
          if(results$power_pmp >= eta){
            N_max <- sum(N[[j]]) - 1
          } else {
            N_min <- sum(N[[j]]) + 1
          }
          pow <- results$power_pmp
        } else if(method == "bf" | method == "BF"){
          if(results$power_bf >= eta){
            N_max <- sum(N[[j]]) - 1
          } else{
            N_min <- sum(N[[j]]) + 1
          }
          pow <- results$power_bf
        }

        # Calculate time metrics
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        avg_time_per_iter <- elapsed / j
        remaining_time <- avg_time_per_iter * (av_it - j)

        # Print progress
        if(is.null(group_sizes)){
          cat(
            sprintf("Iteration %d: N = %d | Power = %.3f | Elapsed: %.1f minutes | Remaining: ~ %.1f minutes \n",
                    j, unlist(N[[j]]), pow, elapsed, remaining_time)
          )
        } else {
          cat(
            sprintf(
              "Iteration %d: N_total = %d | N_g = (%s) | Power = %.3f | Elapsed: %.1f minutes | Remaining: ~ %.1f minutes \n",
              j, sum(N[[j]]), paste(N[[j]], collapse = ", "), pow, elapsed, remaining_time)
          )
        }

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

        # if the current N is evaluated twice, condition is met
        if(length(N) > 2){
          if(sum(N[[j-1]]) == sum(N[[j-2]])){
            stuck <- TRUE
          }
        }
        # if power level is close enough to desired power level, condition is met and the algorithm stops
        if (round(abs(pow - eta), 8) <= tol | isTRUE(stuck)) {
          condition <- TRUE
          total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          if(is.null(group_sizes)){ # output for equal group sizes
            cat(
              sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f) \n",
                      j, total_time, unlist(N[[j]]), pow)
            )
          } else{ # output for unequal group sizes
            cat(
              sprintf("\nConverged in %d iterations (%.1f minutes). Final N_total = %d | N_g = (%s) | Power = %.3f \n",
                      j, total_time, sum(N[[j]]), paste(N[[j]], collapse = ", "), pow)
            )
          }
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
        stuck <- FALSE        # repeated sample size from previous iteration?
        j <- 1                # iteration counter
        pow <- 0              # initialize power
        av_it <- round(log((N.max - N.min + 1), base=2)) # approximation of average numbers of iterations
        N_min <- N.min
        N_max <- N.max
        # print info on sensitivity analysis
        if(i==1){
          cat("\n", "\n", "Sensitivity analysis for b = p/N_eff", "\n")
        } else {
          cat("\n", "\n", "Sensitivity analysis for b = ( p +", i, ")/N_eff", "\n")
        }

        while(condition == F){

          N_mid <- round((N_min + N_max)/2 - .1, digits = 0)         # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
          N_tot <- candidate_N[which.min(abs(candidate_N - N_mid))]

          if (is.null(group_sizes)) { # balanced design
            N[[j]] <- N_tot
          } else { # unbalanced design via ratios
            N_raw <- N_tot * group_sizes
            N[[j]] <- floor(N_raw)
          }

          # set m according to iteration/difference between actual (pow) and desired power (eta)
          if(m>=5000){
            if(j==1 | abs(pow-eta) > .1){ # in the first iteration or if the difference between pow and eta is at least .1, set m to 1000
              current_m <- 1000
            } else {
              current_m <- m
            }
          } else {
            currnet_m <- m
          }

          # generate data and store BFs
          results <- get_power(attrition=attrition, params=params, m=current_m, N=N[[j]],
                               log.grow=log.grow, fraction=i,
                               t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                               cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                               BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

          # check if power is larger or smaller than desired power
          if(method == "bfc" | method == "BFc" | method == "bf_c" | method == "BF_c"){
            if(results$power_bfc >= eta){
              N_max <- sum(N[[j]]) - 1
            } else {
              N_min <- sum(N[[j]]) + 1
            }
            pow <- results$power_bfc
          } else if(method == "pmp" | method == "PMP"){
            if(results$power_pmp >= eta){
              N_max <- sum(N[[j]]) - 1
            } else {
              N_min <- sum(N[[j]]) + 1
            }
            pow <- results$power_pmp
          } else if(method == "bf" | method == "BF"){
            if(results$power_bf >= eta){
              N_max <- sum(N[[j]]) - 1
            } else{
              N_min <- sum(N[[j]]) + 1
            }
            pow <- results$power_bf
          }

          # Calculate time metrics
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          avg_time_per_iter <- elapsed / j
          remaining_time <- avg_time_per_iter * (av_it - j)

          # Print progress
          if(is.null(group_sizes)){
             cat(
              sprintf("Iteration %d: N = %d | Power = %.3f | Elapsed: %.1f minutes | Total remaining: ~ %.1f minutes \n",
                      j, unlist(N[[j]]), pow, elapsed, remaining_time)
            )
          } else {
            cat(
              sprintf(
                "Iteration %d: N_total = %d | N_g = (%s) | Power = %.3f | Elapsed: %.1f minutes | Remaining: ~ %.1f minutes \n",
                j, sum(N[[j]]), paste(N[[j]], collapse = ", "), pow, elapsed, remaining_time)
            )
          }

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

          # if the current N is evaluated twice, condition is met
          if(length(N) > 2){
            if(sum(N[[j-1]]) == sum(N[[j-2]])){
              stuck <- TRUE
            }
          }

          # if power level is close enough to desired power level, condition is met and the algorithm stops
          if(round(abs(pow - eta), 8) <= tol | isTRUE(stuck)) {
            condition <- TRUE
            total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
            if(is.null(group_sizes)){ # output for equal group sizes
              cat(
                sprintf("\nConverged in %d iterations (%.1f minutes). Final N = %d (Power = %.3f) \n",
                        j, total_time, unlist(N[[j]]), pow)
              )
            } else{ # output for unequal group sizes
              cat(
                sprintf("\nConverged in %d iterations (%.1f minutes). Final N_total = %d | N_g = (%s) | Power = %.3f \n",
                        j, total_time, sum(N[[j]]), paste(N[[j]], collapse = ", "), pow)
              )
            }
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
