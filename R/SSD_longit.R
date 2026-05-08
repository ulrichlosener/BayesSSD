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
#' @param group.sizes If group sizes are unequal, vector of proportions of group sizes for each condition. NULL if balanced design is to be used.
#'
#' @return Returns the sample size (number of subjects) necessary to achieve the desired power level `eta`.
#' @export
#' @examples
#' SSD_longit(eta=.8, attrition="weibull", params=c(.5,1),
#' m=100, t.points=c(0,1,2,3,4), var.u0=0.01,
#' var.u1=.1, var.e=.01, cov=0, eff.sizes=c(0, .5, .8),
#' BFthres=5, log.grow=F, seed=NULL,
#' hypothesis="a<b<c", PMPthres=.9, sensitivity=F, tol=.001,
#' N.max=1000)

SSD_longit <- function(eta=.8,
                       hypothesis="a<b<c",
                       eff.sizes=c(0, .5, .8),
                       t.points=c(0,1,2,3,4),
                       m=1000,
                       BFthres=5,
                       PMPthres=.9,
                       method="bfc",
                       attrition="weibull",
                       params=c(.5,1),
                       var.u0=0.033,
                       var.u1=.1,
                       var.e=.026,
                       cov=0,
                       log.grow=F,
                       seed=NULL,
                       sensitivity=F,
                       tol=.01,
                       N.max=1000,
                       N.min=30,
                       group.sizes=NULL) {

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
    if (!is.null(group.sizes)) {
      if (length(group.sizes) != length(eff.sizes)) {stop("Length of 'group.sizes' must match number of conditions.")}
      if (any(group.sizes <= 0)) {stop("'group.sizes' must contain only positive values.")}
    }

    # extract variable names in order of appearance
    cond <- unique(unlist(strsplit(gsub("[^[:alnum:]_]", " ", hypothesis), "\\s+")))
    # check that the number of variables in hypothesis matches length of effect sizes
    if (length(cond) != length(eff.sizes)) {
      stop("Number of effect sizes must match number of unique variables in hypothesis.")
    }

    if(!is.null(seed)) {set.seed(seed)}  # set user-specified seed for reproducibility

    N <- list()
    n_cond <- length(eff.sizes)                              # extract number of conditions

    candidate_N <- seq(from = N.min, to = N.max, by = 1)

    condition <- FALSE                                        # condition initially FALSE until power criterion is reached
    stuck <- FALSE                                            # indicator if algorithm is stuck on repeated sample size
    j <- 1                                                    # iteration counter
    av_it <- round(log(((N.max - N.min + 1)/n_cond), base=2)) # approximation of average numbers of iterations
    pb <- txtProgressBar(min = 0, max = 1, style = 3)         # progress bar
    pow <- list()                                             # initialize power
    prop_simple <- list()                                     # initialize proportion of simplified models
    N_min <- N.min                                            # initialize lower bound of sample sizes
    N_max <- N.max                                            # initialize upper bound of sample sizes

    start_time <- Sys.time()

    if(!sensitivity){
      ################### without sensitivity analysis #########################
      while(!condition){

        N_mid <- round((N_min + N_max)/2 - .1, digits = 0)         # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
        N_tot <- candidate_N[which.min(abs(candidate_N - N_mid))]

        if (is.null(group.sizes)) { # balanced design
          N[[j]] <- N_tot
        } else { # unbalanced design via ratios
          group.sizes <- group.sizes / sum(group.sizes)            # set group.sizes to correct scale if necessary
          N[[j]] <- floor(N_tot * group.sizes)
        }

        # set m according to distance to desired power
        current_m <- if (m >= 5000 && j != 1 && abs(pow[[j-1]] - eta) <= .1) {m} else {min(m, 1000)}

        # generate data and store BFs
        results <- get_power(attrition=attrition, params=params, m=current_m, N=N[[j]],
                             log.grow=log.grow, fraction=1,
                             t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                             cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                             BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)

        # extract correct power
        current_power <- switch(tolower(method),
                                "bfc" = results$power_bfc,
                                "bf_c" = results$power_bfc,
                                "pmp" = results$power_pmp,
                                "bf"  = results$power_bf)

        pow[[j]] <- current_power
        prop_simple[[j]] <- results$prop_simplified

        # update range of N_min and N_max
        if (current_power >= eta) {
          N_max <- sum(N[[j]]) - 1
        } else {
          N_min <- sum(N[[j]]) + 1
        }



        # evaluate power condition
        # if the current N is evaluated twice, or difference between N_min and N_max is 1 or less, condition is met
        if(length(N) > 2){
          if(sum(N[[j-1]]) == sum(N[[j-2]]) || abs(N_max-N_min) <= 1) {
            stuck <- TRUE
          }
        }

        # if power level is close enough to desired power level, condition is met and the algorithm stops
        if (round(abs(pow[[j]] - eta), 8) <= tol || stuck) {
          condition <- TRUE
          total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
        }

        progress <- min(j/av_it, 0.99)
        setTxtProgressBar(pb, progress) # update progress bar
        flush.console()          # clear console
        j <- j+1                 # update iteration number
      }

      # save results
      res <- list(
        evaluations = data.frame(
          N = unlist(N),
          power = round(unlist(pow), 2),
          prop_simplified = unlist(prop_simple)
        ),
        final = list(
          N = N[[j-1]],
          power = round(pow[[j-1]], 3),
          threshold_BF = BFthres,
          threshold_PMP = PMPthres,
          sensitivity = F
        ),
        hypotheses = list(
          hypothesis = hypothesis,
          comparison = method
        ),
        iterations = j-1,
        runtime_minutes = round(total_time)
      )

    } else {

      ##################### with sensitivity analysis ##########################
      run_ssd <- function(i){

        while(!condition){

          N_mid <- round((N_min + N_max)/2 - .1, digits = 0)         # current N (N_mid) is the mid point between N.min and N.max, rounded to the lower number
          N_tot <- candidate_N[which.min(abs(candidate_N - N_mid))]

          if (is.null(group.sizes)) { # balanced design
            N[[j]] <- N_tot
          } else { # unbalanced design via ratios
            N[[j]] <- floor(N_tot * group.sizes)
          }

          # set m according to distance to and desired power (eta)
          current_m <- if (m >= 5000 && j != 1 && abs(pow[[j-1]] - eta) <= .1) {m} else {min(m, 1000)}


          # generate data and store BFs
          results <- get_power(attrition=attrition, params=params, m=current_m, N=N[[j]],
                               log.grow=log.grow, fraction=i,
                               t.points=t.points, var.u0=var.u0, var.u1=var.u1,
                               cov=cov, var.e=var.e, eff.sizes=eff.sizes,
                               BFthres=BFthres, PMPthres=PMPthres, hypothesis=hypothesis)


          # extract correct power
          current_power <- switch(tolower(method),
                                  "bfc" = results$power_bfc,
                                  "bf_c" = results$power_bfc,
                                  "pmp" = results$power_pmp,
                                  "bf"  = results$power_bf)

          pow[[j]] <- current_power
          prop_simple[[j]] <- results$prop_simplified

          # update range of N_min and N_max
          if (current_power >= eta) {
            N_max <- sum(N[[j]]) - 1
          } else {
            N_min <- sum(N[[j]]) + 1
          }

          # if the current N is evaluated twice, condition is met
          if(length(N) > 2){
            if(sum(N[[j-1]]) == sum(N[[j-2]]) || abs(N_max-N_min) <= 1){
              stuck <- TRUE
            }
          }

          # if power level is close enough to desired power level, condition is met and the algorithm stops
          if(round(abs(pow[[j]] - eta), 8) <= tol || stuck) {
            condition <- TRUE
            total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
          }

          progress <- min(j/av_it, 0.99)
          setTxtProgressBar(pb, progress) # update progress bar
          flush.console()
          cat("  of analysis", i, "/ 3")
          j <- j+1                 # update iteration number
        }

        # return result for each b
        list(
          evaluations = data.frame(
            N = unlist(N),
            power = round(unlist(pow), 2),
            prop_simplified = unlist(prop_simple)
          ),
          final = list(
            N = N[[j-1]],
            power = round(pow[[j-1]], 3),
            threshold_BF = BFthres,
            threshold_PMP = PMPthres,
            sensitivity = T
          ),
          hypotheses = list(
            hypothesis = hypothesis,
            comparison = method
          ),
          iterations = j-1,
          runtime_minutes = round(total_time)
        )

      }
      res <- lapply(1:3, run_ssd)
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

  setTxtProgressBar(pb, 1) # set progressbar to maximum at the end of SSD
  close(pb) # close progressbar

  # print and return results
  print_results_SSD_longit(res)
  return(res)
}

# END OF FUNCTION --------------------------------------------------------------
