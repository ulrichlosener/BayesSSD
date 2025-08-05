
#' Calculate survival curves to model attrition
#'
#' @param distributions Can be either a vector of length 1 or of length g, where g is the number of treatment conditions. Specifies the distribution(s) used to model attrition.
#' @param params List with either 1 element or g elements. Each element contains the parameters for the respective survival function.
#' @param t.points Vector containing the measurement occasions.
#' @return Returns the proportion of individuals still remaining in the study at each timepoint.
#'
#' @examples survival(distributions="weibull", params=list(c(.5, 1)), t.points=c(0,1,2,3,4))

survival <- function(distributions, params, t.points) {

  # Convert single distribution/params to list format
  if (is.character(distributions) && length(distributions) == 1) {
    distributions <- list(distributions)
  }
  if (is.numeric(params)) {
    params <- list(params)
  }

  # Validate inputs
  if (!is.list(params) || !all(sapply(params, is.numeric))) {
    stop("params must be a list of numeric vectors")
  }

  # Convert time points to proportional time (0-1)
  time <- t.points/max(t.points)

  # Define distribution functions
  dist_functions <- list(
    weibull = function(pars, time) {
      if (length(pars) < 2) stop("Weibull needs omega, gamma")
      (1-pars[1])^(time^pars[2])
    },
    modified_weibull = function(pars, time) {
      if (length(pars) < 3) stop("Modified Weibull needs omega, gamma, kappa")
      exp(time^pars[2]*exp(pars[3]*(time-1))*log(1-pars[1]))
    },
    linear_exponential = function(pars, time) {
      if (length(pars) < 2) stop("Linear-exponential needs omega, gamma")
      exp((.5*pars[2]+log(1-pars[1]))*time - .5*pars[2]*time^2)
    },
    log_logistic = function(pars, time) {
      if (length(pars) < 2) stop("Log-logistic needs omega, gamma")
      (1-pars[1])/((1-pars[1]) + pars[1]*time^pars[2])
    },
    gompertz = function(pars, time) {
      if (length(pars) < 2) stop("Gompertz needs omega, gamma")
      exp((log(1-pars[1])/(exp(pars[2])-1))*(exp(pars[2]*time)-1))
    },
    nonparametric = function(pars, time) {
      if (time[1] == 0) {
        if (length(pars) < max(time)+1) stop("Need ", max(time)+1, " params for t=0 start")
        pars
      } else {
        if (length(pars) < max(time)) stop("Need ", max(time), " params for t=1 start")
        pars
      }
    }
  )

  # Calculate survival curve(s)
  results <- Map(function(dist, pars) {
    if (!dist %in% names(dist_functions)) {
      stop("Unknown distribution: ", dist)
    }
    dist_functions[[dist]](pars, time)
  }, distributions, params)

  return(results)
}
