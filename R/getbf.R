#' Calculate the BF for a hypothesis with attrition
#'
#' @param N The number of subjects
#' @param attrition The attrition pattern (FALSE for no attrition, otherwise "weibull", "modified_weibull", "linear_exponential", "log_logistic", "gompertz" or "non-parametric")
#' @param params The parameters passed to the survival function specified in "attrition". First parameter is omega, second is gamma.
#' @param hypothesis The hypothesis to be evaluated. Treatment groups are coded as "a", "b", "c", etc.
#' @param t.points The points in time of measurement. Can be non-equidistant.
#' @param var.u0 The intercept variance.
#' @param var.u1 The slope variance.
#' @param cov The covariance between intercept variance and slope variance.
#' @param var.e The residual variance.
#' @param eff.sizes The effect sizes defined as the differences between the regression coefficients of interaction between time and condition.
#' @param fraction The fraction of information used to construct the prior for the Bayes Factor.
#' @param log.grow Use log-linear growth?
#'
#' @return Returns the Bayes Factor or Posterior Model Probabilities for the hypothesis.
#'
#' @examples getbf_mis_mv()

getbf_mis_mv <- function(N, attrition="weibull", params=list(.5,1), hypothesis="a<b<c",
                         t.points, var.u0, var.u1, cov, var.e, eff.sizes=c(0, .5, .8), fraction=1,
                         log.grow){

  n <- length(t.points)  # number of measurement occasions
  # create time variable t
  if(log.grow==F) {t <- rep(t.points, N) # if logarithmic growth is used, take log of t.points
  } else {if(min(t.points)==0) {t <- rep(log(t.points+1), N)} else {t <- rep(log(t.points), N)} # if the first timepoint is zero, we add 1 to all timepoints because log(0) is undefined
  }
  t.prop <- t.points/max(t.points) # create rescaled time variable (from 0 to 1)
  id <- rep(seq_len(N), each=n)  # create ID variable
  treat <- as.character(gl(n=3, k=n, length=N*n, labels=c("a","b","c"))) # create treatment variable
  dat0 <- data.frame(id, treat, t) # combine into data frame
  dat0$treat <- factor(dat0$treat, levels = c("a", "b", "c")) # Forces "a" as reference
  multinorm <- MASS::mvrnorm(n=2*N, mu=c(0,0), matrix(c(var.u0, cov, cov, var.u1), nrow=2, ncol=2)) # draw random effects

  # make params into list of one vector
  if (is.list(params)) {
    params <- list(unlist(params, use.names = FALSE))
  } else {
    params <- list(params)
  }

  # Compute survival curves
  if(any(attrition != F)){
    surviv <- survival(attrition, params, t.points)
    shifted_surviv <- lapply(surviv, function(x) {c(x[-1], NA)}) # drop first element and append NA to compute hazard
  } else if(attrition==F){
    surviv <- list(rep(1, n))
    shifted_surviv <- c(surviv[-1], NA)
  }

  # compute hazard: (S_t - S_{t+1}) / S_t
  hazard <- mapply(function(a, b) {(a - b) / a},
                   surviv,
                   shifted_surviv,
                   SIMPLIFY = FALSE)

  # generate data under the research hypothesis
  u0 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),1], each=n) # random intercepts
  u1 <- rep(multinorm[(nrow(multinorm)/2+1):(nrow(multinorm)),2], each=n) # random slopes
  e <- rnorm(N*n, 0, sqrt(var.e)) # error variance for H1

  # Create treatment dummy variables
  treat_TAU <- as.numeric(treat == "b")
  treat_INT <- as.numeric(treat == "c")

  beta1_WL <- unlist(eff.sizes[1]) * sqrt(var.u1)
  beta2_TAU <- unlist(eff.sizes[2]) * sqrt(var.u1) + beta1_WL # create coefficient of interaction for b from effect size
  beta3_INT <- unlist(eff.sizes[3]) * sqrt(var.u1) + beta1_WL # create coefficient of interaction for c from effect size

  y <- beta1_WL*t + u0 + beta2_TAU*treat_TAU*t + beta3_INT*treat_INT*t + u1*t + e # create outcome variable y

  if(attrition != F){
    dat <- data.frame(dat0, y, hazard=unlist(hazard))
    suppressWarnings(dat$mis <- rbinom(n=nrow(dat), size=1, prob=dat$hazard)) # suppress warning about NAs being produced
    dat <- data.frame(dat %>% group_by(id) %>% mutate(mis = ifelse(cumany(mis == 1), 1, mis)))
    dat$y[which(dat$mis==1)] <- NA
  } else {
    dat <- data.frame(dat0, y)
  }

  # fit MLM to dataset
  model <- lme4::lmer(formula = y ~ t + treat + t:treat + (1 + t | id), data = dat, REML=F, control = lme4::lmerControl(calc.derivs = F))  # fit MLM model under H1
  est <- model@beta[c(2,5,6)] # extract estimates of beta2 and beta3 under H0
  names(est) <- c("a", "b", "c")
  sig_WL <- as.matrix(vcov(model)[2,2])  # extract variance of estimates under H0
  sig_TAU <- as.matrix(vcov(model)[5,5])  # extract variance of estimates under H0
  sig_INT <- as.matrix(vcov(model)[6,6])  # extract variance of estimates under H0

  # calculate N_eff
  n_eff <- get_neff_mis_mv(model=model, N=N, t.points=t.points, surviv=surviv)

  bf_res <- bain(x=est, Sigma=list(sig_WL, sig_TAU, sig_INT), n=unlist(n_eff),
                 hypothesis=hypothesis, group_parameters = 1, joint_parameters = 0)

  bf_c <- bf_res[["fit"]][["BF.c"]][1]
  bfs <- bf_res[["BFmatrix"]]

  PMPc <- bf_res[["fit"]][["PMPc"]][1]

  return(output = list(bf_c=bf_c, PMPc=PMPc))
}
