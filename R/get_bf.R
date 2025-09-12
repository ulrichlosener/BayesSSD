#' Calculate the Bayes Factor and posterior model probabilities
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
#' @importFrom bain bain
#' @export
#' @return Returns the Bayes Factor or Posterior Model Probabilities for the hypothesis.
#' @examples getbf_mis_mv(N=100, attrition="weibull", params=list(.8,1), hypothesis=list("a<b<c","a=b=c"), t.points=c(0,1,2), var.u0=.01, var.u1=.01, cov=0, var.e=.01, eff.sizes=c(0, .5, .8), fraction=1, log.grow=F)

get_bf <- function(N=100, attrition="weibull", params=c(.8,1), hypothesis=list("a<b<c","a=b=c"),
                                 t.points=c(0,1,2), var.u0=.01, var.u1=.002, cov=0, var.e=.01, eff.sizes=c(0.8, -.45, .-.46),
                                 fraction=1, log.grow=F){

  # determine number of conditions from hypotheses
  cond_letters <- unique(unlist(lapply(hypothesis, function(h) {
    unique(unlist(strsplit(gsub("[^a-z]", "", tolower(h)), "")))
  })))
  n_cond <- length(cond_letters)      # extract the number of conditions
  cond_letters <- sort(cond_letters)  # ensure consistent ordering

  # check consistency of number of conditions in arguments
  if(length(eff.sizes) != n_cond) {
    stop("Number of effect sizes must match number of conditions")
  }

  # check maximum number of conditions
  if(length(n_cond) > 10) {
    stop("The maximum number of treatment conditions is 10")
  }

  n <- length(t.points)  # number of measurement occasions

  # create time variable t
  if(log.grow==F) {
    t <- rep(t.points, N)
  } else {
    if(min(t.points)==0) {
      t <- rep(log(t.points+1), N)
    } else {
      t <- rep(log(t.points), N)
    }
  }

  t.prop <- t.points/max(t.points) # create rescaled time variable
  id <- rep(seq_len(N), each=n)  # create ID variable
  treat <- as.character(gl(n=n_cond, k=n, length=N*n, labels=cond_letters))
  dat0 <- data.frame(id, treat, t) # template data frame
  dat0$treat <- factor(dat0$treat, levels = cond_letters) # first letter as reference

  # generate data
  multinorm <- MASS::mvrnorm(n=N, mu=c(0,0), Sigma=matrix(c(var.u0, cov, cov, var.u1), 2, 2)) # draw random effects
  treat_coefs <- eff.sizes[-1] * sqrt(var.u1) + eff.sizes[1] * sqrt(var.u1) # cumulative effect sizes
  treat_dummies <- model.matrix(~ treat - 1, data = dat0)[, -1, drop = FALSE] # create dummies

  components <- list(
    base = eff.sizes[1] * sqrt(var.u1) * t, # baseline effect for condition "a"
    u0 = rep(multinorm[, 1], each=n),       # random intercepts
    u1 = rep(multinorm[, 2], each=n) * t,   # random slopes
    treatments = rowSums(treat_dummies * t %*% t(treat_coefs)), # betas for each condition
    error = rnorm(N * n, 0, sqrt(var.e))    # residual
  )

  y <- Reduce(`+`, components) # add everything together

  # compute survival curves
  if(!identical(attrition, FALSE)){ # with attrition
    surviv <- survival(distributions=attrition, params=params, t.points=t.points)
    shifted_surviv <- lapply(surviv, function(x) {c(x[-1], NA)})
  } else { # without attrition
    surviv <- list(rep(1, n))
    shifted_surviv <- c(surviv[-1], NA)
  }

  # compute hazard (hazard of t0 represents the probability of dropping out between t0 and t1)
  hazard <- mapply(function(a, b) {(a - b) / a},
                   surviv,
                   shifted_surviv,
                   SIMPLIFY = FALSE)

  if(length(hazard)>1){
    names(hazard) <- cond_letters
  }

  # introduce attrition to the simulated data
    dat <- data.frame(dat0, y, hazard = unsplit(hazard, f = dat0$treat))

    # extract subject ids
    ids <- unique(dat$id)

    # for each subject, find the first dropout time (if any)
    first_dropout <- tapply(
      1:nrow(dat),
      dat$id,
      function(idx) {
        hazards <- dat$hazard[idx]
        # Generate dropout decisions for all timepoints at once
        dropout <- rbinom(n-1, 1, hazards[-n])
        # Find the first dropout (returns NA if no dropout)
        which(dropout == 1)[1] + 1  # +1 because hazard[t] affects t+1
      },
      simplify = FALSE
    )

    # mark missing values and consecutive timepoints
    dat$mis <- 0
    dropout_indices <- unlist(
      lapply(ids, function(id) {
        fd <- first_dropout[[as.character(id)]]
        if (!is.na(fd)) {
          idx <- which(dat$id == id)
          seq(idx[fd], max(idx), by = 1)
        }
      })
    ) # remove marked y-values
    dat$mis[dropout_indices] <- 1
    dat$y[dat$mis == 1] <- NA

    simplified <- FALSE # initialize indicator variable for simplified models

  # in case of identifiability issues due to high attrition, simplify the model (constrain random effects to be independent)
  model <- tryCatch({
    lme4::lmer(formula = y ~ t + treat + t:treat + (1 + t | id),
               data = dat,
               REML = F,
               control = lme4::lmerControl(calc.derivs = F))
  }, error = function(e) {
    if (grepl("number of observations.*number of random effects", e$message)) {
      simplified <<- TRUE # indicate when simplidication happens
      tryCatch({
        lme4::lmer(formula = y ~ t + treat + t:treat + (1 | id) + (0 + t | id), # independent random effects
                   data = dat,
                   REML = F,
                   control = lme4::lmerControl(calc.derivs = F))
      })
    }
  })

  # check for rank deficiency
  expected_params <- 2 + (n_cond-1) * 2  # Intercept + t + (k-1) treat + (k-1) t:treat
  if(length(model@beta) < expected_params) { # in case the model is still not identifiable, throw error
    stop("Rank deficiency due to high attrition rate. Consider lowering attrition rate.")
  }

  # Extract estimates from the model
  est_indices <- c(2, seq(2 + n_cond, length.out = n_cond-1)) # positions of relevant estimates
  est <- model@beta[est_indices] # extract estimates
  names(est) <- cond_letters # name them for bain

  # Extract variances using the same indices as for the estimates
  Sigma <- lapply(est_indices, function(i) as.matrix(vcov(model)[i,i]))

  # calculate N_eff
  n_eff <- get_neff(model=model, t.points=t.points, surviv=surviv)

  # evaluate hypotheses
  hyp <- paste(hypothesis, collapse = ";") # put hypotheses in one single string for bain
  n_hyp <- length(hypothesis) # extract the number of hypotheses

  # perform Bayesian hypothesis evaluation using bain package
  bf_res <- bain::bain(x=est, Sigma=Sigma, n=n_eff,
                       hypothesis=hyp, group_parameters = 1, joint_parameters = 0)

  bf_c <- bf_res[["fit"]][["BF.c"]][1] # extract the BF versus the complement hypothesis
  PMPc <- bf_res[["fit"]][["PMPc"]][1] # extract the posterior model probabilities including all hypotheses plus the complement of the union

  bf12 <- NA
  if(n_hyp==2){ # if there are two hypotheses, extract the BF of H1 versus H2
    bf12 <- bf_res[["BFmatrix"]][1,2]
  }

  # return results
  return(list(bf_c=bf_c,
              PMPc=PMPc,
              bf12=bf12,
              simplified=simplified))
}
