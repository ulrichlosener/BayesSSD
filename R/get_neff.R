#' Calculate the effective sample size in a multilevel model with attrition for multiple treatment groups
#'
#' @param model A multilevel model of class 'mermod' (specifically 'lmerMod') with ≥2 treatment conditions
#' @param N Number of subjects in the model
#' @param t.points Vector containing the measurement occasions
#' @param surviv Vector (or list of vectors) containing the proportion of subjects remaining at each timepoint
#'
#' @return Returns the effective sample size for each treatment condition's effects
#' @importFrom lme4 VarCorr getME
#' @importFrom Matrix bdiag Diagonal
#' @importFrom MASS ginv
#' @export
#'
#' @examples
#' # Example with 3 treatment conditions
#' library(lme4)
#' library(Matrix)
#' library(MASS)
#' get_neff_mis_mv(model = lmermodel, N=100, t.points=c(0,1,2,3,4),
#'                surviv=list(c(1,.9,.8,.6,.5), c(1,.85,.7,.5,.3), c(1,.95,.9,.8,.7)))

get_neff <- function(model, surviv) {

  # Extract number of subjects
  N <- as.numeric(lme4::ngrps(model))

  # Get the subject variable name
  subject_name <- names(lme4::ranef(model))

  # Count measurements per subject
  measurements_per_subject <- max(table(model@frame[[subject_name]]))

  # Extract all time points and the number of time points per person
  t.points <- unique(model@frame[["t"]])
  n <- length(t.points)

  # Extract number of treatment conditions
  treat_levels <- unique(model@frame[["treat"]])
  n_cond <- length(treat_levels)

  # Extract model components - optimized extraction
  reVar <- lme4::VarCorr(model)
  sigma2 <- sigma(model)^2
  Zt <- as.matrix(lme4::getME(model, "Zt")[1:2, 1:n, drop = FALSE])  # Correct matrix extraction
  Z <- t(Zt)

  # Create design matrices - optimized using matrix templates
  design_matrices <- vector("list", n_cond)
  base_design <- cbind(intercept = rep(1, n), t = t.points)

  # Create treatment indicator templates
  treat_template <- matrix(0, n, n_cond - 1)
  colnames(treat_template) <- paste0("treat", 1:(n_cond - 1))
  time_treat_template <- treat_template * t.points
  colnames(time_treat_template) <- paste0("t:treat", 1:(n_cond - 1))

  for(i in 1:n_cond) {
    if(i == 1) {
      # Control condition
      design_matrices[[i]] <- cbind(base_design, treat_template, time_treat_template)
    } else {
      # Treatment conditions
      treat_cols <- treat_template
      treat_cols[, i-1] <- 1
      time_treat_cols <- time_treat_template
      time_treat_cols[, i-1] <- t.points
      design_matrices[[i]] <- cbind(base_design, treat_cols, time_treat_cols)
    }
  }

  # Construct D matrix - optimized using direct extraction
  D <- as.matrix(Matrix::bdiag(reVar))

  # Precompute all possible subsets under attrition
  row_indices <- lapply(1:n, function(k) 1:k)

  # Pre-allocate and compute all inverses upfront
  V_inverses <- vector("list", n)
  W_inverses <- vector("list", n)

  for(k in 1:n) {
    rows <- row_indices[[k]]
    Z_sub <- Z[rows, , drop = FALSE]
    V <- Z_sub %*% D %*% t(Z_sub) + sigma2 * diag(k)

    # Handle diagonal matrix creation safely for large matrices
    W <- diag(x = diag(V), nrow = k, ncol = k)  # Safer diagonal matrix creation

    V_inverses[[k]] <- chol2inv(chol(V))  # Fast and  stable inverse
    W_inverses[[k]] <- diag(1/diag(V), nrow = k)  # Direct computation for diagonal
  }

  # Initialize result storage
  p <- ncol(design_matrices[[1]])
  indices <- c(2, seq(2 + n_cond, length.out = n_cond-1)) # positions of relevant estimates

  if(is.list(surviv) && length(surviv) > 1) {
    # Different survival patterns per condition
    res <- vector("list", n_cond)

    for(cond in 1:n_cond) {
      V_sum <- matrix(0, p, p)
      W_sum <- matrix(0, p, p)

      for(k in 1:n) {
        surv_weight <- surviv[[cond]][k]

        X_sub <- design_matrices[[cond]][row_indices[[k]], , drop = FALSE]
        Vinv <- V_inverses[[k]]
        Winv <- W_inverses[[k]]

        V_sum <- V_sum + surv_weight * crossprod(X_sub, Vinv %*% X_sub)
        W_sum <- W_sum + surv_weight * crossprod(X_sub, Winv %*% X_sub)
      }

      var_beta_hat <- tryCatch(solve(V_sum), error = function(e) MASS::ginv(V_sum))
      var_betahat_indep <- tryCatch(solve(W_sum), error = function(e) MASS::ginv(W_sum))

      w <- var_betahat_indep / var_beta_hat # apply formula by Faes et al. to calculate weight w

      # Extract number of non-missing observations for the current condition
      n_obs <- nrow(model@frame[model@frame$treat == treat_levels[cond],])

      N_eff <- w * n_obs # multiply weight w with number of observations in this sepcific condition

      res[[cond]] <- N_eff[indices[cond], indices[cond]] # return only N_eff for one condition
    }

    return(unlist(res))

  } else {
    # Same survival pattern
    V_sum <- matrix(0, p, p)
    W_sum <- matrix(0, p, p)

    # total number of observations
    for(k in 1:n) {
      surv_weight <- surviv[[1]][k]

      for(cond in 1:n_cond) {
        X_sub <- design_matrices[[cond]][row_indices[[k]], , drop = FALSE]
        V_sum <- V_sum + surv_weight * crossprod(X_sub, V_inverses[[k]] %*% X_sub)
        W_sum <- W_sum + surv_weight * crossprod(X_sub, W_inverses[[k]] %*% X_sub)
      }
    }

    var_beta_hat <- tryCatch(solve(V_sum), error = function(e) MASS::ginv(V_sum))
    var_betahat_indep <- tryCatch(solve(W_sum), error = function(e) MASS::ginv(W_sum))

    w <- var_betahat_indep / var_beta_hat # calculate weight w

    n_obs <- nrow(model@frame) # extract actual number of observations (attrition)

    N_eff <- w * n_obs / n_cond # divide N_eff by number of conditions, because they all have the same number of observations

    return(diag(N_eff)[indices])
  }
}


