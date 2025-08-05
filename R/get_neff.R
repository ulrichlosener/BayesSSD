#' Calculate the effective sample size in a multilevel model with attrition for three treatment groups
#'
#' @param model A multilevel model of class 'mermod' (more specifically, an object of subclass 'lmerMod')
#' @param N Number of subjects in the model.
#' @param t.points Vector containing the measurement occasions.
#' @param surviv Vector (or list of vectors) containing the proportion of subjects remaining still in the study at each timepoint as returned by the function "survival".
#'
#' @return Returns the effective sample size for each of the three treatment conditions.
#' @importFrom lme4 VarCorr
#' @importFrom lme4 getME
#' @importFrom Matrix bdiag
#' @importFrom Matrix Diagonal
#' @importFrom MASS ginv
#' @export
#' @examples library(lme4)
#' library(Matrix)
#' library(MASS)
#' get_neff_mis_mv(model = lmermodel, N=100, t.points=c(0,1,2,3,4), surviv=c(1, .9, .8, .6, .5))

get_neff_mis_mv <- function(model, N, t.points, surviv) {

  n <- length(t.points)

  # Extraction of model components
  reVar <- lme4::VarCorr(model)
  sigma2 <- sigma(model)^2
  Zt <- as.matrix(lme4::getME(model, "Zt")[1:2, 1:n, drop = FALSE])
  Z <- t(Zt)

  # Design matrices for 3 conditions
  ones <- rep(1, n)
  zeros <- rep(0, n)

  # 1. Waiting list (reference: both dummies = 0)
  X_wl <- cbind(intercept = ones,
                t = t.points,
                treat1 = zeros,  # Dummy 1: treatment as usual
                treat2 = zeros,  # Dummy 2: intervention
                `t:treat1` = zeros,
                `t:treat2` = zeros)

  # 2. Treatment as usual (treat1 = 1, treat2 = 0)
  X_tau <- cbind(intercept = ones,
                 t = t.points,
                 treat1 = ones,
                 treat2 = zeros,
                 `t:treat1` = t.points,
                 `t:treat2` = zeros)

  # 3. Intervention (treat1 = 0, treat2 = 1)
  X_interv <- cbind(intercept = ones,
                    t = t.points,
                    treat1 = zeros,
                    treat2 = ones,
                    `t:treat1` = zeros,
                    `t:treat2` = t.points)

  # Construct D (variance-covariance matrix of random effects)
  D <- Matrix::bdiag(lapply(reVar, function(x) as(Matrix::bdiag(x), "generalMatrix")))

  # Precompute all possible subsets under attrition
  row_indices <- lapply(n:1, function(k) if(k == n) 1:n else 1:k)

  # Computation of V and W matrices
  V_subsets <- lapply(row_indices, function(rows) {
    Z_sub <- Z[rows, , drop = FALSE]
    as(Z_sub %*% D %*% t(Z_sub) + sigma2 * Matrix::Diagonal(length(rows)), "CsparseMatrix")
  })

  W_subsets <- lapply(V_subsets, function(V) {
    Matrix::Diagonal(x = diag(as.matrix(V)))
  })

  # Computation of inverses
  inverses <- lapply(seq_along(V_subsets), function(i) {
    list(
      V_inv = solve(V_subsets[[i]]),
      W_inv = solve(W_subsets[[i]])
    )
  })

  # Initialize sum matrices
  sum_mat <- sum_mat_indep <- matrix(0, 6, 6)

  # Compute sum matrices under dependence and independence assumption
  if(is.list(surviv) & length(surviv)>1){
    # Backwards loop for each treatment group with unique survival pattern
    res <- lapply(surviv, function(x) {
      for (k in n:1) {
        idx <- n - k + 1
        rows <- row_indices[[idx]]
        inv <- inverses[[idx]]

        # Get subset matrices
        X_wl_sub <- X_wl[rows, , drop = FALSE]
        X_tau_sub <- X_tau[rows, , drop = FALSE]
        X_interv_sub <- X_interv[rows, , drop = FALSE]

        # Take sums for all three conditions
        sum_mat <- sum_mat + x[k] * (
          t(X_wl_sub) %*% inv$V_inv %*% X_wl_sub +
            t(X_tau_sub) %*% inv$V_inv %*% X_tau_sub +
            t(X_interv_sub) %*% inv$V_inv %*% X_interv_sub
        )

        sum_mat_indep <- sum_mat_indep + x[k] * (
          t(X_wl_sub) %*% inv$W_inv %*% X_wl_sub +
            t(X_tau_sub) %*% inv$W_inv %*% X_tau_sub +
            t(X_interv_sub) %*% inv$W_inv %*% X_interv_sub
        )
      }

      # Matrix solving with fallback to pseudoinverse
      var_beta_hat <-  tryCatch(solve(sum_mat), error = function(e) MASS::ginv(sum_mat))
      var_betahat_indep <-  tryCatch(solve(sum_mat_indep), error = function(e) MASS::ginv(sum_mat_indep))

      # Compute N_eff
      w <- var_betahat_indep / var_beta_hat
      N_eff <- w * (N/3) * n # divide N by 3 because of 3 conditions

      return(list(N_eff[2,2], N_eff[5,5], N_eff[6,6]))
    })

    return(list(N_eff_cond1 = res[[1]][1],
                N_eff_cond2 = res[[2]][1],
                N_eff_cond3 = res[[3]][1]
    ))

  } else {
    # same pattern for all three groups
    for (k in n:1) {
      idx <- n - k + 1
      rows <- row_indices[[idx]]
      inv <- inverses[[idx]]

      # Get subset matrices
      X_wl_sub <- X_wl[rows, , drop = FALSE]
      X_tau_sub <- X_tau[rows, , drop = FALSE]
      X_interv_sub <- X_interv[rows, , drop = FALSE]

      # Take sums for all three conditions
      sum_mat <- sum_mat + surviv[[1]][k] * (
        t(X_wl_sub) %*% inv$V_inv %*% X_wl_sub +
          t(X_tau_sub) %*% inv$V_inv %*% X_tau_sub +
          t(X_interv_sub) %*% inv$V_inv %*% X_interv_sub
      )

      sum_mat_indep <- sum_mat_indep + surviv[[1]][k] * (
        t(X_wl_sub) %*% inv$W_inv %*% X_wl_sub +
          t(X_tau_sub) %*% inv$W_inv %*% X_tau_sub +
          t(X_interv_sub) %*% inv$W_inv %*% X_interv_sub
      )
    }

    # Matrix solving with fallback to pseudoinverse
    var_beta_hat <- tryCatch(solve(sum_mat), error = function(e) MASS::ginv(sum_mat))
    var_betahat_indep <- tryCatch(solve(sum_mat_indep), error = function(e) MASS::ginv(sum_mat_indep))

    # Compute N_eff
    w <- var_betahat_indep / var_beta_hat
    N_eff <- w * N * n

    return(list(N_eff_cond1 = N_eff[2,2]/3,
                N_eff_cond2 = N_eff[5,5]/3,
                N_eff_cond3 = N_eff[6,6]/3))
  }
}
