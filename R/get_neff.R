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

get_neff <- function(model, N, t.points, surviv) {

  # Check number of treatment conditions
  treat_levels <- unique(model@frame[["treat"[1]]])
  num_conditions <- length(treat_levels)

  n <- length(t.points)
  ones <- rep(1, n)
  zeros <- rep(0, n)

  # Extract model components
  reVar <- lme4::VarCorr(model)
  sigma2 <- sigma(model)^2
  Zt <- as.matrix(lme4::getME(model, "Zt")[1:2, 1:n, drop = FALSE])
  Z <- t(Zt)

  # Create design matrices for each condition
  design_matrices <- lapply(1:num_conditions, function(i) {
    treat_cols <- matrix(0, nrow = n, ncol = num_conditions - 1)
      if(i > 1) treat_cols[, i-1] <- 1
      colnames(treat_cols) <- paste0("treat", 1:(num_conditions-1))

    time_treat_cols <- matrix(0, nrow = n, ncol = num_conditions - 1)
      if(i > 1) time_treat_cols[, i-1] <- t.points
      colnames(time_treat_cols) <- paste0("t:treat", 1:(num_conditions-1))

    cbind(intercept = ones,
          t = t.points,
          treat_cols,
          time_treat_cols)
  })

  # Construct D matrix
  D <- Matrix::bdiag(lapply(reVar, function(x) as(Matrix::bdiag(x), "generalMatrix")))

  # Precompute all possible subsets under attrition
  row_indices <- lapply(n:1, function(k) if(k == n) 1:n else 1:k)

  # Compute V and W matrices for all subsets
  compute_matrices <- function(rows) {
    Z_sub <- Z[rows, , drop = FALSE]
    V <- as(Z_sub %*% D %*% t(Z_sub) + sigma2 * Matrix::Diagonal(length(rows)), "CsparseMatrix")
    W <- Matrix::Diagonal(x = diag(as.matrix(V)))
    list(V = V, W = W)
  }

  matrices <- lapply(row_indices, compute_matrices)

  # Compute inverses for all subsets
  inverses <- lapply(matrices, function(m) {
    list(V_inv = solve(m$V),
         W_inv = solve(m$W))
  })

  # Initialize sum matrices
  p <- ncol(design_matrices[[1]])
  sum_mat <- sum_mat_indep <- matrix(0, p, p)

  if(is.list(surviv) && length(surviv) > 1) {
    # Different survival patterns for each condition
    process_condition <- function(cond) {
      cond_sums <- Reduce(
        function(acc, k) {
          idx <- n - k + 1
          rows <- row_indices[[idx]]
          inv <- inverses[[idx]]
          X_sub <- design_matrices[[cond]][rows, , drop = FALSE]
          surv_weight <- surviv[[cond]][k]

          list(
            V_sum = acc$V_sum + surv_weight * (t(X_sub) %*% inv$V_inv %*% X_sub),
            W_sum = acc$W_sum + surv_weight * (t(X_sub) %*% inv$W_inv %*% X_sub)
          )
        },
        n:1,
        init = list(V_sum = matrix(0, p, p), W_sum = matrix(0, p, p))
      )

      var_beta_hat <- tryCatch(solve(cond_sums$V_sum), error = function(e) MASS::ginv(cond_sums$V_sum))
      var_betahat_indep <- tryCatch(solve(cond_sums$W_sum), error = function(e) MASS::ginv(cond_sums$W_sum))

      w <- var_betahat_indep / var_beta_hat
      N_eff <- w * (N/num_conditions) * n

      # Extract relevant effects
      effects <- c(
        Time = N_eff[2,2],
          c(
            sapply(1:(num_conditions-1), function(i) c(
              Main = N_eff[2+i,2+i],
              Interaction = N_eff[2+num_conditions-1+i, 2+num_conditions-1+i]
            ))
          )
      )
      effects
    }

    res <- lapply(1:num_conditions, process_condition)
    names(res) <- paste0("Condition_", treat_levels)
    return(res)

  } else {
    # Same survival pattern for all conditions
    sums <- Reduce(
      function(acc, k) {
        idx <- n - k + 1
        rows <- row_indices[[idx]]
        inv <- inverses[[idx]]
        surv_weight <- surviv[[1]][k]

        cond_sums <- lapply(design_matrices, function(X) {
          X_sub <- X[rows, , drop = FALSE]
          list(
            V = surv_weight * (t(X_sub) %*% inv$V_inv %*% X_sub),
            W = surv_weight * (t(X_sub) %*% inv$W_inv %*% X_sub)
          )
        })

        list(
          V_sum = acc$V_sum + Reduce(`+`, lapply(cond_sums, `[[`, "V")),
          W_sum = acc$W_sum + Reduce(`+`, lapply(cond_sums, `[[`, "W"))
        )
      },
      n:1,
      init = list(V_sum = matrix(0, p, p), W_sum = matrix(0, p, p))
    )

    var_beta_hat <- tryCatch(solve(sums$V_sum), error = function(e) MASS::ginv(sums$V_sum))
    var_betahat_indep <- tryCatch(solve(sums$W_sum), error = function(e) MASS::ginv(sums$W_sum))

    w <- var_betahat_indep / var_beta_hat
    N_eff <- w * N * n

    # Format results
    indices <- c(2, seq(2 + num_conditions, length.out = num_conditions-1)) # positions of relevant estimates

    return(N_eff[indices, indices]/num_conditions)
  }
}
