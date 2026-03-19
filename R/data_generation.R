######################## DATA GENERATION ##############################

#"@title Generation of data sets for two treatment-condition cluster randomized trial 
#"@description
#"@arguments
## ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size
## n1: Numeric. Cluster size
## n2: Numeric. Total number of clusters.
## rho: Numeric. Intraclass correlation
## var_u0: Numeric. Between-cluster variance. Variance at the cluster level.
## var_e: Numeric. Within-cluster variance.
## mean_interv: Numeric. Equivalent to the effect size.
## bacth_size: This parameter determines the size of batches used during the fitting of the multilevel model.



gen_CRT_data <- function(ndatasets = ndatasets, n1 = n1, n2 = n2, var_u0 = var_u0,
                         var_e = var_e, mean_interv, batch_size) {
  # Create variables id  of the cluster and condition
  id <- rep(1:n2, each = n1)
  if (n2 %% 2 == 0) {
    condition <- rep(c(0, 1), each = n1 * n2 / 2)
  } else {
    # Odd number of clusters
    # the extra cluster goes to control condition
    half <- floor(n2 / 2)
    condition <- c(rep(0, n1 * half), rep(1, n1 * half), rep(0, n1))
  }
  # Dummy variables for no intercept model
  intervention <- condition
  control <- 1 - intervention
  mean_control <- 0
  marker <- 0
  
  # Tables for results
  output_lmer <- vector(mode = "list", length = ndatasets)
  data_list <- vector(mode = "list", length = ndatasets)
  
  # Data generation ----------------------------------------------------------
  for (iter in seq(ndatasets)) {
    set.seed((iter + 90) * iter)
    u0 <- rnorm(n2, 0, sqrt(var_u0))
    u0 <- rep(u0, each = n1)
    e <- rnorm(n1 * n2, 0, sqrt(var_e))
    resp <- mean_control * control + mean_interv * intervention + u0 + e
    #Data frame
    data_list[[iter]] <- cbind(resp, intervention, control, id)
  }
  data_list <- lapply(data_list, as.data.frame) # Necessary to use lmer
  
  # Multilevel analysis --------------------------------------------------------
  # Batches
  batch_size <- batch_size
  ifelse((ndatasets / batch_size) %% 1 == 0, num_batches <- ndatasets / batch_size,
         num_batches <- (ndatasets / batch_size) + 1)
  for (batch in seq(num_batches)) {
    #Indexes
    start_index <- (batch_size * (batch - 1)) + 1
    end_index <- min(batch * batch_size, ndatasets)
    #Multilevel fitting
    output_lmer[start_index:end_index] <- lapply(data_list[start_index:end_index], fit_lmer)
  }
  marker <- lapply(output_lmer, marker_func) # Mark singularity
  singular_datasets <- Reduce("+", marker) # How many are singular?
  
  estimates <- lapply(output_lmer, fixef)                 # Means
  cov_intervention <- lapply(output_lmer, varcov, 1)      # Covariance
  cov_control <- lapply(output_lmer, varcov, 4)
  cov_list <- Map(list, cov_intervention, cov_control)
  var_u0_data <- unlist(lapply(output_lmer, get_variances, 1))
  var_e_data <- unlist(lapply(output_lmer, get_variances, 2))
  total_var_data <- var_u0_data + var_e_data
  rho_data <- var_u0_data / total_var_data
  # print("Multilevel check")

  rm(id, condition, intervention, control, mean_control, data_list, output_lmer,
     u0, e, resp, batch_size, cov_intervention, cov_control, var_u0_data,
     var_e_data, total_var_data)
  
  return(output <- list("rho_data" = rho_data,
                        "estimates" = estimates,
                        "cov_list" = cov_list,
                        "singularity" = singular_datasets))
}
