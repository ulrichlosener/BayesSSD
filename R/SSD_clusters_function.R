########### Sample Size Determination for Cluster Randomized Trials ############

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description
#'@arguments
#  eff_size: Numeric. Effect size
#  n1: Numeric. Cluster size
#  n2: Numeric. Total number of clusters
#  ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size.
#  rho: Numeric. Intraclass correlation
#  BF_thresh: Numeric. Value of the Bayes factor that is going to be the threshold.
#  eta: Numeric. Probability of exceeding Bayes Factor threshold.
#  fixed: Character. Indicating which sample is fixed (n1 or n2)
#  b_fract: Numeric. Fraction of information used to specify the prior distribution.
#  max: Maximum sample size.
#  bacth_size: This parameter determines the size of batches used during the fitting of the multilevel model.


SSD_crt_null <- function(eff_size, n1 = 15, n2 = 30, ndatasets = 1000, rho, BF_thresh1,
                         BF_thresh0, eta1 = 0.8, eta0 = 0.8, fixed = "n2", b_fract = 3,
                         max = 1000, batch_size = 100) {
    # Libraries ----
    library(lme4)
    library(dplyr)

    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, rho, BF_thresh1, BF_thresh0, 
                     eta1, eta0, b_fract, max, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    if (rho > 1) stop("The intraclass correlation must be standardized and cannot be larger than 1")
    if (rho < 0) stop("The intraclass correlation must be a positive value")
    if (eta1 > 1 | eta0 > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta1 < 0 | eta0 < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if ((b_fract == round(b_fract)) == FALSE) stop("The fraction of information (b) must be an integer")

    #Functions ----------------
    source("data_generation.R")
    source("small_functions.R")
    source("print_results.R")
    source("aafbf.R")

    # Starting values ----------------------------------------------------------
    total_var <- 1
    var_u0 <- rho * total_var       #Between-cluster variance
    var_e <- total_var - var_u0     #Within-cluster variance
    eff_size0 <- 0                  #Effect size for null hypothesis
    conditions_met <- FALSE          #Indication we met the power criteria.
    ultimate_sample_sizes <- FALSE  #Indication that we found the required sample size.
    results_H0 <- matrix(NA, nrow = ndatasets, ncol = 4)
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)

    # Binary search start ------------------------------------------------------
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 5                     # Minimum cluster size
        low <- min_sample                   #lower bound
    }
    high <- max                    #higher bound

    # Hypotheses ----------------------------------------------------------------
    hypothesis1 <- "Intervention>Control"
    null <- "Intervention=Control"
    final_SSD <- vector(mode = "list", length = b_fract)
    type <- "Equality"
    b <- 1
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    singular_warn <- 0

    # Simulation of data and evaluation of condition  ----------------------------------
    while (ultimate_sample_sizes == FALSE) {
        # If H1 is true
        data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                              mean_interv = eff_size,
                                              batch_size = batch_size))
        
        # If H0 is true
        data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                              mean_interv = eff_size0,
                                              batch_size = batch_size))

        while (b < (b_fract + 1)) {
            #Approximated adjusted fractional Bayes factors------------------------------
            n_eff_H1 <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$rho_data)) / 2
            output_AAFBF_H1 <- Map(calc_aafbf, type, data_H1$estimates, data_H1$cov_list, list(b), n_eff_H1)
            
            n_eff_H0 <- ((n1 * n2) / (1 + (n1 - 1) * data_H0$rho_data)) / 2
            output_AAFBF_H0 <- Map(calc_aafbf, type, data_H0$estimates, data_H0$cov_list, list(b), n_eff_H0)
            
            # Results ---------------------------------------------------------------------
            results_H1[, 1] <- unlist(lapply(output_AAFBF_H1, extract_res, 1)) # Bayes factor H1vsH0
            results_H1[, 2] <- unlist(lapply(output_AAFBF_H1, extract_res, 4)) #posterior model probabilities of H1
            results_H1[, 3] <- unlist(lapply(output_AAFBF_H1, extract_res, 2)) # Bayes factor H0vsH1
            results_H1[, 4] <- unlist(lapply(output_AAFBF_H1, extract_res, 3)) #posterior model probabilities of H0
            
            colnames(results_H1) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
            
            results_H0[, 1] <- unlist(lapply(output_AAFBF_H0, extract_res, 1)) # Bayes factor H1vsH0
            results_H0[, 2] <- unlist(lapply(output_AAFBF_H0, extract_res, 4)) #posterior model probabilities of H1
            results_H0[, 3] <- unlist(lapply(output_AAFBF_H0, extract_res, 2)) # Bayes factor H0vsH1
            results_H0[, 4] <- unlist(lapply(output_AAFBF_H0, extract_res, 3)) #posterior model probabilities of H0
            
            colnames(results_H0) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
            
            #Evaluation of condition -------------------------------------------
            # Proportion
            prop_BF10 <- length(which(results_H1[, "BF.10"] > BF_thresh1)) / ndatasets
            prop_BF01 <- length(which(results_H0[, "BF.01"] > BF_thresh0)) / ndatasets
            
            # Evaluation
            ifelse(prop_BF01 > eta0 & prop_BF10 > eta1, conditions_met <- TRUE, conditions_met <- FALSE)
            previous_eta <- current_eta
            if (prop_BF01 < eta0 & prop_BF10 < eta1) {
                current_eta <- min(prop_BF10, prop_BF01)
            } else if (prop_BF01 < eta0 | prop_BF10 < eta1 ) {
                if (prop_BF01 < eta0) {
                    current_eta <- prop_BF01
                } else if (prop_BF10 < eta1) {
                    current_eta <- prop_BF10
                }
            } else if (conditions_met) {
                diff1 <- eta1 - prop_BF10
                diff0 <- eta0 - prop_BF01
                current_eta <- ifelse(diff1 < diff0, prop_BF10, prop_BF01)
                eta <- ifelse(current_eta == prop_BF10, eta1, eta0)
            }
            # print("Bayes factor check!")
            
            # Binary search algorithm ------------------------------------------
            if (conditions_met == FALSE) {
                # print(c("Using cluster size:", n1, "and number of clusters:", n2,
                #         "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10, "b:", b,
                #         "low:", low, "high:", high))
                # print("Increasing sample")
                if (fixed == "n1") {
                    if ((n2 == max) | (n2 > max))    { # If the sample size reaches the maximum
                        final_SSD[[b]] <- list("n1" = n1,
                                               "n2" = n2,
                                               "Proportion.BF01" = prop_BF01,
                                               "Proportion.BF10" = prop_BF10,
                                               "b.frac" = b,
                                               "data_H0" = results_H0,
                                               "data_H1" = results_H1,
                                               "singularity" = cbind(H0 = data_H0$singularity,
                                                                     H1 = data_H1$singularity))
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                    } else {
                        # Increase the number of clusters since eta is too small
                        low <- n2                         #lower bound
                        high <- high                      #higher bound
                        n2 <- round2((low + high) / 2)     #point in the middle
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1) # To ensure number of clusters is even
                        
                        # Adjust higher bound when there is a ceiling effect
                        if (low + n2 == high * 2) {
                            low <- n2                         #lower bound
                            if (previous_high > 0) {
                                high <- previous_high
                            } else {
                                high <- max                   #higher bound
                            }
                            n2 <- round2((low + high) / 2)     #point in the middle
                        }
                    }
                } else if (fixed == "n2") {
                    if ((n1 == max) | (n1 > max))    {# If the sample size reaches the maximum
                        final_SSD[[b]] <- list("n1" = n1,
                                               "n2" = n2,
                                               "Proportion.BF01" = prop_BF01,
                                               "Proportion.BF10" = prop_BF10,
                                               "b.frac" = b,
                                               "data_H0" = results_H0,
                                               "data_H1" = results_H1,
                                               "singularity" = cbind(H0 = data_H0$singularity,
                                                                     H1 = data_H1$singularity))
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                    } else {
                        # Increase the cluster sizes since eta is too small
                        low <- n1                        #lower bound
                        high <- high                     #higher bound
                        n1 <- round2((low + high) / 2)    #point in the middle
                        
                        # Adjust higher bound when there is a ceiling effect
                        if ((low + n1 == high * 2) | (current_eta == previous_eta)) {
                            low <- n1                        #lower bound
                            # Set the higher bound based on the previous high or the maximum
                            if (previous_high > 0 ) {
                                high <- previous_high
                            } else {
                                high <- max
                            }
                            n1 <- round2((low + high) / 2)    #point in the middle
                        }
                    }
                }
                break
            } else if (conditions_met == TRUE) {
                # print(c("Using cluster size:", n1,
                #         "and number of clusters:", n2,
                #         "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
                #         "low: ", low, "high: ", high, "b:", b))
                previous_high <- high
                SSD_object <- list("n1" = n1,
                                   "n2" = n2,
                                   "Proportion.BF01" = prop_BF01,
                                   "Proportion.BF10" = prop_BF10,
                                   "b.frac" = b,
                                   "data_H0" = results_H0,
                                   "data_H1" = results_H1,
                                   "singularity" = cbind(H0 = data_H0$singularity,
                                                         H1 = data_H1$singularity))
                # print("Lowerign sample")
                # print(c("previous:", previous_eta))
                previous_eta <- current_eta
                
                if (fixed == "n1") {
                    # Eta is close enough to the desired eta
                    if (current_eta - eta < 0.1 && n2 - low == 2) {
                        final_SSD[[b]] <- SSD_object
                        singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                        
                    } else if (previous_eta == current_eta && n2 - low == 2) {
                        # If there is no change in eta and the lower bound is close to the middle point
                        final_SSD[[b]] <- SSD_object
                        singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                        
                    } else {
                        # Decreasing to find the ultimate number of clusters
                        low <- low                         #lower bound
                        high <- n2                         #higher bound
                        n2 <- round2((low + high) / 2)      #point in the middle
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                        if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
                        # print("Lowering") # Eliminate later
                        break
                        
                    }
                } else if (fixed == "n2") {
                    # Eta is close enough to the desired eta
                    if (current_eta - eta < 0.1 && n1 - low == 1) {
                        final_SSD[[b]] <- SSD_object
                        singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                        
                    } else if (current_eta == previous_eta && n1 - low == 1) {
                        # If there is no change in eta and the lower bound is close to the middle point
                        final_SSD[[b]] <- SSD_object
                        singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                        b <- b + 1
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                        
                    } else if (current_eta == previous_eta && low + n2 == high * 2) {
                        # Reached the minimum number that meets the Bayesian power condition
                        final_SSD[[b]] <- SSD_object
                        b <- b + 1
                        singular_warn <- c(singular_warn, data_H0$singularity, data_H1$singularity)
                        low <- min_sample
                        previous_eta <- 0
                        previous_high <- 0
                        high <- max
                        next
                        
                    } else {
                        # Decreasing the cluster size to find the ultimate sample size
                        low <- low                         #lower bound
                        high <- n1                         #higher bound
                        n1 <- round2((low + high) / 2)      #point in the middle
                        # print("Lowering") # Eliminate later
                        break
                    }
                }
            } # Finish condition met
            
            # print(c("low:", low, "n2:", n2, "n1:", n1, "h:", high, "b:", b)) # Eliminate
        } # Finish while loop b
        
        # print(c("b fraction:", b))
        
        # Break loop
        if (b == b_fract + 1) {
            ultimate_sample_sizes <- TRUE
        }
        
        rm(data_H0, data_H1)
    } # Finish while loop ultimate_sample_size
    
    final_SSD[[b_fract + 1]] <- list(null, hypothesis1)
    final_SSD[[b_fract + 2]] <- list(BF_thresh0, BF_thresh1)
    final_SSD[[b_fract + 3]] <- list(eta0, eta1)
    
    # Final output -----
    print_results(final_SSD)
    if (any(singular_warn > 0)) warning("At least one of the fitted models is singular. For more information about singularity see help('isSingular').
                               The number of models that are singular can be found in the output object.")
    invisible(final_SSD)
}

# Test
# nulla <- SSD_crt_nullv2(eff_size = 0.5, ndatasets = 100, rho = 0.1, BF_thresh = 3, fixed = "n1",
#                       b_fract = 3)
# a <- SSD_crt_null(eff_size = 0.6, n1 = 15, n2 = 40, ndatasets = 15, rho = 0.05,
#                   BF_thresh1 = 3, BF_thresh0 = 1, eta1 = 0.8, eta0 = 0.7, fixed = "n1",
#                   b_fract = 2, max = 300, batch_size = 15)
