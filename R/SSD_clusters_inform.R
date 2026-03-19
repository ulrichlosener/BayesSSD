############# Sample Size Determination for Cluster Randomized Trials evaluating  #############
############################### informative hypotheses ######################################

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description 
#'@arguments
## eff.size: Numeric. Effect size
## n1: Numeric. Cluster size
## n2: Numeric. Total number of clusters.
## n.datasets: Numeric. Number of data sets that the user wants to generate to determine the sample size.
## rho: Numeric. Intraclass correlation
## BF.thresh: Numeric. Value of the Bayes factor that is going to be the threshold.
## eta: Numeric. Probability of exceeding Bayes Factor threshold.
## fixed: Character. Indicating which sample is fixed (n1 or n2)
## max: Maximum sample size.
## bacth_size: This parameter determines the size of batches used during the fitting of the multilevel model.

SSD_crt_inform <- function(eff_size, n1 = 15, n2 = 30, ndatasets = 1000, rho, BF_thresh,
                           eta = 0.8, fixed = 'n2', max = 1000, batch_size = 100) {
    # Libraries -----------------
    library(lme4)
    library(dplyr)
    
    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, rho, BF_thresh, eta, max, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    if (rho > 1) stop("The intraclass correlation must be standardized and cannot be larger than 1")
    if (rho < 0) stop("The intraclass correlation must be a positive value")
    if (eta > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    
    # Functions
    source('data_generation.R')
    source("small_functions.R")
    source("print_results.R")
    source("aafbf.R")
    
    # Starting values ----------------------------------------------------------
    total_var <- 1
    var_u0 <- rho * total_var       #Between-cluster variance
    var_e <- total_var - var_u0     #Within-cluster variance
    condition_met <- FALSE          #Indication we met the power criteria.
    ultimate_sample_sizes <- FALSE            #Indication that we found the required sample size.
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    
    # Binary search start ------------------------------------------------------
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum cluster size
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 5                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    }
    high <- max                    #higher bound
    
    #Hypotheses ----------------------------------------------------------------
    hypothesis1 <- "Intervention1 > Intervention2"
    hypothesis2 <- "Intervention1 < Intervention2"
    final_SSD <- vector(mode = "list", 3)
    type <- "Inequalities"
    previous_high <- 0
    previous_eta <- 0
    current_eta <- 0
    
    # Simulation of data and evaluation of condition  ----------------------------------
    while (ultimate_sample_sizes == FALSE) {
        # If H1 is true
        data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                              mean_interv = eff_size, 
                                              batch_size = batch_size))
        
        #Approximated adjusted fractional Bayes factors------------------------------
        n_eff_H1 <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$rho_data)) / 2
        output_AAFBF_H1 <- Map(calc_aafbf, type, data_H1$estimates, data_H1$cov_list, list(b), n_eff_H1)
        
        # Results ---------------------------------------------------------------------
        results_H1[, 1] <- unlist(lapply(output_AAFBF_H1, extract_res, 1)) # Bayes factor H1vsH2
        results_H1[, 2] <- unlist(lapply(output_AAFBF_H1, extract_res, 4)) #posterior model probabilities of H1
        results_H1[, 3] <- unlist(lapply(output_AAFBF_H1, extract_res, 2)) # Bayes factor H2vsH1
        results_H1[, 4] <- unlist(lapply(output_AAFBF_H1, extract_res, 3)) #posterior model probabilities of H2
        
        colnames(results_H1) <- c("BF.12", "PMP.1", "BF.21","PMP.2")
        
        #Evaluation of condition -------------------------------------------
        # Proportion
        prop_BF12 <- length(which(results_H1[, "BF.12"] > BF_thresh)) / ndatasets
        # Evaluation
        ifelse(prop_BF12 > eta, condition_met <- TRUE, condition_met <- FALSE)
        current_eta <- prop_BF12
        
        # Binary search algorithm ------------------------------------------
        if (condition_met == FALSE) {
            # print(c("Using cluster size:", n1, "and number of clusters:", n2,
            #         "prop_BF12: ", prop_BF12))
            if (fixed == "n1") {
                # Increase the number of clusters since eta is too small
                low <- n2                         #lower bound
                high <- high                      #higher bound
                n2 <- round((low + high) / 2)     #point in the middle
                ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1) # Ensure number of clusters is even
                
                # Adjust higher bound when there is a ceiling effect
                if (high + n2 == high * 2) {
                    low <- n2                         #lower bound
                    # Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0) {
                        high <- previous_high
                    } else {
                        high <- max
                    }
                    n2 <- round((low + high) / 2)     #point in the middle
                }
            } else if (fixed == "n2") {
                # Increase the cluster sizes since eta is too small
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round((low + high) / 2)    #point in the middle
                
                # Adjust higher bound when there is a ceiling effect
                if (high + n1 == high * 2) {
                    low <- n1                        #lower bound
                    # Set the higher bound based on the previous high or the maximum
                    if (previous_high > 0) {
                        high <- previous_high
                    } else {
                        high <- max
                    }
                    n1 <- round((low + high) / 2)    #point in the middle
                }
            }
        } else if (condition_met == TRUE) {
            # print(c("Using cluster size:", n1,
            #         "and number of clusters:", n2,
            #         "prop_BF12: ", prop_BF12,
            #         "low: ", low, "high: ", high))
            previous_high <- high
            if (fixed == "n1") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1) {
                    if (n2 - low == 2) {
                        ultimate_sample_sizes <- TRUE
                    } else {
                        # Decreasing with small steps to find the ultimate number of clusters
                        low <- low                         #lower bound
                        n2 <- n2 - 2
                        high <- (n2 * 2) - low
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                        if (n2 < 30) warning("The number of groups is less than 30.
                                                 This may cause problems in convergence and singularity.")
                    }
                    
                } else if (previous_eta == current_eta && n2 - low == 2) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    ultimate_sample_sizes = TRUE
                    
                } else {
                    # Decreasing to find the ultimate number of clusters
                    low <- low                         #lower bound
                    high <- n2                         #higher bound
                    n2 <- round((low + high) / 2)      #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                    if (n2 < 30) warning("The number of groups is less than 30.
                                             This may cause problems in convergence and singularity.")
                }
            } else if (fixed == "n2") {
                # Eta is close enough to the desired eta
                if (current_eta - eta < 0.1) {
                    if (n1 - low == 2) {
                        ultimate_sample_sizes <- TRUE
                        
                    } else {
                        # Decreasing with small steps to find the ultimate number of clusters
                        low <- low                     #lower bound
                        n1 <- n1 - 1 
                        high <- (n1*2) - low
                        if (n2 < 30) warning("The number of groups is less than 30.
                                                 This may cause problems in convergence and singularity.")
                        # print("Lowering with baby steps") # Eliminate late
                    }
                    
                } else if (current_eta == previous_eta && n1 - low == 1) {
                    # If there is no change in eta and the lower bound is close to the middle point
                    ultimate_sample_sizes <- TRUE
                    
                } else {
                    # Decreasing the cluster size to find the ultimate sample size
                    low <- low                         #lower bound
                    high <- n1                         #higher bound
                    n1 <- round((low + high) / 2)      #point in the middle
                    # print("Lowering") # Eliminate later
                }
            }
            
        } # Condition met
        # Break loop
        # If the sample size reaches the maximum
        previous_eta <- current_eta
        if (n2 == max) {
            break
        } else if (n1 == max) {
            break
        }
    } # Finish while loop ultimate_sample_size
    SSD_object <- list("n1" = n1,
                       "n2" = n2,
                       "Proportion.BF12" = prop_BF12,
                       "data_H1" = results_H1,
                       "HYpotheses" = list(hypothesis1, hypothesis2),
                       "BF.threshold" = BF_thresh,
                       "Constraint" = type)
    
    # Final output -----
    print_results(SSD_object)
    invisible(SSD_object)
}


# Test -------------------------------------------------------------------------
# start.time <- Sys.time()
# a <- SSD_crt_inform(eff.size = 0.5, n.datasets = 100, rho = 0.1, BF.thresh = 3, fixed = 'n1')
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# a <- SSD_crt_inform(eff_size = 0.6, n1 = 15, n2 = 40, ndatasets = 10, rho = 0.05, 
#                     BF_thresh = )
# 

# 
