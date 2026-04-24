############################## PRINT RESULTS #################################
#'@title print results for SSD_longit
#'@description print results
#'@param result object containing the result of SSD_longit
#'@export

print_results_SSD_longit <- function(result) {
  title <- "Final sample size"
  cat(paste("\n", title, "\n", sep = ""))
  row <- paste(rep("=", nchar(title)), collapse = "")
  cat(row, "\n")
  if (result$hypotheses$comparison == "bfc"  |
      result$hypotheses$comparison == "BFc"  |
      result$hypotheses$comparison == "bf_c" |
      result$hypotheses$comparison == "BF_c") {   # Print for BFc
    cat("Hypotheses:", "\n")
    cat("    H1:", result$hypotheses$hypothesis, "\n")
    cat("    Hc: not H1\n")
    cat("Using  a total of N =", sum(result$final$N), "Participants\n")
    cat("P (BF_1c >", result$final$threshold, "| H1) = ", res$final$power, "\n")
  } else {                                        # NOT YET IMPLEMENTED
    n_object <- length(object_result)
    b_number <- n_object - 3
    results_matrix <- matrix(NA, nrow = b_number, ncol = 5)
    results_matrix[, 1] <- seq(b_number)
    object_result_b <- object_result[1:b_number]
    results_matrix[, 2] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 1)), b_number))) #n1
    results_matrix[, 3] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 2)), b_number))) #n2
    results_matrix[, 4] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 3)), b_number)), 3) #BF_01
    results_matrix[, 5] <- round(as.numeric(head(unlist(lapply(object_result_b, `[[`, 4)), b_number)), 3) #BF_10
    colnames(results_matrix) <- c("b", "n1", "n2", paste("P(BF.01 >", object_result[[(n_object - 1)]][[1]], "| H0) > ", object_result[[n_object]][[1]], sep = " "),
                                  paste("P(BF.10 >", object_result[[(n_object - 1)]][[2]], "| H1) > ", object_result[[n_object]][[2]], sep = " "))

    cat("Hypotheses:", "\n")
    cat("    H0:", object_result[[b_number + 1]][[1]], "\n")
    cat("    H1:", object_result[[b_number + 1]][[2]], "\n")

    cat("***********************************************************************", "\n")
    print(format(results_matrix, justify = "centre"))
    cat("***********************************************************************", "\n")
    cat("n1: Cluster sizes", "\n")
    cat("n2: Number of clusters", "\n")
  }
}
