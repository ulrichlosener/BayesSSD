#' Print results for the function SSD_longit
#'
#'@title print results for SSD_longit
#'@description print and format results from the SSD_longit function.
#'@param result list containing the result of SSD_longit
#'@export

print_results_SSD_longit <- function(result) {

  ############################### without sensitivity analysis ###############################
  if(length(result) != 3) {
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")
    # print for BFc
    if (tolower(result$hypotheses$comparison) == "bfc"  |
        tolower(result$hypotheses$comparison) == "bf_c") {
      cat("Hypotheses:\n")
      cat("    H1:", result$hypotheses$hypothesis, "\n")
      cat("    Hc: not H1\n")
      cat("Using a total of N =", sum(result$final$N), "Participants\n")
      cat("Power = P(BF_1c >", result$final$threshold_BF, "| H1) =", result$final$power, "\n")
      # print for BF
    } else if (tolower(result$hypotheses$comparison) == "bf") {
      cat("Hypotheses:\n")
      cat("    H1:", result$hypotheses$hypothesis[[1]], "\n")
      cat("    H2:", result$hypotheses$hypothesis[[2]], "\n")
      cat("Using a total of N =", sum(result$final$N), "Participants\n")
      cat("Power = P(BF_12 >", result$final$threshold_BF, "| H1) =", result$final$power, "\n")
      # print for PMP
    } else if(tolower(result$hypotheses$comparison) == "pmp") {
      cat("Hypotheses:\n")
      cat(paste0(
        sprintf("    H%d: %s", seq_along(result$hypotheses$hypothesis), result$hypotheses$hypothesis),
        collapse = "\n"
      ), "\n")
      cat("Using a total of N =", sum(result$final$N), "Participants\n")
      cat("Power = P(PMP_H1 >", result$final$threshold_PMP, "| H1) =", result$final$power, "\n")
    }

  ############################### with sensitivity analysis ##################################
  } else {
    title <- "Final sample size with sensitivity analysis"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    row2 <- paste(rep("-", nchar(title)), collapse = "")
    cat(row, "\n")

      # print for BFc
    if (tolower(result[[1]]$hypotheses$comparison) == "bfc"  |
        tolower(result[[1]]$hypotheses$comparison) == "bf_c") {
      cat("Hypotheses:\n")
      cat("    H1:", result[[1]]$hypotheses$hypothesis, "\n")
      cat("    Hc: not H1\n")
      cat(row2, "\n")
        # for b1
      cat("For b = 1/N_eff:\n")
      cat("Using a total of N =", sum(result[[1]]$final$N), "Participants\n")
      cat("Power = P(BF_1c >", result[[1]]$final$threshold_BF, "| H1) =", result[[1]]$final$power, "\n")
      cat(row2, "\n")
        # for b2
      cat("For b = 2/N_eff:\n")
      cat("Using a total of N =", sum(result[[2]]$final$N), "Participants\n")
      cat("Power = P(BF_1c >", result[[2]]$final$threshold_BF, "| H1) =", result[[2]]$final$power, "\n")
      cat(row2, "\n")
        # for b3
      cat("For b = 3/N_eff:\n")
      cat("Using a total of N =", sum(result[[3]]$final$N), "Participants\n")
      cat("Power = P(BF_1c >", result[[3]]$final$threshold_BF, "| H1) =", result[[3]]$final$power, "\n")
      cat(row2, "\n")
        # print for BF
    } else if (tolower(result[[1]]$hypotheses$comparison) == "bf") {
      cat("Hypotheses:\n")
      cat("    H1:", result[[1]]$hypotheses$hypothesis[[1]], "\n")
      cat("    H2:", result[[1]]$hypotheses$hypothesis[[2]], "\n")
      cat(row2, "\n")
        # for b1
      cat("For b = 1/N_eff:\n")
      cat("    Using a total of N =", sum(result[[1]]$final$N), "Participants\n")
      cat("    Power = P(BF_12 >", result[[1]]$final$threshold_BF, "| H1) =", result[[1]]$final$power, "\n")
      cat(row2, "\n")
        # for b2
      cat("For b = 2/N_eff:\n")
      cat("    Using a total of N =", sum(result[[2]]$final$N), "Participants\n")
      cat("    Power = P(BF_12 >", result[[2]]$final$threshold_BF, "| H1) =", result[[2]]$final$power, "\n")
      cat(row2, "\n")
        # for b3
      cat("For b = 3/N_eff:\n")
      cat("    Using a total of N =", sum(result[[3]]$final$N), "Participants\n")
      cat("    Power = P(BF_12 >", result[[3]]$final$threshold_BF, "| H1) =", result[[3]]$final$power, "\n")
      cat(row2, "\n")
      # print for PMP
    } else if(tolower(result[[1]]$hypotheses$comparison) == "pmp") {
      cat("Hypotheses:\n")
      cat(paste0(
        sprintf("    H%d: %s", seq_along(result[[1]]$hypotheses$hypothesis), result[[1]]$hypotheses$hypothesis),
        collapse = "\n"
      ), "\n")
      cat(row2, "\n")
      # for b1
      cat("For b = 1/N_eff:\n")
      cat("    Using a total of N =", sum(result[[1]]$final$N), "Participants\n")
      cat("    Power = P(PMP_H1 >", result[[1]]$final$threshold_PMP, "| H1) =", result[[1]]$final$power, "\n")
      cat(row2, "\n")
      # for b2
      cat("For b = 2/N_eff:\n")
      cat("    Using a total of N =", sum(result[[2]]$final$N), "Participants\n")
      cat("    Power = P(PMP_H1 >", result[[2]]$final$threshold_PMP, "| H1) =", result[[2]]$final$power, "\n")
      cat(row2, "\n")
      # for b3
      cat("For b = 3/N_eff:\n")
      cat("    Using a total of N =", sum(result[[3]]$final$N), "Participants\n")
      cat("    Power = P(PMP_H1 >", result[[3]]$final$threshold_PMP, "| H1) =", result[[3]]$final$power, "\n")
      cat(row2, "\n")
    }
  }
}
