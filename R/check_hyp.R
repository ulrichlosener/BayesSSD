#' Check if hypothesis is true according to 'eff.sizes'
#'
#' @param hypothesis The research hypothesis of interest, which must be true
#' @param eff.sizes The vector of true effect sizes
#'
#' @returns TRUE if the hypothesis is true or FALSE if it is not true
#' @export
#'
#' @examples check_hyp("a<b<c", c(0, .5, .8)) # returns TRUE

check_hyp <- function(hypothesis, eff.sizes) {

  # eliminate whitespace
  hypothesis <- gsub("\\s+", "", hypothesis)

  # split constraints by ';' or '&'
  constraints <- unlist(strsplit(hypothesis, "\\s*[;&]\\s*"))

  # extract variable names in order of appearance
  vars <- unique(unlist(strsplit(gsub("[^[:alnum:]_]", " ", hypothesis), "\\s+")))
  # check that the number of variables in hypothesis matches length of effect sizes
  if (length(vars) != length(eff.sizes)) {
    stop("Number of effect sizes must match number of unique variables in hypothesis.")
  }
  names(eff.sizes) <- vars

  # helper to check one constraint
  check_one <- function(constr) {
    # split into tokens (variable names and operators)
    tokens <- unlist(strsplit(constr, "(?<=\\w)(?=[<>=])|(?<=[<>=])(?![=])|(?<=<>)(?=\\w)", perl = TRUE))
    tokens <- tokens[tokens != ""]

    # alternate between variable and operator
    for (i in seq(1, length(tokens) - 2, by = 2)) {
      left <- tokens[i]
      op <- tokens[i + 1]
      right <- tokens[i + 2]
      a <- eff.sizes[left]
      b <- eff.sizes[right]

      ok <- switch(op,
                   "<"  = a < b,
                   ">"  = a > b,
                   "="  = abs(a - b) < .Machine$double.eps^0.5,
                   stop(paste("Unknown operator:", op))
      )
      if (!ok) return(FALSE) # if hypothesis is not true
    }
    TRUE # is hypothesis is true
  }
  all(sapply(constraints, check_one)) # apply to all constraints
}

# END OF FUNCTION --------------------------------------------------------------
