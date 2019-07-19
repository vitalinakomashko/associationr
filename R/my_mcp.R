my_mcp <- function (linfct,
                    interaction_average = FALSE,
                    covariate_average = FALSE) {
  # comment out library loading
  # library(multcomp)
  linfct <- lapply(linfct, function(x) {
    if (is.numeric(x) && !is.matrix(x)) {
      return(matrix(x, nrow = 1))
    } else {
      return(x)
    }
  })
  if (is.null(names(linfct)))
    stop(sQuote("linfct"), " doesn't have a ", sQuote("names"),
         " attribute")
  classes <- sapply(linfct, function(x) inherits(x, "matrix") ||
                      inherits(x, "character"))
  if (length(linfct) == 1 && linfct[[1]] == "Means") {
    class(linfct) <- "means"
    return(linfct)
  }
  attr(linfct, "interaction_average") <- interaction_average
  attr(linfct, "covariate_average") <- covariate_average
  if (all(classes)) {
    class(linfct) <- "mcp"
    return(linfct)
  }
  stop("Arguments don't consist of either matrices or characters")
}
