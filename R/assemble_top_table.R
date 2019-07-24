assemble_top_table <- function(lm.obj,
                       contrast.matrix,
                       PRIMARY_COV){
  linfct <- list(t(contrast.matrix))
  names(linfct)[1] <- PRIMARY_COV
  # comment out library loading
  # library(multcomp)
  # limfct - a specification of the linear hypotheses to be tested.
  # Linear functions can be specified by either the matrix of
  # coefficients or by symbolic descriptions of one or more linear
  # hypotheses. Multiple comparisons in AN(C)OVA models are specified
  # by objects returned from function mcp
  t <- multcomp::glht(lm.obj, linfct = my_mcp(linfct))

  if (ncol(contrast.matrix) == 1) {
    tstat <- summary(t, multcomp::univariate())$test$tstat[1]
    df <- summary(t, multcomp::univariate())$df
    if (df > 0) {
      pvals <- 2 * pt(abs(tstat), df, lower.tail = FALSE)
    } else {
      pvals <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
    }
    return(data.frame(ID = NA,
                      logFC = summary(t, multcomp::univariate())$test$coefficients[1],
                      AveExpr = NA,
                      t = summary(t, multcomp::univariate())$test$tstat[1],
                      P.Value = pvals,
                      row.names = NULL))
  } else {
    # Ftest
    return(data.frame(ID = NA,
                      t(summary(t, test = multicomp::Ftest())$test$coefficients),
                      AveExpr = NA,
                      F = summary(t, test = multicomp::Ftest())$test$fstat[1, 1],
                      P.Value = summary(t, test = multicomp::Ftest())$test$pvalue[1, 1],
                      row.names = NULL))
  }
}
