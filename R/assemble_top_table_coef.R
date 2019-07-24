assemble_top_table_coef <- function(lm.obj,
                                    indx_coef){
  if (class(lm.obj) == "lmerMod") {
    K <- diag(length(coef(lm.obj)[[1]]))[indx_coef, , drop = FALSE]
    rownames(K) <- names(coef(lm.obj)[[1]])[indx_coef]
    t <- multcomp::glht(lm.obj, linfct = K)
  } else {
    K <- diag(length(coef(lm.obj)))[indx_coef, , drop = FALSE]
    rownames(K) <- names(coef(lm.obj))[indx_coef]
    t <- multcomp::glht(lm.obj, linfct = K)
  }

  if (length(indx_coef) == 1) {
    tstat <- summary(t, multcomp::univariate())$test$tstat[1]
    df <- summary(t, multcomp::univariate())$df
    if (df > 0) {
      pvals <- 2 * pt(abs(tstat), df, lower.tail = FALSE)
    } else {
      pvals <- 2 * pnorm(abs(tstat), lower.tail = FALSE)
    }
    return(
      data.frame(ID = NA,
                 logFC = summary(t, multcomp::univariate())$test$coefficients[1],
                 AveExpr = NA,
                 t = summary(t, multcomp::univariate())$test$tstat[1],
                 P.Value = pvals,
                 row.names = NULL)
    )
  } else {
    # Ftest
    return(
      data.frame(ID = NA,
                 t(summary(t, test = multicomp::Ftest())$test$coefficients),
                 AveExpr = NA,
                 F = summary(t, test = multicomp::Ftest())$test$fstat[1, 1],
                 P.Value = summary(t, test = multicomp::Ftest())$test$pvalue[1, 1],
                 row.names = NULL)
    )
  }
}
