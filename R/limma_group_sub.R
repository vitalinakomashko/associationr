
limma_group_sub <- function(indx_coef,
                            design,
                            contrast.matrix,
                            contrast.group,
                            NEED_FOR_PLOT,
                            PRIMARY_COV,
                            ADJUST_COVS,
                            BLOCK_COV,
                            GROUP_COV,
                            pheno,
                            temp_mat,
                            gene.ann,
                            OUTLOC,
                            B.num){
  # fit with linear model
  if (!is.null(BLOCK_COV)) {
    message("mixed effect modeling")
    # limma::duplicateCorrelation? ----
    # estimate the correlation between duplicate spots (regularly
    # spaced replicae spots on the same array. Why is it used here?
    corfit <- limma::duplicateCorrelation(temp_mat,
                                          design,
                                          block = pheno[[BLOCK_COV]])
    # are we looking here for "consensus.correlation" value?
    fit <- limma::lmFit(temp_mat,
                        design,
                        block = pheno[[BLOCK_COV]],
                        correlation = corfit$consensus)
  } else {
    # fit linear model for each gene give a series of arrays
    fit <- limma::lmFit(temp_mat,design)
  }
  # compute moderated t, F statistics and log odds of differential
  # expression
  fit <- limma::eBayes(fit)


  # get residual
  # back primary variable information
  if (colnames(design)[1] %in% c("(Intercept)","X.Intercept.")) {
    back_coef <- c(1,indx_coef)
  } else {
    back_coef <- indx_coef
  }
  temp_mat_adjusted <- residuals(fit,temp_mat)
  temp_mat_adjusted <- temp_mat_adjusted +
    t(design[, back_coef] %*% t(fit$coefficients[, back_coef]))
  # fitting with all
  fitted.value <- t(fit$design %*% t(fit$coefficients))
  # fitting with adjusted covs
  adjust.value <- t(fit$design[,-back_coef] %*%
                      t(fit$coefficients[, -back_coef]))
  # adjusting by adjust covs
  fitted.value_adjusted <- fitted.value - adjust.value

  # ATTENTION ----
  # function side effect - file system modification
  if (is.null(contrast.matrix)){
    # numeric test
    saveRDS(fit, file = paste0(OUTLOC,
                               "/fit_",
                               PRIMARY_COV,
                               ".rds"))
    write.table(design, file = paste0(OUTLOC,
                                      "/design_",
                                      PRIMARY_COV,
                                      ".txt"),
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
    # extract a table of the top-ranked genes from a linear
    # model fit
    DS <- limma::topTable(fit,
                          genelist = rownames(temp_mat),
                          number = fit$coefficients %>% nrow(),
                          p.value = 1,
                          lfc = log2(1),
                          coef = indx_coef,
                          adjust = "BH")
    suffix <- NULL
    p_group <- unique(pheno[[GROUP_COV]])

    # ATTENTION ----
    # limma_group_plot doesn't return anything
    limma_group_plot(DS,
                     temp_mat,
                     temp_mat_adjusted,
                     fitted.value,
                     fitted.value_adjusted,
                     gene.ann,
                     pheno,
                     PRIMARY_COV,
                     ADJUST_COVS,
                     BLOCK_COV,
                     GROUP_COV,
                     suffix,
                     p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)

    return(list(data = temp_mat,
                data.adjusted = temp_mat_adjusted,
                fitted = fitted.value,
                fitted.adjusted = fitted.value_adjusted))
  }

  # given a lienar model fit to microarray data, compute
  # estimated coefficients and standard errors for a
  # given set of contrasts
  fit.const <- limma::contrasts.fit(fit,contrast.matrix)
  fit.const <- limma::eBayes(fit.const) # this is basically ANOVA


  # ATTENTION ----
  # side effect - file system modification
  # save result
  saveRDS(contrast.matrix, file = paste0(OUTLOC,
                                         "/contrast.matrix_",
                                         PRIMARY_COV,
                                         ".rds"))
  saveRDS(fit.const, file = paste0(OUTLOC,
                                   "/fit_",
                                   PRIMARY_COV,
                                   ".rds"))
  write.table(design,
              file = paste0(OUTLOC,
                            "/design_",
                            PRIMARY_COV,
                            ".txt"),
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)


  if (ncol(contrast.matrix) == 1) {
    message("just two group comparision")
    suffix <- contrast.group
    DS <- limma::topTable(fit.const,
                          genelist = rownames(temp_mat),
                          number = fit.const$coefficients %>% nrow(),
                          p.value = 1,
                          lfc = log2(1),
                          coef = suffix,
                          adjust = "BH")

    p_group <- unique(unlist(strsplit(suffix, "-")))
    # ATTENTION ----
    # limma_group_plot doesn't return anything
    limma_group_plot(DS,
                     temp_mat,
                     temp_mat_adjusted,
                     fitted.value,
                     fitted.value_adjusted,
                     gene.ann,
                     pheno,
                     PRIMARY_COV,
                     ADJUST_COVS,
                     BLOCK_COV,
                     GROUP_COV,
                     suffix,
                     p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
  } else {
    message("do anova")
    suffix <- "ANOVA"
    DS <- limma::topTable(fit.const,
                          genelist = rownames(temp_mat),
                          number = fit.const$coefficients %>% nrow(),
                          p.value = 1,
                          lfc = log2(1),
                          coef = colnames(contrast.matrix),
                          adjust = "BH")
    p_group <- unique(unlist(strsplit(contrast.group, "-")))
    # ATTENTION ----
    # limma_group_plot doesn't return anything
    limma_group_plot(DS,
                     temp_mat,
                     temp_mat_adjusted,
                     fitted.value,
                     fitted.value_adjusted,
                     gene.ann,
                     pheno,
                     PRIMARY_COV,
                     ADJUST_COVS,
                     BLOCK_COV,
                     GROUP_COV,
                     suffix,
                     p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)

    message("do each of comparison")

    # this is a different suffix than used above in limma_group_plot
    for (suffix in unique(contrast.group)){

      DS <- limma::topTable(fit.const,
                            genelist = rownames(temp_mat),
                            number = fit.const$coefficients %>% nrow(),
                            p.value = 1, lfc = log2(1),
                            coef = colnames(contrast.matrix)[contrast.group == suffix],
                            adjust = "BH")
      p_group <- unique(unlist(strsplit(suffix, "-")))
      # ATTENTION ----
      # limma_group_plot doesn't return anything
      limma_group_plot(DS,
                       temp_mat,
                       temp_mat_adjusted,
                       fitted.value,
                       fitted.value_adjusted,
                       gene.ann,
                       pheno,
                       PRIMARY_COV,
                       ADJUST_COVS,
                       BLOCK_COV,
                       GROUP_COV,
                       suffix,
                       p_group,
                       NEED_FOR_PLOT,
                       OUTLOC)
    }
  }

  return(list(data = temp_mat,
              data.adjusted = temp_mat_adjusted,
              fitted = fitted.value,
              fitted.adjusted = fitted.value_adjusted))
}
