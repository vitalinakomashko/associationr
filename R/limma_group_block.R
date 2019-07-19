limma_group_block <- function(expdata,
                              PRIMARY_COV,
                              ADJUST_COVS,
                              BLOCK_COV,
                              INTERACT_COV = NULL,
                              spline_func = NULL,
                              OUTLOC,
                              B.num = 1000,
                              eBayes = TRUE,
                              SAVE_ADJUSTED = FALSE,
                              PARAM,
                              DiffVar = FALSE){
  # comment out library loading
  # library(limma)
  # library(splines)
  #
  # ATTENTION ----
  # function side effect: modification of the file system
  # make output directory
  dir.create(paste0(OUTLOC), recursive = TRUE, showWarnings = FALSE)

  # Variable should be there in data
  NEED_FOR_PLOT <- c(PRIMARY_COV, ADJUST_COVS, BLOCK_COV, INTERACT_COV)

  # ATTENTION ----
  # data checking should be taken out of the this function and possibly
  # main function
  # check data is ready
  expdata <- sanity_check(expdata, PRIMARY_COV, NEED_FOR_PLOT)

  # extract data components
  pheno <- expdata$sample.ann
  temp_mat <- expdata$data
  gene.ann <- expdata$gene.ann

  # spline_function
  # if spline_func argument is not provided the assign PRIMARY_COV to the
  # spline_func. What are the possible values for the spline function?
  # is it TRUE/FALSE or a name of a variable?
  if (is.null(spline_func)) {
    spline_func <- PRIMARY_COV
  } else {
    spline_func <- paste0(spline_func, "(", PRIMARY_COV, ")" )
  }


  # ATTENTION ----
  # design_obj can be possibly taken out of the function
  # limma_group_block
  # this will be important for limma only, because INTERACT_COV
  # is not used for lm
  if (!is.null(INTERACT_COV)) {
    design_obj <- make_design_obj_interaction(pheno,
                                              spline_func,
                                              PRIMARY_COV,
                                              INTERACT_COV,
                                              ADJUST_COVS)
    GROUP_COV <- INTERACT_COV
  } else {
    design_obj <- make_design_obj(pheno,
                                  spline_func,
                                  PRIMARY_COV,
                                  ADJUST_COVS)
    GROUP_COV <- PRIMARY_COV
  }


  # ATTENTION ----
  # PARAM$VOOM is not documented in the example
  #
  # check if PARAM$VOOM is present and set to TRUE
  if (!is.null(PARAM$VOOM) && PARAM$VOOM) {
    message("voom normalization")
    # Normalize using voom, WITH adjusting for known and HV covariates:
    temp_mat <- edgeR::calcNormFactors(edgeR::DGEList(counts = temp_mat),
                                       method = PARAM$calcNormFactors.method) %>%
      limma::voom(., design = design_obj$design,
                  normalize.method = PARAM$normalize.method,
                  plot = FALSE)
  }

  # ATTENTION ----
  # Shinya didn't provide the function diffvar_group_sub therefore
  # this block needs to be rewritten and DiffVar variable might need
  # to be taken out from the list of input parameters to limma_group_block.
  # run limma block
  if (DiffVar==TRUE) {
    temp_mat_adjusted <- diffvar_group_sub(design_obj$indx_coef,
                                           design_obj$design,
                                           design_obj$contrast.matrix,
                                           design_obj$contrast.group,
                                           NEED_FOR_PLOT,
                                           PRIMARY_COV,
                                           ADJUST_COVS,
                                           BLOCK_COV,
                                           GROUP_COV,
                                           pheno,
                                           temp_mat,
                                           gene.ann,
                                           OUTLOC)
  } else if (eBayes) {
    temp_mat_adjusted <- limma_group_sub(design_obj$indx_coef,
                                         design_obj$design,
                                         design_obj$contrast.matrix,
                                         design_obj$contrast.group,
                                         NEED_FOR_PLOT,
                                         PRIMARY_COV,
                                         ADJUST_COVS,
                                         BLOCK_COV,
                                         GROUP_COV,
                                         pheno,
                                         temp_mat,
                                         gene.ann,
                                         OUTLOC)
  } else {
    temp_mat_adjusted <- lm_group_sub(design_obj$indx_coef,
                                      design_obj$design,
                                      design_obj$contrast.matrix,
                                      design_obj$contrast.group,
                                      spline_func,
                                      NEED_FOR_PLOT,
                                      PRIMARY_COV,
                                      ADJUST_COVS,
                                      BLOCK_COV,
                                      GROUP_COV,
                                      pheno,
                                      temp_mat,
                                      gene.ann,
                                      OUTLOC)
  }

  # ATTENTION ----
  # limma_group_block returns TRUE and has a side effect of
  # modifying the file system. Should return the generated object.
  if (SAVE_ADJUSTED) {
    saveRDS(temp_mat_adjusted,
            file = paste0(OUTLOC,"/data_adjusted.rds"))
  }
  return(TRUE)
}
