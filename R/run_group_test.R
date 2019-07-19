run_group_test <- function(PARAM, expdata){
  # comment out library loading
  # library(limma)
  # drop NA: find rows which don't have missing values in the primary phenotype
  indx <- !is.na(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]])
  # in the data, retain only the samples without missing values in the
  # primary phenotype
  expdata$data <- expdata$data[, indx]
  # in the sample annotation data frame also keep the samples without missing
  # values in the primary phenotype
  expdata$sample.ann <- expdata$sample.ann[indx, , drop=FALSE]

  # primary covariate has a variation?
  # stop execution if the primary phenotype has fewer or equal to one unique
  # value
  # ATTENTION ----
  # Replace message and return() with stop(message)
  if (length(unique(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]])) <= 1) {
    message("phenotype has no variation")
    return()
  }

  # ATTENTION ----
  # PARAM$ignore_sample_size is not documented in the example.
  # At this point the directory is not created, so doesn't need
  # to be unlink. Replace message and return() with stop(message), remove
  # unlink. Add "!" to message() in else.
  #
  # if the primary phenotype is a factor and the smallest number of samples
  # for a level is less than 2, the check if PARAME$ignore_sample_size is present
  # and has a value (TRUE or FALSE).
  if (is.factor(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]]) &&
      min(table(expdata$sample.ann[[PARAM[["PRIMARY_COVS"]]]])) < 2) {
    if(is.null(PARAM$ignore_sample_size) || !PARAM$ignore_sample_size) {
      message("minimum phenotype group should have 3 samples at least")
      unlink(list.files(PARAM$OUTLOC,full.names = T))
      return()
    } else {
      message("sample size less than 2")
    }
  }

  # ATTENTION ----
  # Replace message and return() with stop(message)
  # if number of sample is fewer than the number of covariates to adjust for
  # plus 2 then exit.
  if ((expdata$sample.ann%>%nrow) < (2+length(PARAM[["ADJUST_COVS"]]))) {
    message("sample size too small")
    return()
  }

  # ATTENTION ----
  # possibly need to check for NULL values for each. Also need to check
  # that it doesn't lead to an empty vector for ADJUST_COVS after removing
  # INTERACT_COV.
  #
  # removes "INTERACT_COV" from "ADJUST_COVS"
  # make interaction covariates and adjusted covariate unique
  PARAM[["ADJUST_COVS"]] <- setdiff(PARAM[["ADJUST_COVS"]],
                                    PARAM[["INTERACT_COV"]])


  # ATTENTION ----
  # limma_group_block sets parameters based on PARAM list, but also
  # passes the entire list as a separate parameter
  #
  # of the TEST.METHOD is LIMMA, the function is called with
  # eBayes == TRUE. The other difference might have to do with the
  # INTERACT_COV. It might be used in LIMMA and might not be used in LM,
  # but this needs to be confirmed.
  if (PARAM$TEST.METHOD == "LIMMA") {
    limma_group_block(expdata,
                      PRIMARY_COV = PARAM[["PRIMARY_COVS"]],
                      ADJUST_COVS = PARAM[["ADJUST_COVS"]],
                      BLOCK_COV = PARAM[["BLOCK_COV"]],
                      INTERACT_COV = PARAM[["INTERACT_COV"]],
                      spline_func = PARAM[["spline_func"]],
                      OUTLOC = PARAM[["OUTLOC"]],
                      eBayes = TRUE,
                      SAVE_ADJUSTED = PARAM[["SAVE_ADJUSTED"]],
                      PARAM = PARAM)
  } else {
    # ATTENTION ----
    # logic is not clear
    # if the TEST.METHOD is not LIMMA (assuming LM) and
    # INTERACT_COV is provided then the function exits.
    if (!is.null(PARAM[["INTERACT_COV"]])) {
      return()
    }
    limma_group_block(expdata,
                      PRIMARY_COV = PARAM[["PRIMARY_COVS"]],
                      ADJUST_COVS = PARAM[["ADJUST_COVS"]],
                      BLOCK_COV = PARAM[["BLOCK_COV"]],
                      INTERACT_COV = PARAM[["INTERACT_COV"]],
                      spline_func = PARAM[["spline_func"]],
                      OUTLOC = PARAM[["OUTLOC"]],
                      eBayes = FALSE,
                      SAVE_ADJUSTED = PARAM[["SAVE_ADJUSTED"]],
                      PARAM = PARAM)
  }
}
