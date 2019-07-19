make_design_obj_interaction <- function(pheno,
                                        spline_func,
                                        PRIMARY_COV,
                                        INTERACT_COV,
                                        ADJUST_COVS) {
  # make design matrix
  if (is.numeric(pheno[[PRIMARY_COV]])) {
    design <- model.matrix(as.formula(paste0("~1+",
                                             paste(c(paste0(spline_func,
                                                            "*",
                                                            INTERACT_COV),
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    #colnames(design)
    indx_coef <- which(regexpr(spline_func,
                               colnames(design),
                               fixed = TRUE) > 0 &
                         regexpr(":",
                                 colnames(design),
                                 fixed = TRUE) > 0 &
                         regexpr(INTERACT_COV,
                                 colnames(design),
                                 fixed = TRUE) > 0)
    if (is.numeric(pheno[[INTERACT_COV]]) | length(indx_coef) == 1) {
      colnames(design) <- make.names(colnames(design))
      contrast.matrix <- NULL
      contrast.group <- NULL
    } else {
      cnst.obj <- make_interaction_cnst(pheno,
                                        design,
                                        indx_coef,
                                        INTERACT_COV)
      colnames(design) <- make.names(colnames(design))
      contrast.matrix <- limma::makeContrasts(contrasts = cnst.obj$cnst,
                                              levels = design)
      contrast.group <- cnst.obj$cnst_group
    }
  } else {
    p_cof_levels <- levels(pheno[[PRIMARY_COV]])
    if (length(p_cof_levels) != 2) {
      stop("now only can handle two classes")
    }
    design <- model.matrix(as.formula(paste0("~0+",
                                             paste(c(paste0(spline_func,
                                                            "*",
                                                            INTERACT_COV),
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef <- which(regexpr(spline_func,
                               colnames(design),
                               fixed = TRUE) > 0 &
                         regexpr(":",
                                 colnames(design),
                                 fixed = TRUE) > 0 &
                         regexpr(INTERACT_COV,
                                 colnames(design),
                                 fixed = TRUE) > 0)
    colnames(design) <- gsub(PRIMARY_COV, "", colnames(design))

    if (is.numeric(pheno[[INTERACT_COV]])) {
      colnames(design) <- make.names(colnames(design))
      contrast.matrix <- NULL
      contrast.group <- NULL
    } else {
      colnames(design) <- gsub(INTERACT_COV, "", colnames(design))
      colnames(design) <- make.names(colnames(design))

      if (length(indx_coef) == 1) {
        contrast.matrix <- NULL
        contrast.group <- NULL
      } else {
        # make constraints
        cnst <- colnames(design)[indx_coef]
        for (ii in 1:(length(indx_coef) - 1)){
          for (jj in (ii + 1):length(indx_coef)){
            cnst <- c(cnst, paste(colnames(design)[indx_coef[jj]],
                                  "-",
                                  colnames(design)[indx_coef[ii]],
                                  sep=""))
          }
        }

        contrast.matrix <- limma::makeContrasts(contrasts = cnst,
                                                levels = design)
        colname_contrast.matrix <- gsub(
          paste0(p_cof_levels[length(p_cof_levels)], "." ),
          "",
          colnames(contrast.matrix))
        colname_contrast.matrix[1 : length(indx_coef)] <-
          paste0(colname_contrast.matrix[1 : length(indx_coef)],
                 "-",
                 levels(pheno[[INTERACT_COV]])[1])
        colnames(contrast.matrix) <- colname_contrast.matrix
        contrast.group <- colnames(contrast.matrix)
      }
    }
  }

  colnames(design) <- gsub(INTERACT_COV, "", colnames(design))
  colnames(contrast.matrix) <- gsub(INTERACT_COV, "", colnames(contrast.matrix))
  rownames(contrast.matrix) <- gsub(INTERACT_COV, "", rownames(contrast.matrix))
  #contrast.group = gsub(INTERACT_COV,"",contrast.group)
  return(list(indx_coef = indx_coef,
              design = design,
              contrast.matrix = contrast.matrix,
              contrast.group = contrast.group))
}
