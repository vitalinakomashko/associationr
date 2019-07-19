make_design_obj <- function(pheno,
                            spline_func,
                            PRIMARY_COV,
                            ADJUST_COVS){

  # make design matrix
  if (is.numeric(pheno[[PRIMARY_COV]])) {
    design <- model.matrix(as.formula(paste0("~1+",
                                             paste(c(spline_func,
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef = which(regexpr(PRIMARY_COV,colnames(design),fixed = T)>0)

    contrast.matrix <- NULL
    contrast.group <- NULL
  } else {
    design <- model.matrix(as.formula(paste0("~0+",
                                             paste(c(spline_func,
                                                     ADJUST_COVS),
                                                   collapse = "+"))),
                           data = pheno)
    indx_coef <- which(regexpr(PRIMARY_COV,
                               colnames(design),
                               fixed = TRUE) > 0)
    colnames(design) <- gsub(PRIMARY_COV, "", colnames(design))

    # make constraints
    cnst <- c()
    for (ii in 1:(length(indx_coef)-1)){
      for (jj in (ii+1):length(indx_coef)){
        cnst <- c(cnst,
                  paste0(colnames(design)[indx_coef[jj]],
                         "-",
                         colnames(design)[indx_coef[ii]]))
      }
    }

    #print(cnst)
    contrast.matrix <- limma::makeContrasts(contrasts = cnst, levels = design)
    contrast.group <- colnames(contrast.matrix)
  }
  return(list(indx_coef = indx_coef,
              design = design,
              contrast.matrix = contrast.matrix,
              contrast.group = contrast.group))
}
