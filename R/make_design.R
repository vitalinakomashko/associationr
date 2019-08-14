# Functions related to creation of design objects


#' Create a design object
#' @param sample_annotation Data frame with the sample annotation
#' @return design object
#' @export

create_design_object <- function(sample_annotation, param_list){
  if ("interact_covs" %in% names(param_list)) {
    design_obj <- make_design_obj_interaction(sample_annotation,
                                              spline_func,
                                              PRIMARY_COV,
                                              INTERACT_COV,
                                              ADJUST_COVS)
    # this is already done in verify_input_data_parameters
    # GROUP_COV <- INTERACT_COV
  } else {
    design_obj <- make_design_obj(sample_annotation,
                                  spline_func,
                                  PRIMARY_COV,
                                  ADJUST_COVS)
    # this is already done in verify_input_data_parameters
    # GROUP_COV <- PRIMARY_COV
  }
  return(design_obj)
}



#' Create design with interactions
#' @param sample_annotation Data frame with the sample annotation
#' @param param_list Named list with the analysis parameters
#' @return design object
#' @export

make_design_obj_interaction <- function(sample_annotation,
                                        param_list){
  # make design matrix
  primary_covs <- param_list[["primary_covs"]]
  spline_fun <- param_list[["spline_fun"]]
  interact_covs <- param_list[["interact_covs"]]
  adjust_covs <- param_list[["adjust_covs"]]
  if (is.numeric(sample_annotation[[primary_covs]])) {
    design_form <- paste0("~1+",
                          paste(c(paste0(spline_fun, "*",
                                         interact_covs), adjust_covs),
                                collapse = "+"))
    design <- model.matrix(as.formula(design_form), data = sample_annotation)
    #colnames(design)
    indx_coef <- which(regexpr(spline_fun,
                               colnames(design),
                               fixed = TRUE) > 0 &
                         regexpr(":",
                                 colnames(design),
                                 fixed = TRUE) > 0 &
                         regexpr(interact_covs,
                                 colnames(design),
                                 fixed = TRUE) > 0)
    if (is.numeric(sample_annotation[[interact_covs]]) | length(indx_coef) == 1) {
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


#' Create design without interactions
#' @param sample_annotation Data frame with the sample annotation
#' @param
#' @return design object
#' @export
make_design_obj <- function(sample_annotation,
                            ){
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


#' Make interaction cnst (?)
make_interaction_cnst <- function(pheno,
                                  design,
                                  indx_coef,
                                  INTERACT_COV){

  # make contrast matrix
  ###
  cnst <- colnames(design)[indx_coef]
  cnst_group <- rep(NA,length(cnst))

  i_levels <- levels(pheno[[INTERACT_COV]])
  cnst_group <- do.call(rbind, strsplit(cnst, ":"))[,2]
  cnst_group <- paste0(cnst_group,
                       "-",
                       paste0(INTERACT_COV, i_levels[1]))

  cnstA <- make.names(cnst)
  cnst_groupA <- cnst_group

  ###
  itable <- do.call(rbind,
                    strsplit(colnames(design)[indx_coef],
                             ":"))

  u_plevel <- unique(itable[ , 1])
  p_position <- which(itable[ , 1] %in% u_plevel[1])
  cnst <- apply(t(combn(p_position, 2)), 1, function(x){
    paste0(make.names(colnames(design)[indx_coef[x]]),
           collapse = "-")
  })
  cnst <- sapply(1:length(u_plevel), function(x){
    gsub(u_plevel[1], u_plevel[x], cnst, fixed = TRUE)
  })

  cnst_group <- cnst
  for (x in make.names(unique(itable[,1]))){
    cnst_group <- gsub(x,"", cnst_group, fixed = TRUE)
  }
  cnst_group <- gsub("^\\.", "", cnst_group)
  cnst_group <- gsub("-\\.", "-", cnst_group)
  cnst <- c(cnst,cnstA)
  cnst_group <- c(cnst_group, cnst_groupA)
  cnst_group <- gsub(INTERACT_COV, "", cnst_group)
  return(list(cnst = cnst,
              cnst_group = cnst_group))
}
