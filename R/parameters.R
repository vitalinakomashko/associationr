#' Read analysis parameters from a YAML file
#'
#' @param file_path Character string with a path to a YAML file.
#'
#' @return Named list with parameters for analysis.
#'
#' @export

read_params <- function(file_path){
  if (!file.exists(file_path)) {
    stop("Attempting to read a yaml file. File `file_path` is not found.",
         call. = FALSE)
  } else {
    params <- yaml::read_yaml(file_path)
    # verify minimum required parameters
    params <- verify_parameters(params)
    return(params)
  }
}



#' Verify validity of provided analysis parameters
#'
#' \code{verify_parameters} checks for the length of the list,
#' presence of the mandatory parameters (\strong{primary_covs} and
#' \strong{test_method}), allowed test method values, removes NULL parameters
#' and checks "partner" parameters (for \strong{voom} it looks for
#' \strong{voom_calc_norm_factors_method} and \strong{voom_normalize_method}).
#'
#' @param param_list List with parameters.
#'
#' @return List with parameters.

verify_parameters <- function(param_list){

  # check list length
  if (length(param_list) == 0) {
    stop("Expecting at least three parameters for analysis in `param_list`, ",
         "however, zero provided.",
         call. = FALSE)
  }

  # check if any of the parameters have been set to NULL and remove them
  if ("NULL" %in% param_list) {
    param_list[[which(param_list == "NULL")]] <- NULL
  }

  # check mandatory
  mandatory <- c("primary_covs", "test_method", "adjust_covs")
  for(i in mandatory){
    if (!(i %in% names(param_list))) {
      stop("Expecting mandatory parameter ", i, " in `param_list`, ",
           "however, did not find it.",
           call. = FALSE)
    }
  }

  # check that test_method is set to either limma or lm
  test_method_allowed <- c("limma", "lm")
  if (param_list$test_method != "limma" & param_list$test_method != "lm") {
    stop("Parameter `test_method` in `param_list` can be 'limma' or 'lm', ",
         "however, found ", param_list$test_method, ".",
         call. = FALSE)
  }

  # if voom is present verify that norm_factors_method and
  # voom_normalize_method are also present
  if ("voom" %in% names(param_list)) {
    if (!(all(c("norm_factors_method", "voom_normalize_method") %in%
              names(param_list$voom)))) {
      stop("Parameter `voom` is present, however, ",
           "either `norm_factors_method` or `voom_normalize_method` ",
           "are not present.",
           call. = FALSE)
    }
    norm_factors_method <- c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
    if (!(any(norm_factors_method %in% param_list$voom$norm_factors_method))) {
      stop("Unexpected value for 'norm_factor_method' in `param_list`, ",
           "should be one of the following: 'TMM', 'TMMwsp', 'RLE', ",
           "'upperquartile' or 'none'. See `edgeR::calcNormFactors()` ",
           "for documentation.",
           call. = FALSE)
    }
    voom_normalize_method <- c("none", "scale", "quantile", "cyclicloess")
    if(!any(voom_normalize_method %in% param_list$voom$voom_normalize_method)) {
      stop("Unexpected value for 'voom_normalize_method' in `param_list`, ",
           "should be one of the following: 'none', 'scale', 'quantile' or ",
           "'cyclicloess'. See `limma::voom()` for documentation.",
           call. = FALSE)
    }
  }
  return(param_list)
}
