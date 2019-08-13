#' Verify provided data and its structure.
#'
#' @param expdata Named list with the data components: data, gene_ann and
#' sample_ann.
#'
#' @return Named list with the data components: data, gene_ann and sample_ann.
#'
#' @export
verify_expdata <- function(expdata){

  # verify that expdata is a list
  if (class(expdata) != "list") {
    stop("Expecting parameter `expdata` to be a list.",
         call. = FALSE)
  }

  # verify that the number of elements is 3
  if (length(expdata) != 3) {
    stop("Expecting three elements in `expdata`.",
         call. = FALSE)
  }

  # verify presence of all elements
  needed <- c("data", "gene_ann", "sample_ann")
  for(i in needed){
    if (!(i %in% names(expdata))) {
      stop("Expecting mandatory element ", i, " in `expdata`, ",
           "however, did not find it.",
           call. = FALSE)
    }
  }

  if(!is.matrix(expdata$data)){
    stop("Expecting `data` element in `expdata` to be of class matrix.",
         call. = FALSE)
  }

  # verify that `data` has rownames
  if (is.null(rownames(expdata$data))) {
    stop("Expecting non-null row names in the element `data` ",
         "of `expdata`.",
         call. = FALSE)
  }

  # verify that `data` has column names
  if (is.null(colnames(expdata$data))) {
    stop("Expecting non-null column names in the element `data` ",
         "of `expdata`.",
         call. = FALSE)
  }

  # check class of gene_ann component
  if (!is.data.frame(expdata$gene_ann)) {
    stop("Element `gene_ann` in `expdata` expected to be an instance of class ",
         "data.frame.",
         call. = FALSE)
  }

  # verify that the column names in `gene_ann` are set
  if (is.null(colnames(expdata$gene_ann))){
    stop("Expecting non-null column names in the element `gene_ann` ",
         "of `expdata`.",
         call. = FALSE)
  }

  # check class of sample_ann component
  if (!is.data.frame(expdata$sample_ann)) {
    stop("Expecting element `sample_ann` in `expdata` to be an instance of ",
         "class data.frame.",
         call. = FALSE)
  }

  # verify that the column names in `gene_ann` are set
  if (is.null(colnames(expdata$sample_ann))){
    stop("Expecting non-null column names in the element `sample_ann` ",
         "of `expdata`.",
         call. = FALSE)
  }

  # verify number of columns in `data` == number of rows in `sample_ann`
  if (nrow(expdata$sample_ann) != ncol(expdata$data)) {
    stop("The number of columns in `data` should be equal to ",
         "the number of rows in `sample_ann` of `expdata`.",
         call. = FALSE)
  }

  # verify match of row names in `sample_ann` and column names in `data`
  if (all(rownames(expdata$sample_ann) != colnames(expdata$data))) {
    stop("The column names in `data` should match the row names ",
         "in `sample_ann` of `expdata`.",
         call. = FALSE)
  }

  # verify number of rows in `data` == number of rows in `gene_ann`
  if (nrow(expdata$data) != nrow(expdata$gene_ann)) {
    stop("The number of rows in `data` should be equal to ",
         "the number of rows `gene_ann` of `expdata`.",
         call. = FALSE)
  }

  # verify match of row names in `gene_ann` and row names in `data`
  if (all(rownames(expdata$gene_ann) != rownames(expdata$data))) {
    stop("The row names in the element `data` should match the row names ",
         "`gene_ann` of `expdata`.",
         call. = FALSE)
  }

  return(expdata)
}




#' Verify provided parameters vs the data provided.
#'
#' The function performs checks essential for running the analysis: presence
#' of for primary_covs, adjust_covs and others in the columns of
#' strong{expdata$sample_ann}; checks the classes of the columns in
#' strong{expdata$sample_ann} needed for analysis; removes all samples with the
#' missing values in the column with the primary covariate; check for the number
#' of unique values in the primary covariate; if primary covariate is a
#' categorical variable, then checks for the smallest number of samples per
#' category; compares number of samples and the number of adjustment variables;
#' adjust the spline function and creates group covariates depending on the
#' presence of adjustment and interaction covariates.
#'
#' @param expdata Named list with data.
#' @param param_list Named list with the analysis parameters.
#'
#' @return Named list with strong{expdata} and strong{param_list}.
#' @export

 verify_input_data_parameters <- function(expdata, param_list, verbose = TRUE) {

  # here check for plotting variables
  # verify that variables don't have missing values in sample_ann
  # checking NA
  # if (any(colSums(is.na(expdata$sample.ann[ , NEED_FOR_PLOT, drop = FALSE])) > 0)) {
  #   stop("having NAs")
  # }

  # checking pheno type existence
  # if (any(!(NEED_FOR_PLOT %in% colnames(expdata$sample.ann)))) {
  #   stop("cannot find nessesary phenotypes")
  # }

  # checking phenotype class
  # indx <- !sapply(expdata$sample.ann[ , NEED_FOR_PLOT,
  #                                     drop = FALSE], class) %in%
  #   c("factor","numeric","integer")
  # if (sum(indx) > 0) {
  #   stop("factor or numeric")
  # }


  # verify that the variable names in param_list are present in sample_ann
  possible_variables <- c("primary_covs", "adjust_covs", "interact_covs",
                          "block_covs", "spline_fun")

  provided_params <- param_list[names(param_list) %in% possible_variables]
  for(i in names(provided_params)){
    # all is for the cases where there are multiple values per key in the yaml
    if(!all(provided_params[[i]] %in% colnames(expdata$sample_ann))){
      stop("Values in the `param_list` not found ",
           "in the column names of `expdata$sample_ann`.",
           call. = FALSE)
    }
  }


  # verify their type: should be numeric, factor or integer
  for(i in names(provided_params)){
    k <- provided_params[[i]]
    if(!(is.integer(expdata$sample_ann[[k]]) |
       is.numeric(expdata$sample_ann[[k]]) |
       is.factor(expdata$sample_ann[[k]]))){
      stop("Expected class of ", i, " in `expdata$sample_ann` is integer, ",
           "numeric or factor.",
           call. = FALSE)
    }
  }

  primary_covs <- param_list[["primary_covs"]]

  # find samples which don't have missing values in the main phenotype
  # primary_covs we subset the data to have only those samples.
  indx <- !is.na(expdata$sample_ann[[primary_covs]])
  # in the data, retain only the samples without missing values in the
  # primary phenotype
  expdata$data <- expdata$data[, indx]
  # in the sample annotation data frame also keep the samples without missing
  # values in the primary phenotype
  expdata$sample_ann <- expdata$sample_ann[indx, , drop = FALSE]



  # require that the main phenotype (primary_covs) has more than 1 unique value
  primary_cov_unique_values <- length(unique(expdata$sample_ann[[primary_covs]]))
  if (primary_cov_unique_values == 1) {
    stop("Primary phenotype has only 1 unique value, need at least 2.",
         call. = FALSE)

  }


  # if the primary phenotype is a categorical variable and the smallest number of
  # samples in per category is less than 2 then we check the parameter
  # ignore_sample_size.
  # if ignore_sample_size is absent or FALSE then the function exits. Otherwise,
  # provide a message that the sample size is small, but proceed.
  if (is.factor(expdata$sample_ann[[primary_covs]])) {
    min_number_samples <- expdata$sample_ann %>%
      dplyr::group_by_at(primary_covs) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::slice(1)
    if (min_number_samples$n <= 3){
      if("ignore_sample_size" %in% names(param_list)){
        if(param_list$ignore_sample_size){
          if(verbose){
            message("Min number of samples per category in the primary ",
                    "covariate is fewer than 3.")
            }
        } else {
          stop("Min number of samples per category in the primary covariate ",
               "should be at least 3.",
               call. = FALSE)

        }
      } else {
        stop("Min number of samples per category in the primary covariate ",
             "should be at least 3 or ignore_sample_size parameter should be ",
             "set to TRUE.",
             call. = FALSE)
      }
    }
  }

  # verify that values in the primary_covs variable are compatible with limma;
  # if not then create valid names using make.names()
  if (is.factor(expdata$sample_ann[[primary_covs]])) {
    primary_cov_valid_names <- make.names(
      as.character(
        expdata$sample_ann[[primary_covs]]
        )
      )
    if (!all(primary_cov_valid_names ==
             as.character(expdata$sample_ann[[primary_covs]]))) {
      if (verbose) {
        message("Character vectors in the column ", primary_covs,
                "of `expdata$sample_ann` are not compatible with limma.",
                "Will use `make.names()` to create syntactically valid names.")
      }
      expdata$sample_ann[[primary_covs]] <-
        factor(make.names(expdata$sample_ann[[primary_covs]]))
    }
  }


  # if the number of samples in the data is fewer than the number of adjustment
  # parameters plus 2 then exit the calculations.
  if (nrow(expdata$sample_ann) < (length(param_list[["adjust_covs"]]) + 2)) {
    stop("The number of samples is fewer than the number of ",
         "adjustment variables plus 2.",
         call. = FALSE)
  }

  # remove interaction covariates from the adjustment covariates
  if ("interact_covs" %in% names(param_list)) {
    param_list[["adjust_covs"]] <- setdiff(param_list[["adjust_covs"]],
                                           param_list[["interact_covs"]])
  }


  # adjust the spline_fun parameter
  if (!("spline_fun" %in% names(param_list))) {
    param_list[["spline_fun"]] <- param_list[["primary_covs"]]
  } else {
    param_list[["spline_fun"]] <- paste0(param_list[["spline_fun"]],
                                         "(",
                                         param_list[["primary_covs"]],
                                         ")")
  }


  # populate group_covs in param_list
  if ("interact_covs" %in% names(param_list)) {
    param_list[["group_covs"]] <- param_list[["interact_covs"]]
  } else {
    param_list[["group_covs"]] <- param_list[["primary_covs"]]
  }

  return(list(expdata = expdata, param_list = param_list))
}
