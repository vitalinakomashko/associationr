#' Verify provided data and its structure.
#'
#' @param expdata Named list with the data components: data, gene_ann and
#' sample_ann.
#'
#' @return Named list with the data components: data, gene_ann and sample_ann.
#'
#' @export
check_expdata <- function(expdata){

  # verify that expdata is a list
  if (class(expdata) != "list") {
    stop("Expecting parameter `expdata` to be a list.",
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

  # check class of data component
  if (!(is.matrix(expdata$data) |
        class(expdata$data) == "EList" |
        class(expdata$data) == "DGEList")) {
    stop("Expecting element `data` in `expdata` to be an instance of a class ",
         "matrix or EList, or DGEList.",
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
    stop("Expecting the number of columns in the element `data` and ",
         "the number of rows in the element `sample_ann` of `expdata` to be ",
         "the same.",
         call. = FALSE)
  }

  # verify match of row names in `sample_ann` and column names in `data`
  if (all(rownames(expdata$sample_ann) != colnames(expdata$data))) {
    stop("Expecting column names in the element `data` and ",
         "row names in the element `sample_ann` of `expdata` to match.",
         call. = FALSE)
  }

  # verify number of rows in `data` == number of rows in `gene_ann`
  if (nrow(expdata$data) != nrow(expdata$gene_ann)) {
    stop("Expecting the number of rows in the element `data` and ",
         "the number of rows in the element `gene_ann` of `expdata` to be ",
         "the same.",
         call. = FALSE)
  }

  # verify match of row names in `gene_ann` and row names in `data`
  if (all(rownames(expdata$gene_ann) != rownames(expdata$data))) {
    stop("Expecting row names in the element `data` and ",
         "row names in the element `gene_ann` of `expdata` to match.",
         call. = FALSE)
  }

  return(expdata)
}




#' Verify
sanity_check <- function(expdata,
                         PRIMARY_COV,
                         NEED_FOR_PLOT) {
  # checking expdata
  if (!check_expdata(expdata)){
    stop()
  }

  # checking NA
  if (any(colSums(is.na(expdata$sample.ann[ , NEED_FOR_PLOT, drop = FALSE])) > 0)) {
    stop("having NAs")
  }

  # checking pheno type existence
  if (any(!(NEED_FOR_PLOT %in% colnames(expdata$sample.ann)))) {
    stop("cannot find nessesary phenotypes")
  }

  # checking phenotype class
  indx <- !sapply(expdata$sample.ann[ , NEED_FOR_PLOT, drop = FALSE], class) %in%
    c("factor","numeric","integer")
  if (sum(indx) > 0) {
    stop("factor or numeric")
  }

  # redefine factor levels
  for (i in 1:ncol(expdata$sample.ann)) {
    if (is.factor(expdata$sample.ann[ , i])) {
      expdata$sample.ann[ , i] <- factor(expdata$sample.ann[ , i])
    }
  }

  if (!is.numeric(expdata$sample.ann[[PRIMARY_COV]])) {
    if (all.equal(make.names(expdata$sample.ann[[PRIMARY_COV]]),
                  as.character(expdata$sample.ann[[PRIMARY_COV]])) != TRUE){

      print("variable name is not compatible with limma")
      print("renaming...")
      expdata$sample.ann[[PRIMARY_COV]] =
        make.names(expdata$sample.ann[[PRIMARY_COV]]) %>% as.factor()
      # stop()
    }
  }
  return(expdata)
}
