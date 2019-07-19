check_expdata_easy <- function(expdata){
  # check expression data format
  # nessesary
  needed <- c("data","gene.ann","sample.ann")

  indx <- needed %in% names(expdata)
  if (!all(indx)) {
    stop(paste0("NOT found ",
                paste0(needed[!indx], collapse=" "),
                sep = " "))
  }

  # check class
  if (!(is.matrix(expdata$data) |
        class(expdata$data) == "EList" |
        class(expdata$data) == "DGEList")) {
    stop("expression data should be a matrix or EList")
  }

  if (!is.data.frame(expdata$gene.ann)) {
    stop("gene annotation data should be a data.frame")
  }

  if (!is.data.frame(expdata$sample.ann)) {
    stop("sample annotation data should be a data.frame")
  }

  # check sample number
  nsample <- dim(expdata$data)[2]
  if (dim(expdata$sample.ann)[1] != nsample) {
    stop("sample size does not match")
  }

  # check sample name
  sampleid <- colnames(expdata$data)
  if (!all(rownames(expdata$sample.ann) == sampleid)) {
    stop("sample name not match")
  }


  # check gene number
  ngene <- dim(expdata$data)[1]
  if (dim(expdata$gene.ann)[1] != ngene) {
    stop("gene size does not match")
  }

  # check gene name
  geneid = rownames(expdata$data)
  if (!all(rownames(expdata$gene.ann) == geneid)) {
    stop("gene name does not match")
  }

  return(TRUE)
}


sanity_check <- function(expdata,
                         PRIMARY_COV,
                         NEED_FOR_PLOT) {
  # checking expdata
  if (!check_expdata_easy(expdata)){
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
