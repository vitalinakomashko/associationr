test_that("it is a list", {
  expdata <- data.frame(matrix(rnorm(10*10), nrow = 10, ncol = 10))
  expect_error(verify_expdata(mtcars), "list")
})


test_that("check for the number of elements", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expect_error(verify_expdata(expdata), "three elements")
})


test_that("check for mandatory elements", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$sample.ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "mandatory element")
  expdata <- list()
  expdata$dat <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "mandatory element")
  expdata <- list()
  expdata$dat <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "mandatory element")
})


test_that("check that data is a matrix or DGElist", {
  # data is a data frame instead of a matrix
  expdata <- list()
  expdata$data <- data.frame(matrix(rnorm(10*10), nrow = 10, ncol = 10))
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "matrix")

})


test_that("check for row names in data", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "non-null row names")
})


test_that("check for column names in data matrix", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata), "non-null column names")
})


test_that("gene_ann is a data frame", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- list()
  expect_error(verify_expdata(expdata),
               "Element `gene_ann` in `expdata` expected to be an instance of class data.frame.")
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expect_error(verify_expdata(expdata),
               "Element `gene_ann` in `expdata` expected to be an instance of class data.frame.")
})


test_that("gene_ann has column names", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  colnames(expdata$gene_ann) <- NULL
  expect_error(verify_expdata(expdata),
               "non-null column names in the element `gene_ann`")
})

test_that("sample_ann is a data frame", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- list()
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata),
               "Expecting element `sample_ann` in `expdata` to be an instance of class data.frame.")
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata),
               "Expecting element `sample_ann` in `expdata` to be an instance of class data.frame.")
})

test_that("sample_ann has column names", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  colnames(expdata$sample_ann) <- NULL
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  expect_error(verify_expdata(expdata),
               "Expecting non-null column names in the element `sample_ann`")
})

test_that("nrow(expdata$sample_ann) == ncol(expdata$data)", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*20), nrow = 10, ncol = 20)
  rownames(expdata$data) <- paste0("g", 1:10)
  colnames(expdata$data) <- paste0("s", 1:20)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s", 1:10),
                                   Group = factor(c(rep(1,5), rep(0,5))),
                                   Cov1 = rnorm(10))
  rownames(expdata$sample_ann) <- expdata$sample_ann$SAMPLE.ID
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  rownames(expdata$gene_ann) <- expdata$gene_ann$Name
  expect_error(verify_expdata(expdata),
               "The number of columns in `data` should be equal to the number of rows in `sample_ann` of `expdata`")
})

test_that("rownames(expdata$sample_ann) == colnames(expdata$data)", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  rownames(expdata$data) <- paste0("g", 1:10)
  colnames(expdata$data) <- paste0("sample", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s", 1:10),
                                   Group = factor(c(rep(1,5), rep(0,5))),
                                   Cov1 = rnorm(10))
  rownames(expdata$sample_ann) <- expdata$sample_ann$SAMPLE.ID
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  rownames(expdata$gene_ann) <- expdata$gene_ann$Name
  expect_error(verify_expdata(expdata),
               "The column names in `data` should match the row names in `sample_ann` of `expdata`.")
})


test_that("nrow(expdata$data) == nrow(expdata$gene_ann)", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*20), nrow = 20, ncol = 10)
  rownames(expdata$data) <- paste0("gene", 1:20)
  colnames(expdata$data) <- paste0("s", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s", 1:10),
                                   Group = factor(c(rep(1,5), rep(0,5))),
                                   Cov1 = rnorm(10))
  rownames(expdata$sample_ann) <- expdata$sample_ann$SAMPLE.ID
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  rownames(expdata$gene_ann) <- expdata$gene_ann$Name
  expect_error(verify_expdata(expdata),
               "The number of rows in `data` should be equal to the number of rows `gene_ann` of `expdata`.")
})

test_that("rownames(expdata$data) == rownames(expdata$gene_ann)", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  rownames(expdata$data) <- paste0("protein", 1:10)
  colnames(expdata$data) <- paste0("s", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s", 1:10),
                                   Group = factor(c(rep(1,5), rep(0,5))),
                                   Cov1 = rnorm(10))
  rownames(expdata$sample_ann) <- expdata$sample_ann$SAMPLE.ID
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  rownames(expdata$gene_ann) <- expdata$gene_ann$Name
  expect_error(verify_expdata(expdata),
               "The row names in the element `data` should match the row names `gene_ann` of `expdata`.")
})
