test_that("param values in sample_ann columns", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1,5),rep(0,5))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "group")
  expect_error(verify_input_data_parameters(expdata, params),
               paste0("Values in the `param_list` not found ",
                      "in the column names"))
  params <- list(primary_covs = "Group", adjust_covs = c("Cov2", "Cov1"))
  expect_error(verify_input_data_parameters(expdata, params),
               paste0("Values in the `param_list` not found ",
                      "in the column names"))
})


test_that("column types", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = letters[1:10],
                                   Cov1 = rnorm(10),
                                   stringsAsFactors = FALSE)
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_error(verify_input_data_parameters(expdata, params),
               "Expected classes of columns")
})


test_that("remove samples with missing values", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1, 4),
                                                  rep(2, 4),
                                                  rep(NA, 2))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expdata$data <- expdata$data[, !is.na(expdata$sample_ann$Group)]
  expdata$sample_ann <- expdata$sample_ann[!is.na(expdata$sample_ann$Group), ]
  expect_equal(verify_input_data_parameters(expdata, params)$expdata$data,
               expdata$data)
})



test_that("need more than 1 unique value in the primary phenotype", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1, 8),
                                                    rep(NA, 2))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_error(verify_input_data_parameters(expdata, params),
               "Primary phenotype has only 1 unique value")
})


test_that("min number of samples", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1, 8), rep(2,2))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_error(verify_input_data_parameters(expdata, params),
                 paste0("Min number of samples per category in the primary covariate ",
                 "should be at least 3 or ignore_sample_size parameter should be ",
                 "set to TRUE."))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1",
                 ignore_sample_size = TRUE)
  expect_message(verify_input_data_parameters(expdata, params),
               paste0("Min number of samples per category in the primary ",
                      "covariate is fewer than 3."))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1",
                 ignore_sample_size = FALSE)
  expect_error(verify_input_data_parameters(expdata, params),
               paste0("Min number of samples per category in the primary covariate ",
                      "should be at least 3."))
})


test_that("factor primary covariate can have only 2 levels", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*20), nrow = 10, ncol = 20)
  colnames(expdata$data) <- paste0("s", 1:20)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:20),
                                   Group = factor(c(rep(1, 6),
                                                    rep(3, 6),
                                                    rep(4, 8))),
                                   Cov1 = rnorm(20))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_error(verify_input_data_parameters(expdata, params),
               paste0("Currently the package can handle a primary covariate with ",
                      "no more than two classes."))
})


test_that("limma compatible values in primary_covs", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1, 5), rep(2, 4), NA)),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_equal(verify_input_data_parameters(expdata, params)$expdata$sample_ann$Group,
                          factor(c(rep("X1", 5), rep("X2", 4))))
  expect_message(verify_input_data_parameters(expdata, params),
                 "not compatible with limma")
})


test_that("number of samples fewer than the number of adjust_covs", {
  expdata <- list()
  expdata$data <- matrix(rnorm(10*10), nrow = 10, ncol = 10)
  colnames(expdata$data) <- paste0("s", 1:10)
  rownames(expdata$data) <- paste0("g", 1:10)
  expdata$sample_ann <- data.frame(SAMPLE.ID = paste0("s",1:10),
                                   Group = factor(c(rep(1, 3), rep(2, 3),
                                                    rep(NA, 4))),
                                   Cov1 = rnorm(10),
                                   Cov2 = rnorm(10, 2, 0.3),
                                   Cov3 = rbinom(10, 2, 0.3),
                                   Cov4 = rchisq(10, 2),
                                   Cov5 = rnorm(10),
                                   Cov6 = factor(1:10),
                                   Cov7 = rexp(10),
                                   Cov8 = rgeom(10, 0.2),
                                   Age = sample(60:85, 10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = paste0("Cov", 1:6),
                 ignore_sample_size = TRUE)
  expect_error(verify_input_data_parameters(expdata, params),
               paste0("The number of samples is fewer than the number of ",
                      "adjustment variables plus 2."))
})
