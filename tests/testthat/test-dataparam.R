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
               "`param_list` is expected, but not found")
  params <- list(primary_covs = "Group", adjust_covs = "Cov2")
  expect_error(verify_input_data_parameters(expdata, params),
               "`param_list` is expected, but not found")
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
               "Expected class")
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
                                   Group = factor(c(1, rep(NA, 9))),
                                   Cov1 = rnorm(10))
  expdata$gene_ann <- data.frame(Name = paste0("g", 1:10))
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1")
  expect_error(verify_input_data_parameters(expdata, params),
               "Primary phenotype has only 1 unique value")
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1",
                 ignore_sample_size = TRUE)
  expect_message(verify_input_data_parameters(expdata, params),
                 "Min number of samples per category in the primary covariate is fewer than 3.")
})

