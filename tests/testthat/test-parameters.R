test_that("verify length", {
  params <- list()
  expect_error(verify_parameters(params))
  params <- list(primary_covs = "Group")
  expect_error(verify_parameters(params))
})


test_that("missing mandatory parameters", {
  params <- list(primary_covs = "Group", adjust_covs = "Cov1")
  expect_error(verify_parameters(params))
})


test_that("remove NULL", {
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "NULL")
  expect_equal(verify_parameters(params),
               list(primary_covs = "Group",
                    test_method = "limma"))
})


test_that("verify test method", {
  params <- list(primary_covs = "Group",
                 test_method = "PCA",
                 adjust_covs = "NULL")
  expect_error(verify_parameters(params))
})


test_that("voom parameters", {
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "NULL",
                 voom = TRUE)
  expect_error(verify_parameters(params))
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "NULL",
                 voom = TRUE,
                 norm_factors_method = "method")
  expect_error(verify_parameters(params))
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "NULL",
                 voom = TRUE,
                 voom_normalize_method = "method")
  expect_error(verify_parameters(params))
})


test_that("verify file path", {
  file_path <- ""
  expect_error(read_params(file_path))
})
