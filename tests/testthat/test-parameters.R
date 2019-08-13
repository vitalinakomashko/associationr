test_that("verify length", {
  params <- list()
  expect_error(verify_parameters(params), "at least three parameters")
})


test_that("missing mandatory parameters", {
  params <- list(primary_covs = "Group", adjust_covs = "Cov1")
  expect_error(verify_parameters(params),
               "Expecting mandatory parameter")
  params <- list(primary_covs = "Group",
                 adjust_covs = "Cov1",
                 test_method = "NULL")
  expect_error(verify_parameters(params),
               "Expecting mandatory parameter")
})


test_that("remove NULL", {
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 something_else = "NULL")
  expect_equal(verify_parameters(params),
               list(primary_covs = "Group",
                    test_method = "limma",
                    adjust_covs = "Group1"))
})


test_that("verify test method", {
  params <- list(primary_covs = "Group",
                 test_method = "PCA",
                 adjust_covs = "Group1")
  expect_error(verify_parameters(params), "'limma' or 'lm'")
})


test_that("voom parameters", {
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 voom = list())
  expect_error(verify_parameters(params),
               "`norm_factors_method` or `voom_normalize_method` are not present.")
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 voom = list(norm_factors_method = "method",
                             voom_normalize_method = "method"))
  expect_error(verify_parameters(params), "Unexpected value")
})


test_that("verify file path", {
  file_path <- ""
  expect_error(read_params(file_path))
})


test_that("verify ignore_sample_size is TRUE or FALSE if provided", {
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 ignore_sample_size = "boop")
  expect_error(verify_parameters(params),
               "If parameter `ignore_sample_size` is provided it should be set to TRUE or FALSE.")
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 ignore_sample_size = TRUE)
  expect_equal(verify_parameters(params), params)
  params <- list(primary_covs = "Group",
                 test_method = "limma",
                 adjust_covs = "Group1",
                 ignore_sample_size = FALSE)
  expect_equal(verify_parameters(params), params)
})
