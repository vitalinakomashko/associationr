get_fitted_values <- function(lm.obj,
                              back_coef){

  #fitting with all
  fitted_value <- t(matrix(fitted(lm.obj)))

  if (class(coef(lm.obj)) == "coef.mer") {
    mat1 <- model.matrix(lm.obj)[, names(lme4::fixef(lm.obj))[back_coef]]
    mat2 <- lme4::fixef(lm.obj)[back_coef]
    fitted_value_adjusted <- t(mat1 %*% mat2)
  } else {
    mat1 <- model.matrix(lm.obj)[, back_coef]
    mat2 <- coef(lm.obj)[back_coef]
    fitted_value_adjusted <- t(mat1 %*% mat2)
  }

  temp_mat_adjusted <- fitted_value_adjusted + residuals(lm.obj)

  return(list(raw_value_adjusted = temp_mat_adjusted,
              fitted_value = fitted_value,
              fitted_value_adjusted = fitted_value_adjusted))
}
