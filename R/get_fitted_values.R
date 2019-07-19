get_fitted_values <- function(lm.obj,
                              back_coef){

  #fitting with all
  fitted.value <- t(matrix(fitted(lm.obj)))

  if (class(coef(lm.obj)) == "coef.mer") {
    #fitted.value_adjusted = rowSums(model.matrix(lm.obj)[,colnames(coef(lm.obj)[[1]])[back_coef]]*
    #                                  coef(lm.obj)[[1]][,back_coef])
    fitted.value_adjusted <- t(model.matrix(lm.obj)[ , names(lme4::fixef(lm.obj))[back_coef]] %*%
                                 lme4::fixef(lm.obj)[back_coef])
  }else{
    fitted.value_adjusted <- t(model.matrix(lm.obj)[, back_coef] %*% coef(lm.obj)[back_coef])
  }

  temp_mat_adjusted <- fitted.value_adjusted + residuals(lm.obj)

  return(list(raw.value_adjusted = temp_mat_adjusted,
              fitted.value = fitted.value,
              fitted.value_adjusted = fitted.value_adjusted))
}
