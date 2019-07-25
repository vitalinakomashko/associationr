# lm_group_sub definition ----
lm_group_sub <- function(indx_coef,
                         design,
                         contrast.matrix,
                         contrast.group,
                         spline_func,
                         NEED_FOR_PLOT,
                         PRIMARY_COV,
                         ADJUST_COVS,
                         BLOCK_COV,
                         GROUP_COV,
                         pheno,
                         temp_mat,
                         gene.ann,
                         OUTLOC,
                         B.num){
  # comment out library loading
  # library(lme4)
  orig_name <- rownames(temp_mat)
  rownames(temp_mat) <- paste0("X", 1:nrow(temp_mat)) # get fomula can handle
  #fit_ds = cbind(pheno[,c(PRIMARY_COV,ADJUST_COVS,BLOCK_COV),drop=FALSE],t(temp_mat))
  pheno_use <- pheno[ , c(PRIMARY_COV, ADJUST_COVS, BLOCK_COV), drop=FALSE]
  mat_use <- t(temp_mat)
  rownames(temp_mat) <- orig_name # back to original ones


  # ATTENTION ----
  # pretty significant chunk of code, may want to take it out and create a
  # separate function
  DS.ordinary <- lapply(1:nrow(temp_mat), function(j){
    #print(j)
    y <- paste0("X", j)
    tmp.data <- cbind(mat_use[ , j, drop = FALSE], pheno_use) %>%
      #dplyr::select_(fit_ds,.dots=c(y,PRIMARY_COV,ADJUST_COVS,BLOCK_COV))%>%
      dplyr::mutate(r.n. = rownames(.)) %>%
      dplyr::filter(rowSums(is.na(.)) == 0)

    indx <- apply(tmp.data, 2, function(x) length(unique(x))) == 1
    ADJUST_COVS_rm <- colnames(tmp.data)[indx]


    if (length(unique(tmp.data[!is.na(tmp.data[,1]),2])) == 1 |
        length(unique(tmp.data[!is.na(tmp.data[,1]),1])) < 3) {
      return(NULL)
    }

    if (is.null(contrast.matrix)) {
      # numeric test
      suffix <- "quant"

      if (is.null(BLOCK_COV)) {
        lm.obj <- lm(as.formula(paste0(y,paste(c("~1",
                                                 spline_func,
                                                 setdiff(ADJUST_COVS,
                                                         ADJUST_COVS_rm)),
                                               collapse = "+"))),
                     data = tmp.data,
                     x = TRUE)
        indx_coef <- which(regexpr(PRIMARY_COV,
                                   names(lm.obj$coefficients),
                                   fixed = TRUE) > 0)
      } else {
        lm.obj <- lme4::lmer(as.formula(paste0(y,
                                               paste(c("~1",
                                                       spline_func,
                                                       setdiff(ADJUST_COVS,
                                                               ADJUST_COVS_rm),
                                                       paste0("(1|",
                                                              BLOCK_COV,")")),
                                                     collapse = "+"))),
                             data = tmp.data,
                             REML = FALSE)
        indx_coef <- which(regexpr(PRIMARY_COV,
                                   names(coef(lm.obj)[[1]]),
                                   fixed = TRUE) > 0)
      }
      DS.ordinary <- list()
      DS.ordinary[[suffix]] <- assemble_top_table_coef(lm.obj, indx_coef)
      DS.ordinary[[suffix]]$ID <- orig_name[j]
      DS.ordinary[[suffix]]$AveExpr <- mean(tmp.data[[y]])

      # get residual
      back_coef <- c(1, indx_coef)
      fitted <- get_fitted_values(lm.obj, back_coef)
      colnames(fitted$raw_value_adjusted) <- tmp.data$r.n.
      colnames(fitted$fitted_value) <- tmp.data$r.n.
      colnames(fitted$fitted_value_adjusted) <- tmp.data$r.n.
      return(list(DS.ordinary = DS.ordinary,
                  fitted = fitted))
    }

    if (is.null(BLOCK_COV)) {
      lm.obj <- lm(as.formula(paste0(y,
                                     paste(c("~0",
                                             spline_func,
                                             setdiff(ADJUST_COVS,
                                                     ADJUST_COVS_rm)),
                                           collapse = "+"))),
                   data = tmp.data,
                   x = TRUE)
      indx_coef <- which(regexpr(PRIMARY_COV,
                                 names(lm.obj$coefficients),
                                 fixed = TRUE) > 0)
    } else {
      lm.obj <- lme4::lmer(as.formula(paste0(y,
                                             paste(c("~0",
                                                     spline_func,
                                                     setdiff(ADJUST_COVS,
                                                             ADJUST_COVS_rm),
                                                     paste0("(1|",
                                                            BLOCK_COV,")")),
                                                   collapse = "+"))),
                           data = tmp.data,
                           REML = FALSE)
      indx_coef <- which(regexpr(PRIMARY_COV,
                                 names(lme4::fixef(lm.obj)),
                                 fixed = TRUE) > 0)
    }


    fitted <- get_fitted_values(lm.obj,indx_coef)
    colnames(fitted$raw.value_adjusted) <- tmp.data$r.n.
    colnames(fitted$fitted.value) <- tmp.data$r.n.
    colnames(fitted$fitted.value_adjusted) <- tmp.data$r.n.

    DS.ordinary <- list()
    if (ncol(contrast.matrix) == 1) {
      suffix <- contrast.group
      # single contrast
      DS.ordinary[[suffix]] <- assemble_top_table(lm.obj,
                                          contrast.matrix[indx_coef, , drop = FALSE],
                                          PRIMARY_COV)
      DS.ordinary[[suffix]]$ID <- orig_name[j]
      DS.ordinary[[suffix]]$AveExpr <- mean(tmp.data[[y]])
    } else {
      #ANOVA
      suffix <- "ANOVA"
      DS.ordinary[[suffix]] <- assemble_top_table(lm.obj,
                                          contrast.matrix[indx_coef, , drop = FALSE],
                                          PRIMARY_COV)
      DS.ordinary[[suffix]]$ID <- orig_name[j]
      DS.ordinary[[suffix]]$AveExpr <- mean(tmp.data[[y]])
      for (suffix in unique(contrast.group)) {
        # each contrast
        DS.ordinary[[suffix]] <- assemble_top_table(lm.obj,
                                            contrast.matrix[indx_coef,
                                                            colnames(contrast.matrix)[contrast.group== suffix],
                                                            drop=FALSE],
                                            PRIMARY_COV)
        DS.ordinary[[suffix]]$ID <- orig_name[j]
        DS.ordinary[[suffix]]$AveExpr <- mean(tmp.data[[y]])
      }
    }
    return(list(DS.ordinary = DS.ordinary,
                fitted = fitted))
  }) # end of lapply call.

  names(DS.ordinary) <- orig_name
  indx <- !sapply(DS.ordinary, is.null) #remove NA
  DS.ordinary <- DS.ordinary[indx]
  temp_mat <- temp_mat[indx,]

  # recover nessesary information for plotting
  temp_mat_adjusted <- lapply(1 : length(DS.ordinary),
                              function(x){
                                data.frame(SAMPLE.ID =
                                             colnames(DS.ordinary[[x]]$fitted$raw.value_adjusted),
                                           ID = names(DS.ordinary)[x],
                                           value = t(DS.ordinary[[x]]$fitted$raw.value_adjusted))
                              })
  temp_mat_adjusted <- do.call(rbind, temp_mat_adjusted)
  temp_mat_adjusted <- dplyr::group_by(temp_mat_adjusted, SAMPLE.ID, ID) %>%
    tidyr::spread(SAMPLE.ID,value)
  r.n. <- temp_mat_adjusted$ID
  temp_mat_adjusted <- as.matrix(temp_mat_adjusted[,-1])
  rownames(temp_mat_adjusted) <- r.n.
  indx <- match(rownames(temp_mat), rownames(temp_mat_adjusted))
  temp_mat_adjusted <- temp_mat_adjusted[indx, ]
  indx <- match(colnames(temp_mat),colnames(temp_mat_adjusted))
  temp_mat_adjusted <- temp_mat_adjusted[,indx]

  fitted.value <- lapply(1:length(DS.ordinary),
                         function(x){
                           data.frame(SAMPLE.ID = colnames(DS.ordinary[[x]]$fitted$fitted.value),
                                      ID = names(DS.ordinary)[x],
                                      value = t(DS.ordinary[[x]]$fitted$fitted.value))
                         })
  fitted.value <- do.call(rbind,fitted.value)
  fitted.value <- dplyr::group_by(fitted.value, SAMPLE.ID, ID) %>%
    tidyr::spread(SAMPLE.ID,value)
  r.n. <- fitted.value$ID
  fitted.value <- as.matrix(fitted.value[ , -1])
  rownames(fitted.value) <- r.n.
  indx <- match(rownames(temp_mat), rownames(fitted.value))
  fitted.value <- fitted.value[indx, ]
  indx <- match(colnames(temp_mat), colnames(fitted.value))
  fitted.value <- fitted.value[ , indx]

  fitted.value_adjusted <- lapply(1 : length(DS.ordinary),
                                  function(x){
                                    data.frame(SAMPLE.ID = colnames(DS.ordinary[[x]]$fitted$fitted.value_adjusted),
                                               ID = names(DS.ordinary)[x],
                                               value = t(DS.ordinary[[x]]$fitted$fitted.value_adjusted))
                                  })
  fitted.value_adjusted <- do.call(rbind,fitted.value_adjusted)
  fitted.value_adjusted <- dplyr::group_by(fitted.value_adjusted, SAMPLE.ID, ID) %>%
    tidyr::spread(SAMPLE.ID, value)
  r.n. <- fitted.value_adjusted$ID
  fitted.value_adjusted <- as.matrix(fitted.value_adjusted[ , -1])
  rownames(fitted.value_adjusted) <- r.n.
  indx <- match(rownames(temp_mat), rownames(fitted.value_adjusted))
  fitted.value_adjusted <- fitted.value_adjusted[indx, ]
  indx <- match(colnames(temp_mat), colnames(fitted.value_adjusted))
  fitted.value_adjusted <- fitted.value_adjusted[ , indx]

  for (suffix in names(DS.ordinary[[1]]$DS.ordinary)){
    DS <- lapply(DS.ordinary, function(x){
      return(x$DS.ordinary[[suffix]])
    })
    DS <- do.call(rbind,DS) %>% dplyr::arrange(P.Value)

    if (suffix == "ANOVA") {
      p_group <- unique(unlist(strsplit(contrast.group, "-")))
    } else if (suffix == "quant") {
      suffix <- NULL
      p_group <- unique(pheno[[GROUP_COV]])
    } else {
      p_group <- unique(unlist(strsplit(suffix, "-")))
    }
    # ATTENTION ----
    # limma_group_plot doesn't return anything
    limma_group_plot(DS,
                     temp_mat,
                     temp_mat_adjusted,
                     fitted.value,
                     fitted.value_adjusted,
                     gene.ann,
                     pheno,
                     PRIMARY_COV,
                     ADJUST_COVS,
                     BLOCK_COV,
                     GROUP_COV,
                     suffix,
                     p_group,
                     NEED_FOR_PLOT,
                     OUTLOC)
  }
  return(list(data = temp_mat,
              data.adjusted = temp_mat_adjusted,
              fitted = fitted.value,
              fitted.adjusted = fitted.value_adjusted))
}

