# limma_group_plot definition ----
limma_group_plot <- function(DS,
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
                             OUTLOC){

  #output
  indx <- match(DS[, 1], rownames(gene.ann))
  DS <- cbind(DS, gene.ann[indx, ,drop = FALSE])
  colnames(DS)[1] <- "Name"
  tryCatch({
    #    write.table(DS,file=paste(OUTLOC,"/stat_",PRIMARY_COV,"_",suffix,".txt",sep=""),
    #                quote=FALSE,sep="\t",row.names=FALSE)
  }, error = function(e) {}
  )
  # ATTENTION ----
  # function side effect - file system modification
  saveRDS(DS,
          file = paste0(OUTLOC,
                        "/stat_",
                        PRIMARY_COV,
                        "_",
                        suffix,
                        ".rds"))
  tryCatch({
    # make pvalue histgram
    # ATTENTION ----
    # function side effect - file system modification
    pdf(paste0(OUTLOC,
               "/pvalue_dist_",
               PRIMARY_COV,
               "_",
               suffix,
               ".pdf"))
    hist(DS$P.Value)
    dev.off()
  }, warning = function(w) {
  }, error = function(e) {
  }, finally = {
  })

  # plot top 100
  indx <- DS$P.Value < 0.05
  indx[is.na(indx)] <- FALSE
  SIG_VARS <- DS[1:min(6, nrow(DS)), 1]
  #SIG_VARS = SIG_VARS[1:6]


  # plot raw data and adjusted data
  # also spline fit curve

  # take out data for plot
  # unadjusted
  if (class(temp_mat)=="EList") {
    temp_mat <- temp_mat$E
  }
  ds <- cbind(t(temp_mat[SIG_VARS, , drop = FALSE]),
              pheno[ , NEED_FOR_PLOT, drop = FALSE])
  ds <- reshape2::melt(ds, id.vars = NEED_FOR_PLOT)

  ds_fitted <- cbind(t(fitted.value[SIG_VARS, , drop=FALSE]),
                     pheno[ , NEED_FOR_PLOT, drop = FALSE])
  ds_fitted <- reshape2::melt(ds_fitted,
                              id.vars = NEED_FOR_PLOT,
                              value.name = "fit")
  # combine
  ds <- merge(ds,
              ds_fitted,
              by = c(NEED_FOR_PLOT, "variable"))
  ds$Type <- "raw"

  # adjusted
  ds_adjusted <- cbind(t(temp_mat_adjusted[SIG_VARS, , drop = FALSE]),
                       pheno[ , NEED_FOR_PLOT, drop = FALSE])
  ds_adjusted <- reshape2::melt(ds_adjusted,
                                id.vars = NEED_FOR_PLOT)

  ds_fitted <- cbind(t(fitted.value_adjusted[SIG_VARS, , drop = FALSE]),
                     pheno[ , NEED_FOR_PLOT, drop = FALSE])
  ds_fitted <- reshape2::melt(ds_fitted, id.vars = NEED_FOR_PLOT,
                              value.name = "fit")
  # combine
  ds_adjusted <- merge(ds_adjusted,
                       ds_fitted,
                       by = c(NEED_FOR_PLOT, "variable"))
  ds_adjusted$Type <- "adjusted"

  ds_adjusted$variable <- factor(as.character(ds_adjusted$variable),
                                 levels = SIG_VARS)
  ds$variable <- factor(as.character(ds$variable),
                        levels=SIG_VARS)
  # filter group
  ds <- ds[ds[[GROUP_COV]] %in% p_group, ]
  ds_adjusted <- ds_adjusted[ds_adjusted[[GROUP_COV]] %in% p_group, ]

  # comment out library loading here
  # library(ggthemes)
  if (is.numeric(ds[[PRIMARY_COV]]) | is.numeric(ds[[GROUP_COV]])) {
    g_obj <- list()
    g_obj[[1]] <- ggplot2::ggplot(ds_adjusted,
                                  ggplot2::aes_string(PRIMARY_COV,
                                                      "value",
                                                      group = BLOCK_COV,
                                                      color = GROUP_COV)) +
      ggplot2::facet_wrap(~ variable,
                          scales = "free_y",
                          ncol = 3) +
      ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(y = fit),
                         alpha = 0.3,
                         size = 2)
    g_obj[[2]] <- ggplot2::ggplot(ds,
                                  ggplot2::aes_string(PRIMARY_COV,
                                                      "value",
                                                      group = BLOCK_COV,
                                                      color = GROUP_COV)) +
      ggplot2::facet_wrap( ~ variable,
                           scales = "free_y",
                           ncol = 3) +
      ggplot2::geom_point() +
      ggplot2::geom_line(ggplot2::aes(y = fit),
                         alpha = 0.3,
                         size = 2)

    plot.nrow <- ceiling(length(unique(ds$variable)) / 3)
    plot.ncol <- min(3, length(unique(ds$variable)))

    pdf(paste0(OUTLOC,
               "/plot_",
               PRIMARY_COV,
               "_",
               suffix,
               ".pdf"),
        width = plot.ncol * 3 + 2,
        height = plot.nrow * 2 + 0.5)
    plot(g_obj[[1]])
    plot(g_obj[[2]])
    dev.off()
    return()
  }

  # for classes, plot two version
  g_obj = list()
  g_obj[[1]] <- ggplot2::ggplot(ds_adjusted,
                                ggplot2::aes_string(PRIMARY_COV,
                                                    "value",
                                                    fill = GROUP_COV)) +
    ggplot2::facet_wrap( ~ variable,
                         scales = "free_y",
                         ncol = 3) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggthemes::scale_fill_stata() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1,
                                                       colour = "grey10"))
  g_obj[[2]] <- ggplot2::ggplot(ds,
                                ggplot2::aes_string(PRIMARY_COV,
                                                    "value",
                                                    fill = GROUP_COV)) +
    ggplot2::facet_wrap( ~ variable,
                         scales = "free_y",
                         ncol = 3) +
    ggplot2::geom_boxplot() +
    ggplot2::theme_bw() +
    ggthemes::scale_fill_stata() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60,
                                                       hjust = 1,
                                                       colour = "grey10"))

  plot.nrow <- ceiling(length(unique(ds$variable)) / 3)
  plot.ncol <- min(3, length(unique(ds$variable)))
  n_groups <- length(p_group)

  pdf(paste0(OUTLOC,
             "/plot_",
             PRIMARY_COV,
             "_",
             suffix,
             ".pdf"),
      width = plot.ncol*(1 + n_groups * 0.5) + 2,
      height = plot.nrow * 2 + 0.5)
  plot(g_obj[[1]])
  plot(g_obj[[2]])
  dev.off()
}
