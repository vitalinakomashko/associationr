make_interaction_cnst <- function(pheno,
                                  design,
                                  indx_coef,
                                  INTERACT_COV){

  # make contrast matrix
  ###
  cnst <- colnames(design)[indx_coef]
  cnst_group <- rep(NA,length(cnst))

  i_levels <- levels(pheno[[INTERACT_COV]])
  cnst_group <- do.call(rbind, strsplit(cnst, ":"))[,2]
  cnst_group <- paste0(cnst_group,
                       "-",
                       paste0(INTERACT_COV, i_levels[1]))

  cnstA <- make.names(cnst)
  cnst_groupA <- cnst_group

  ###
  itable <- do.call(rbind,
                    strsplit(colnames(design)[indx_coef],
                             ":"))

  u_plevel <- unique(itable[ , 1])
  p_position <- which(itable[ , 1] %in% u_plevel[1])
  cnst <- apply(t(combn(p_position, 2)), 1, function(x){
    paste0(make.names(colnames(design)[indx_coef[x]]),
           collapse = "-")
  })
  cnst <- sapply(1:length(u_plevel), function(x){
    gsub(u_plevel[1], u_plevel[x], cnst, fixed = TRUE)
  })

  cnst_group <- cnst
  for (x in make.names(unique(itable[,1]))){
    cnst_group <- gsub(x,"", cnst_group, fixed = TRUE)
  }
  cnst_group <- gsub("^\\.", "", cnst_group)
  cnst_group <- gsub("-\\.", "-", cnst_group)
  cnst <- c(cnst,cnstA)
  cnst_group <- c(cnst_group, cnst_groupA)
  cnst_group <- gsub(INTERACT_COV, "", cnst_group)
  return(list(cnst = cnst,
              cnst_group = cnst_group))
}
