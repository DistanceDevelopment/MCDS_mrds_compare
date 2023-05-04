test.models <- function(dist.data, truncation, transect = "line",
                        cutpoints = NULL, convert.units){
  # Fit all the following models to the data passed in:
  #  - uniform with 0, 1, 2, 3 cosine adjustments
  #  - half-normal with 0, 1, 2 cosine adjustments
  #  - half-normal with 0, 1, 2 Hermite polynomial
  #  - hazard rate with 0, 1, 2 simple polynomial adjustments
  models <- data.frame(key = c(rep("unif", 4),
                               rep("hn", 6),
                               rep("hr", 3)),
                       adj = c(rep("cos",7),
                               rep("herm", 3),
                               rep("poly",3)),
                       nadj = c(0,1,2,3, rep(0:2,3)))
  
  # Results storage
  lnl_R <- lnl_MCDS <- optimizer <- p_R <- p_MCDS <- Nhat_R <- Nhat_MCDS <- NULL
  
  for(i in seq(along = models$key)){
    # Fit model using only R optimizer
    fit_R <- try(ds(dist.data,
                truncation = truncation,
                transect = transect,
                formula = ~1,
                key = models$key[i],
                adjustment = models$adj[i],
                nadj = models$nadj[i],
                cutpoints = cutpoints,
                optimizer = "R"))
    # Store values
    if(!is(fit_R, "try-error")){
    lnl_R[i] <- fit_R$ddf$lnl
    p_R[i] <- length(fit_R$ddf$ds$aux$ddfobj$xmat$distance)/fit_R$ddf$Nhat
    Nhat_R[i] <- fit_R$ddf$Nhat
    }else{
      lnl_R[i] <- p_R[i] <- Nhat_R[i] <- NA
    }
    
    # Fit model using only MCDS optimizer
    fit_MCDS <- try(ds(dist.data,
                   truncation = truncation,
                   transect = transect,
                   formula = ~1,
                   key = models$key[i],
                   adjustment = models$adj[i],
                   nadj = models$nadj[i],
                   cutpoints = cutpoints,
                   optimizer = "MCDS"))
    # Store values
    if(!is(fit_MCDS, "try-error")){
      lnl_MCDS[i] <- fit_MCDS$ddf$lnl
      p_MCDS[i] <- length(fit_MCDS$ddf$ds$aux$ddfobj$xmat$distance)/fit_MCDS$ddf$Nhat
      Nhat_MCDS[i] <- fit_MCDS$ddf$Nhat  
    }else{
      lnl_MCDS[i] <- p_MCDS[i] <- Nhat_MCDS[i] <- NA
    }
    
    
    # Fit model using both
    fit_both <- try(ds(dist.data,
                   truncation = truncation,
                   transect = transect,
                   formula = ~1,
                   key = models$key[i],
                   adjustment = models$adj[i],
                   nadj = models$nadj[i],
                   cutpoints = cutpoints,
                   optimizer = "both"))
    # Store values
    if(!is(fit_both, "try-error")){
      optimizer[i] <- fit_both$ddf$optimise
    }else{
      optimizer[i] <- NA
    }
  }
  
  # Add results 
  models$lnl_R <- lnl_R
  models$lnl_MCDS <- lnl_MCDS
  models$optimizer <- optimizer
  models$p_R <- round(p_R,2)
  models$p_MCDS <- round(p_MCDS,2)
  models$Nhat_R <- round(Nhat_R,2)
  models$Nhat_MCDS <- round(Nhat_MCDS,2)

  return(models)
}


test.cov.models <- function(dist.data, truncation, transect = "line",
                        cutpoints = NULL, convert.units, models){
  # Fit all models that are passed in
  
  # Results storage
  lnl_R <- lnl_MCDS <- optimizer <- p_R <- p_MCDS <- Nhat_R <- Nhat_MCDS <- NULL
  
  for(i in seq(along = models)){
    # Fit model using only R optmizer
    fit_R <- try(ds(dist.data,
                    truncation = truncation,
                    transect = transect,
                    formula = models[[i]],
                    cutpoints = cutpoints,
                    optimizer = "R"))
    # Store values
    if(!is(fit_R, "try-error")){
      lnl_R[i] <- fit_R$ddf$lnl
      p_R[i] <- length(fit_R$ddf$ds$aux$ddfobj$xmat$distance)/fit_R$ddf$Nhat
      Nhat_R[i] <- fit_R$ddf$Nhat
    }else{
      lnl_R[i] <- p_R[i] <- Nhat_R[i] <- NA
    }
    
    # Fit model using only MCDS optimizer
    fit_MCDS <- try(ds(dist.data,
                       truncation = truncation,
                       transect = transect,
                       formula = models[[i]],
                       cutpoints = cutpoints,
                       optimizer = "MCDS"))
    # Store values
    if(!is(fit_MCDS, "try-error")){
      lnl_MCDS[i] <- fit_MCDS$ddf$lnl
      p_MCDS[i] <- length(fit_MCDS$ddf$ds$aux$ddfobj$xmat$distance)/fit_MCDS$ddf$Nhat
      Nhat_MCDS[i] <- fit_MCDS$ddf$Nhat  
    }else{
      lnl_MCDS[i] <- p_MCDS[i] <- Nhat_MCDS[i] <- NA
    }
    
    
    # Fit model using both
    fit_both <- try(ds(dist.data,
                       truncation = truncation,
                       transect = transect,
                       formula = models[[i]],
                       cutpoints = cutpoints,
                       optimizer = "both"))
    # Store values
    if(!is(fit_both, "try-error")){
      optimizer[i] <- fit_both$ddf$optimise
    }else{
      optimizer[i] <- NA
    }
  }
  
  # Add results
  results <- data.frame(models = as.character(models))
  results$lnl_R <- lnl_R
  results$lnl_MCDS <- lnl_MCDS
  results$optimizer <- optimizer
  results$p_R <- round(p_R,2)
  results$p_MCDS <- round(p_MCDS,2)
  results$Nhat_R <- round(Nhat_R,2)
  results$Nhat_MCDS <- round(Nhat_MCDS,2)
  
  return(results)
}