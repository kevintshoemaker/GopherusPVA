# SCRIPT: Make dataframe of replicate-level results for analysis --------------------------------

#   NOTE: KL changed this so that it is sourced in DT_RunSimulation.R Does not need to be run separately...
#    this script uses the poparray objects and replicate info to 
#    build a data frame in which each unique replicate/scenario/clim is a row and key response and pred vars are columns

# KL: feel free to parallelize this if it might speed things up...

allrepsdf <- NULL
row = 1
for(row in 1:nrow(allruns)){   # 
  thisscen<-allruns$scen[row] 
  thisrep <- allruns$rep[row]
  thisscen_name <- paste0("scen_",thisscen)
 
  thisscen_specs <- params$allscenarios[thisscen,]
  
  thisclim <- thisscen_specs$clim_scenario
  
  thispar <- readRDS(sprintf("%s/replist_%s.rds",run_name, thisrep))
  
  thismatfunc <- GenerateMatrixSampler(thispar = thispar, nstages=params$NSTAGES, 
                                         stagenames = names(params$INITABUND),
                                         thisdat = params$paraminfo,
                                         thisscen = thisscen_specs
  )
  
  temp <- rast("setup_files/OutputTifs/rcp26_1.tif")
  base_climvars <- temp[1]
  temp <- lapply(params$clim_var_names,function(t) base_climvars[[t]] <<- 0)
  base_climvars$spateff <- 0
   
  basemat <- thismatfunc(base_climvars,0)
    
  # start building up dataframe
  temp2 <- data.frame(
    scenario = thisscen
  )
  temp2$rep <- thisrep
  temp2 <- cbind(temp2,thisscen_specs)
    
  # GROWTH (and size at maturity)
  temp2$GRk_int <- thispar$ANC$ANC_growth.a.int$this
  temp2$GRk_siteprec <- thispar$ANC$ANC_growth.a.siteprec$this
  temp2$GRk_tmax  <- thispar$ANC$ANC_growth.k.tmax$this
  temp2$GRk_p  <- thispar$ANC$ANC_growth.k.precip$this
  temp2$GRk_tmin  <- thispar$ANC$ANC_growth.k.tmin$this
  temp2$GRa_int <- thispar$ANC$ANC_growth.k.int$this
  temp2$GRa_siteprec <- thispar$ANC$ANC_growth.k.siteprec$this
  temp2$GRa_tmax  <- thispar$ANC$ANC_growth.a.tmax$this
  temp2$GRa_p  <- thispar$ANC$ANC_growth.a.precip$this
  temp2$GRa_tmin  <- thispar$ANC$ANC_growth.a.tmin$this
  temp2$GRt0 <- thispar$ANC$ANC_growth.t0$this
  temp2$SAM <- thispar$MA$size.at.maturity$this
  
  # PHI_A
  temp2$PHI_A_base <- plogis(thispar$PHI_A$PHI_A_intercept$this+0.5)
  temp2$PHI_A_p <- thispar$PHI_A$PHI_A_precip$this * thisscen_specs$phi.a_clim
  temp2$PHI_A_t <- thispar$PHI_A$PHI_A_tmax$this * thisscen_specs$phi.a_clim
  temp2$PHI_A_syeff <- thispar$PHI_A$PHI_A_siteyear.sd$this  # magnitude of random effect term

  # PHI_J
  temp2$PHI_J_base <- max(0.5,min(0.99,thispar$PHI_J$juvsurv$this))
  
  # CS
  temp2$CS_base <- thispar$CS$CS_intercept$this
  temp2$CS_t <- thispar$CS$CS_springtmean$this*thisscen_specs$cs_clim
  temp2$CS_bs <- thispar$CS$CS_bodysize$this
  temp2$CS_p <- thispar$CS$CS_precip$this*thisscen_specs$cs_clim
  
  # HS
  temp2$HS_base <- thispar$HS$hs.mn$this
  temp2$HS_t <- thispar$HS$nest.hs.slope$this * thispar$HS$HS_nest_tempeff$this * thisscen_specs$hs_clim   # effect of nest temp on HS
  #temp2$nestT <- thispar$HS$HS_nest_tempeff$this * thisscen_specs$hs_clim  # effect of ambient T on nest T
  
  # NS
  temp2$NS_base <- max(0.001,min(.99, thispar$NS$nestsucc$this))
  
  # PF
  temp2$PF_t  <- thispar$HS$HS_nest_tempeff$this * thispar$PF$sex.tmean$this * thisscen_specs$pf_clim
  # NOTE: baseline is set at 0.5
  
  # dispersal
  temp2$dispersal <- thispar$gamma$dispersal$this

  
  # PHI_H
  temp2$PHI_H_base <- max(0.001,min(.999, thispar$PHI_H$hatchsurv$this))
  
  # PR
  temp2$PR_base <- max(0.001,min(.999, plogis(thispar$PR$PR_intercept$this + thispar$PR$PR_priorrep$this*0.9 )))
  temp2$PR_p  <- thispar$PR$PR_precip$this * thisscen_specs$pr_clim
  temp2$PR_bs  <- thispar$PR$PR_bodysize$this
  temp2$PR_t  <- thispar$PR$PR_temprange$this * thisscen_specs$pr_clim
  
  # Matrix properties (really properties of the base matrix for this replicate, not reflective of the scenario)
  temp2$lambda = popbio::lambda(basemat)   # to compare lambda... not perfect, but provides some kind of comparison potential
  temp3 <- elasticity(basemat)
  temp2$elast_j1 <- temp3[2]
  temp2$elast_j8 <- temp3[params$NSTAGES*7+9]
  temp2$elast_a <- temp3[params$NSTAGES*params$NSTAGES]
  temp2$elast_f <- temp3[params$NSTAGES*(params$NSTAGES-1)+1]
  
  poparray <- readRDS(sprintf("%s/poparray_scen_%s_%s_rep_%s.rds",poppath,thisscen,thisclim,thisrep)) # matrix: years, stages
  
  thisabunds <- poparray 
  thistotabunds <- apply(thisabunds,1,sum)
  temp2$final_abund = thistotabunds[nyears+1]
  temp2$final_abund2 <- mean(thistotabunds[(nyears-9):(nyears+1)])  # KTS: added mean of last 10 years
  temp2$init_abund <- thistotabunds[1]
  temp2$init_abund2 <- mean(thistotabunds[1:10])
  temp2$min_abund <- min(thistotabunds)  
  temp2$delta_abund <- thistotabunds[nyears+1]-thistotabunds[1]   # positive number is growth, negative is decline
  temp2$delta_abund2 <- temp2$final_abund2-temp2$init_abund2
  temp2$frac_change <- temp2$delta_abund/thistotabunds[1]
  temp2$frac_change2 <- temp2$delta_abund2/temp2$init_abund2
  temp2$final_adult <- thisabunds[nyears+1,params$NSTAGES]
  temp2$final_adult2 <- mean(thisabunds[(nyears-9):(nyears+1),params$NSTAGES])
  temp2$final_hatch <- thisabunds[nyears+1,1]
  temp2$final_hatch2 <- mean(thisabunds[(nyears-9):(nyears+1),1])
  temp2$final_agerat <- sum(thisabunds[nyears+1,(1:params$NSTAGES-1)])/temp2$final_adult
  
  allrepsdf <- rbind(allrepsdf,temp3)
}

head(allrepsdf) 

write.csv(allrepsdf , sprintf("%s/allreps_results.csv",run_name ) )



# End script --------------------


