# SCRIPT: Make dataframe of replicate-level results for analysis --------------------------------

#    this script uses the poparray objects and replicate info to 
#    build a data frame in which each unique replicate/scenario/clim is a row and key response and pred vars are columns

# KL: feel free to parallelize this if it might speed things up...

# KTS: changed to have only one version- with standalone switch turned OFF


STANDALONE=T

if(STANDALONE){
  
  # Clear workspace --------------------------------------------------------
  
  rm(list=ls())   # clear workspace
  
  # Load Functions --------------------------------------------------------
  
  source('GT_functions.R')   # load functions from file
  
    run_name<-"GT_06_04_23" #set the run name
  globalscratch=T
  if(globalscratch){
  poppath<-sprintf("/globalscratch/kjloope/%s/poparray",run_name)
  } else{poppath<-sprintf("%s/poparray",run_name)}
  # Load Packages -----------------------------------------------------------
  
  suppressWarnings(loadPackages())
  
  # Set up the Workspace for simulations -----------------------------
  # make sure reps and years match the actual simuations run
  
  params <- prepWorkspace(nreps=120,nyears = 94)
  
  # Read in list of scens and reps -------------------------------
  
  allruns <- read.csv(file=sprintf("%s/allruns.csv",run_name))
  
  # Prepare for making data frame of replicates ----------------------
  
  nreps <- max(allruns$rep)
  nyears <- params$NYEARS
  nclims <- length(params$clim_scenarios)
  
  allcmods <- params$clim_scenarios
  allscens <- params$allscenarios
  
}



allreps <- list()

 myCluster <- makeCluster(spec=60,type="PSOCK")    # PREPARE FOR RUNNING PARALLEL
  registerDoParallel(myCluster)
  
  #r=1
  allreps <- foreach(row=1:nrow(allruns),.packages = c('terra','raster','tryCatchLog','mvtnorm','popbio'))%dopar% {
#row = 1
#for(row in 1:nrow(allruns)){   # 
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
  
  temp <- rast(sprintf("setup_files/OutputTifs/%s_1.tif",params$clim_scenarios[1]))
  base_climvars <- temp[1]
  temp <- lapply(params$clim_var_names,function(t) base_climvars[[t]] <<- 0)
  base_climvars$baseT <- 20
  base_climvars$YSB <- 1
   
  basemat <- thismatfunc(base_climvars)
    
  # start building up dataframe
  temp2 <- data.frame(
    scenario = thisscen
  )
  temp2$rep <- thisrep
  temp2 <- cbind(temp2,thisscen_specs)
    
  # age at maturity (MA)  (not included since this will never change across replicates or scenario)
  # temp2$MA_base <- max(9,min(nstages,round(thispar$MA$ma.intercept.mn$this + thispar$MA$ma.tempeff.mn$this*base_climvars[["baseT"]]  ) ) )
  
  # PHI_A
  temp2$PHI_A_base <- thispar$PHI_A$adultsurv.mn$this
  temp2$PHI_A_burn <- thispar$PHI_A$phi.a.burneff.mn$this * thisscen_specs$phi.a_burn

  # PHI_J
  temp2$PHI_J_base <- thispar$PHI_J$juvsurv$this
  temp2$PHI_J_burn <- thispar$PHI_A$phi.a.burneff.mn$this * thisscen_specs$phi.j_burn
  
  # CS
  temp2$CS_base <- max(0.01,thispar$CS$cs.intercept.mn$this + thispar$CS$cs.tempeff.mn$this * base_climvars[["baseT"]] + thispar$CS$CS_intercept$this )
  temp2$CS_t <- thispar$CS$CS_junetmax$this*thisscen_specs$cs_clim
  
  # HS
  temp2$HS_base <- thispar$HS$hs.mn$this
  temp2$HS_t <- thispar$HS$HS_nest_tempeff$this*thisscen_specs$hs_clim
  temp2$HS_p <- thispar$HS$HS_nest_precipeff$this*thisscen_specs$hs_clim
  
  # NS
  temp2$NS_base <- thispar$NS$nestsucc$this
  
  # PF
  temp2$PF_base <- params$paraminfo$PF$PF_clim_eff$func(deltaTjune = 0,
                                                        deltaPjune = 0,
                                                        hs.tempeff.draw= 0,
                                                        hs.precipeff.draw = 0,
                                                        t.mean.ref = params$paraminfo$PF$t.mean.ref,
                                                        p.mean.ref = params$paraminfo$PF$p.mean.ref,
                                                        sex.tmean.pmean.top = params$paraminfo$PF$PF_clim_eff_model$obj)["perc.fem"] 
  temp2$PF_t <- thispar$PF$PF_nest_tempeff$this * thisscen_specs$pf_clim
  temp2$PF_p <- thispar$PF$PF_nest_precipeff$this * thisscen_specs$pf_clim
  
  # dispersal
  temp2$dispersal <- thispar$gamma$dispersal$this

  # PHI_H
  temp2$PHI_H_base <- max(0.001,min(0.99,thispar$PHI_H$hatchsurv$this))
  
  # PR
  temp2$PR_base <- plogis(thispar$PR$PR_intercept$this+1.8) 
  temp2$PR_t <- thispar$PR$PR_aprmaytmax$this * thisscen_specs$pr_clim
  
  # Matrix properties (really properties of the base matrix for this replicate, not reflective of the scenario)
  temp2$lambda = popbio::lambda(basemat)   # to compare lambda... not perfect, but provides some kind of comparison potential
  temp3 <- elasticity(basemat)
  temp2$elast_j1 <- temp3[2]
  temp2$elast_j8 <- temp3[params$NSTAGES*7+9]
  temp2$elast_a <- temp3[params$NSTAGES*params$NSTAGES]
  temp2$elast_f <- temp3[params$NSTAGES*(params$NSTAGES-1)+1]
  
  poparray <- readRDS(sprintf("%s/poparray2_scen_%s_%s_rep_%s.rds",poppath,thisscen,thisclim,thisrep)) # matrix: years, stages
  
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
  
  #allreps[[row]] <<- temp2
  return(temp2)
}

stopCluster(myCluster)
allrepsdf<-do.call(rbind,allreps)
#head(allrepsdf) 
saveRDS(allreps,file=sprintf("%s/allreps_results.rds",run_name ))
write.csv(allrepsdf , sprintf("%s/allreps_results.csv",run_name ) )



# End script --------------------


