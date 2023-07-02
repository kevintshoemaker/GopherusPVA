
# DESERT TORTOISE SIMULATION FUNCTIONS ----------

# note: many functions are shared with analogous "GT_functions_*.R" source. If key 
#   updates are made, the corresponding update should also be made to the GT script
#   and vice versa!


# FUNCTION Explore_CHB  ----------------------------------
#    this script returns a data frame useful for visualizing critical habitat breadth
#    for a given scenario and replicate

# params=params;haspop=haspop;thispar=thispar;scen=scen;rep=rep;thisyears=1:89;thisclim=thisclim;sample=250
Explore_CHB <- function(params,haspop,thispar,scen,rep,thisyears=1:89,thisclim,sample=250,yeareff=-1.5){ 
  
  allpars <- names(thispar)
  thisscen <- params$allscenarios[scen,]
  
  spateff <- terra::rast(params$zerorast)
  
  thismatfunc <- GenerateMatrixSampler(thispar = thispar, nstages=params$NSTAGES, 
                                       stagenames = params$ST_NAMES,
                                       thisdat = params$paraminfo,
                                       thisscen = thisscen,vitals=T) 
  
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$growth_clim) clim_vars_used <- c(clim_vars_used,"feb_may_tmax_","nov_apr_pr_","ann_tmin_")
  if(thisscen$phi.a_clim) clim_vars_used <- c(clim_vars_used,"sepx2_mar_pr_","mar_nov_tmax_")
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"spr_tmean_","mar_mar_pr_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_aug_tmean_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"july_tmean_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"mar_mar_pr_","mar_apr_mtr_")
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used,"spateff"))
  
  toret=NULL
  
  yr=1
  for(yr in thisyears){ 
    # yeareff = -2   # no need to account for good years and bad years here...
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",thisclim,yr) ) 
    thisyr_brick <- c(thisyr_brick,spateff)
    names(thisyr_brick)[terra::nlyr(thisyr_brick)] <- "spateff"
    
    thissamp <- sample(haspop,sample,replace=T)
    
    toret2 <- thisyr_brick[thissamp]
    temp <- sapply(clim_vars_unused,function(t) toret2[[t]] <<- 0   )
    
    VR_list <- lapply(1:nrow(toret2),function(t) as.data.frame(thismatfunc(climvars=toret2[t,],yeareff=yeareff)) )
    VR_list <- do.call(rbind,VR_list)
    toret2 <- cbind(toret2,VR_list) 
    
    toret2$year = yr
    toret <- rbind(toret,toret2)
    
  }
  # toret
  forPCA <- toret[,clim_vars_used]
  # pca <- princomp(forPCA)
  
  # Maximum Likelihood Factor Analysis
  # entering raw data and extracting 3 factors,
  # with varimax rotation
  fit <- factanal(forPCA, 3, rotation="varimax",scores="regression")
  print(fit, digits=2, cutoff=.3, sort=TRUE)
  # plot factor 1 by factor 2
  load <- fit$loadings[,c(1,2,3)]
  
  toret <- cbind(toret,fit$scores)
  
  # plot(load,type="n") # set up plot
  # text(load,labels=names(forPCA),cex=.7) # add variable names
  
  return(list(toret,load))   # matrix(out,ncol=nstages)
} 


# FUNCTION ForecastVitalRates  ----------------------------------
#    this script returns vital rates over time for a given pixel and climate scenario over a given time frame
# params;thispar;scen;rep;thisyears=1:89;thispixel;thisclim


ForecastVitalRates <- function(params,thispar,scen,rep,thisyears=1:89,thispixel,thisclim){ 
  
  allpars <- names(thispar)
  thisscen <- params$allscenarios[scen,]

  spateff <- terra::rast(params$zerorast)
  
  thismatfunc <- GenerateMatrixSampler(thispar = thispar, nstages=params$NSTAGES, 
                                       stagenames = params$ST_NAMES,
                                       thisdat = params$paraminfo,
                                       thisscen = thisscen,vitals=T) 
  
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$growth_clim) clim_vars_used <- c(clim_vars_used,"feb_may_tmax_","nov_apr_pr_","ann_tmin_")
  if(thisscen$phi.a_pr) clim_vars_used <- c(clim_vars_used,"nov_apr_pr_norm_")
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"spr_tmean_","mar_mar_pr_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_aug_tmean_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"july_tmean_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"mar_mar_pr_","mar_apr_mtr_")
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used,"spateff"))
  
  toret=NULL
  
  yr=1
  for(yr in thisyears){ 
    yeareff = rnorm(1)   # account for good years and bad years...
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",thisclim,yr) ) 
    thisyr_brick <- c(thisyr_brick,spateff)
    names(thisyr_brick)[terra::nlyr(thisyr_brick)] <- "spateff"
    
    climvars <- thisyr_brick[thispixel]
    #  tempdf <- thisyr_brick[haspop]
    temp <- sapply(clim_vars_unused,function(t) climvars[[t]] <<- 0   )
    
    toret2 <- list()
    toret2$year = yr
    
    toret2 <- c(toret2,thismatfunc(climvars,yeareff))
    toret <- rbind(toret,as.data.frame(toret2))
  }
  # toret
  
  return(toret)   # matrix(out,ncol=nstages)
} 


# FUNCTION "prepWorkspace" --------------------------
#    prepares the workspace with key data for DT simulation
# nreps=4; nyears=89
prepWorkspace <- function(nreps,nyears){
  # Set up workspace 
  
  params <- load_paraminfo()    # load parameter information and key R objects from file (based on google sheets info)
  
  # Prepare spatial data layers needed for sim, and store in 'params'
  params <- prepspatial(params,fn="setup_files/basehabmap.tif",rangemap="setup_files/DESERT_FOCALAREA.shp",
                        totabund=(212343/2))    # prepare spatial data
  # str(params)
  # NOTE: initial abund based on allison and mcluckie 2018...
  #    total abund is 336,393 - 124,050 =  336393 - 124050  = 212343
  #    this is only adults, so we now generate "unseen" juveniles when spreading across age classes
  
  
  # set number of reps
  params$nreps <- nreps    # should represent the total number of replicates...
  
  # Finalize Parameter list 
  
  params <- FinalizeParams(params,nyears=nyears,ddist=1000)
  
  # Make Scenarios 
  
  params$allscenarios <- makescenlist(params)   # KTS: add climate in with these scenarios
  params$nscenarios <- nrow(params$allscenarios)   # 256 scenarios (64 scens, 4 climates)
  
  
  # Final preparation of climate data 
  # Break up climate data into one-year chunks (for memory management in cluster)
  # only need to do this once- can take some time...
  
  temp <- prepclimdat(params,climyears=nyears,folder="setup_files/OutputTifs",DO=F)
  
  # Do burnin simulation (if needed) -----------------------------
  # y=5;dir="setup_files/initabund";useold=F    # KTS: now reverted to raster # KTS: now burns in separately for each climate scenario
  # originit <- rast(params$INITABUND)
  # sum(terra::global(originit,"sum",na.rm=T))   # 398355 originally
  params <- doBurnin(params,y=50,dir="setup_files/initabund",useold=T)
  # sum(terra::global(terra::rast(params$INITABUND[1]),"sum",na.rm=T))  # RCP2.6: 333622 
  # sum(terra::global(terra::rast(params$INITABUND[2]),"sum",na.rm=T))  # RCP4.5: 337117
  # sum(terra::global(terra::rast(params$INITABUND[3]),"sum",na.rm=T))  # RCP4.5: 348482
  # sum(terra::global(terra::rast(params$INITABUND[4]),"sum",na.rm=T))  # RCP8.5: 355795
  return(params)
  
}


# FUNCTION "loadPackages" -------------------------

loadPackages <- function(){
  suppressPackageStartupMessages(require(popbio))
  #library(pracma)
  suppressPackageStartupMessages(require(mvtnorm))
  suppressPackageStartupMessages(require(foreach))
  suppressPackageStartupMessages(require(doParallel))
  suppressPackageStartupMessages(require(terra))
  suppressPackageStartupMessages(require(raster))
  suppressPackageStartupMessages(require(abind))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(geodata))
  # install.packages("devtools")
  # suppressPackageStartupMessages(require(devtools))
  # install_github("aryoda/tryCatchLog")
  suppressPackageStartupMessages(require(tryCatchLog))
}


# FUNCTION "prepclimdat" -------------------------
   # chunks up climate data into one-year chunks for improved memory management
   # hopefully this will help enable simulations to be run with smaller grain size

prepclimdat <- function(params,climyears=89,folder="setup_files/OutputTifs", DO=F){
  if(DO){
    clim=1
    for(clim in 1:length(params$clim_scenarios)){
      
      this_climmod = params$clim_scenarios[clim]
      
      climvar_list <- list()    # KTS: use zeros matrix if climate effect is OFF
      v=1
      for(v in 1:length(params$clim_var_names)){
        thisvar <- params$clim_var_names[v]
        thisfile <- sprintf("%s/%s%s.tif",folder,thisvar,this_climmod)   # now read in as geotiffs
        climvar_list[[thisvar]] <- terra::rast(thisfile)
      }
      
      yr=1
      for(yr in 1:climyears){
        
        # make brick of climate data for this year
        
        thisyr_brick <- do.call("c",lapply(params$clim_var_names,function(t) climvar_list[[t]][[yr]]))
        names(thisyr_brick) <- params$clim_var_names
        towrite <- sprintf("%s/%s_%s.tif",folder,this_climmod,yr)
        terra::writeRaster(thisyr_brick,towrite,overwrite=T)
        
      }
    }
  }
  
}

# FUNCTION "ismov_func"  -----------------------
   # for determining who moves and who doesn't- used in terra::app
   # KTS: not currently used (after moving back to raster package)
# x <- rpois(10,2);prob=0.1
ismov_func <- function(x,prob){   # convert to cpp for speed?  NOTE: only works rn if one dispersal stage
  this <- x
  thisndx <- which(!is.na(x))
  if(length(thisndx)>0){
    rn <- rbinom(length(thisndx),x[thisndx],prob)
    # rn[rn==0] <- NA  
    this[thisndx] <- rn
  }
  this
}
# ismov_func(x,prob)
  

# FUNCTION "makereplist" --------------------
#  this function makes a list of replicates (unique set of draws representing parameter values for a given replicate)
# nreps=5
#  KTS: this function now stores samples as RDS files for better parallelization
makereplist <- function(params,nreps,run_name){
  
  if(!file.exists(sprintf("%s/replist_1.rds",run_name))){
    replist <- list()
    onerast <- terra::rast(params$onerast)
    
    # r=1
    for(r in 1:nreps){
      # sample a set of parameters randomly drawn from posterior distributions, etc.
      # use this same set of parameters for all climate models and all scenarios!
      replist[[r]] <- sample_params(params,baseline=F)
      saveRDS(replist[[r]],file=sprintf("%s/replist_%s.rds",run_name,r))
      
      this_spateff <- terra::app(onerast,rnorm)  # spatial variation
      fn <- sprintf("%s/spateff_%s.tif",run_name,r)
      terra::writeRaster(this_spateff,fn,overwrite=TRUE)
    }
    
  }
  
  # saveRDS(replist,file="replist.rds")
  # return(replist)
}

# FUNCTION "makescenlist" ------------------
# this function makes the list of scenarios for appending to the params object

makescenlist <- function(params){
  scenlist <- list()
  scenlist$phi.a_clim <- c(T,F)   # does adult survival depend on climate (TODO: add effect on juv surv as well?)
  scenlist$cs_clim <- c(T,F)      # does clutch size change with climate change?
  scenlist$pr_clim <- c(T,F)      # does repro prob depend directly on climate?
  scenlist$hs_clim <- c(T,F)      # does hatching success depend directly on climate
  scenlist$pf_clim <- c(T,F)      # does sex ratio depend directly on climate
  scenlist$growth_clim <- c(T,F)  # does growth and MA depend directly on climate? 
  # scenlist$growth_pr <- c(T,F)        # does growth and MA depend on precip? 
  scenlist$clim_scenario <- params$clim_scenarios
  allscenarios <- do.call(expand.grid,scenlist)    # 64 scenarios to start off... 
  
  allscenarios$clim_scenario <- as.character(allscenarios$clim_scenario)
  # nrow(allscenarios)
  
  allscenarios$index<-1:nrow(allscenarios)
  #what scenario has everything true?
  # subset(allscenarios,phi.a_pr==T & cs_clim==T & pr_clim==T & hs_clim==T & pf_clim==T  & growth_clim==T)
  #scen_1
  
  #what scenario has everything false?
  # subset(allscenarios,phi.a_pr==F & cs_clim==F & pr_clim==F & hs_clim==F & pf_clim==F  & growth_clim==F)
  #scen_64
  
  #what scenario has everything true EXCEPT adult survival ~! clim?
  # subset(allscenarios,phi.a_pr==F & cs_clim==T & pr_clim==T & hs_clim==T & pf_clim==T & growth_clim==T)
  #scen_2
  
  #what scenario has everything true EXCEPT clutchsize & PR ~ clim?
  # subset(allscenarios,phi.a_pr==T & cs_clim==F & pr_clim==F & hs_clim==T & pf_clim==T &  growth_clim==T)
  #scen_7
  
  #what scenario has everything true EXCEPT hatch success ~ clim?
  # subset(allscenarios,phi.a_pr==T & cs_clim==T & pr_clim==T & hs_clim==F & pf_clim==T &  growth_clim==T)
  #scen_9
  
  #what scenario has everything true EXCEPT sex ratio ~ clim?
  # subset(allscenarios,phi.a_pr==T & cs_clim==T & pr_clim==T & hs_clim==T & pf_clim==F & growth_clim==T)
  #scen_17
  
  #what scenario has everything true EXCEPT growth ~ clim?
  # subset(allscenarios,phi.a_pr==T & cs_clim==T & pr_clim==T & hs_clim==T & pf_clim==T & growth_clim==F)
  #scen_33
  
  
  
  NSCENARIOS <- nrow(allscenarios)
  return(allscenarios)
}


# FUNCTION "prepspatial" ------------------------
# this function reads in desert tortoise initial habitat suitability map, spreads the total 
#  abundance across all cells and writes a habitat map to file as geotiff
#  this function also ensures that key spatial data are bundled in 'params'

# fn="setup_files/basehabmap.tif";rangemap="setup_files/DESERT_FOCALAREA.shp";totabund_a=(212343/2)
prepspatial <- function(params,fn="setup_files/basehabmap.tif",rangemap="setup_files/DESERT_FOCALAREA.shp",
                        totabund_a=(212343/2)){
  
  # NOTE: terra doesn't serialize well, so raster map objects always should be read in from file
     #  current strategy: use terra for climate rasters etc, and use raster for main population structure
  
  habmap <-  terra::rast(fn) 
  habmap <-  habmap/(terra::global(habmap,"max",na.rm=T)["layer", "max"])   # Note: base.habmap map comes from the data bundle
  # should be max of 1. Min should be close to zero
  
  params$CRS <- terra::crs(habmap,proj=T)
  params$RESOLUTION <- terra::res(habmap)[1]     # in m per side
  params$CELL_AREA_KM <- (params$RESOLUTION/1000)^2   # in km2
  
  # assume anything less than 0.1 is non-habitat... (arbitrary- work with Ken to make this more realistic)
  habmap <- terra::ifel(habmap<0.1,NA,habmap)
  
  # clip habmap to rangemap (derived from the recovery units)
  range <- terra::vect(rangemap)
  range <- terra::project(range,habmap)
  # plot(habmap)
  # plot(range,add=T)
  
  habmap <- terra::mask(habmap,range)
  
  # write habmap to file and record the name
  params$HABMAP <- "setup_files/hapmap.tif"
  # unlink(params$HABMAP)
  terra::writeRaster(habmap,params$HABMAP,overwrite=TRUE)   # write map to working directory
  
  # plot(habmap)
  
  # generate initial abundance
  temp3 <- habmap
  temp3 <- terra::ifel(temp3<0.25,0,temp3)
  temp3 <- temp3 / terra::global(temp3,"sum",na.rm=T)["layer","sum"]
  # plot(temp3)
  
  initabund_a <- temp3 * totabund_a   # female adults only  [for now: just total abundance- break down into ages later!]
  # plot(datalist$initabund)
  
  ## NOTE: habitat/density relationships from Ken Nussear:
    # Encounter Rate (ER) = -0.007814 + 0.1275 * HS    
    # Density (torts per km2) = -0.1201 + 116.1 * ER  (see email from Ken 4.24.03 in Outlook)
  
  # params$initabund_ad <- "initabund_ad.tif"
  # terra::writeRaster(initabund,params$initabund_adfem,overwrite=TRUE)
  
  # derived spatial params objects
  habmap2 <- terra::ifel(is.na(habmap),0,habmap)   # same, but with zeros instead of NAs
  habmap3 <- terra::ifel(habmap2==0,0.0001,habmap2)    # needed for dispersal algorithm
  habmap3[habmap3==0] <- 0.0001
  temp <- terra::rast(habmap)
  zerorast <- terra::init(temp,0) # all zeros
  countrast <- terra::init(temp,fun="cell")  # raster value is cell index
  # params$countrast2 <- raster::setValues(raster(params$HABMAP),1:length(raster(params$HABMAP)@data@values))
  countrast <- terra::ifel(is.na(habmap),NA,countrast)  # make non-habitat cells NA to ensure no dispersal to those cells
  # plot(params$countrast)
  onerast <- terra::init(temp,1)
  
  # write to file and store filenames
  
  params$HABMAP2 <- "setup_files/habmap2.tif"
  # unlink(params$HABMAP2)
  terra::writeRaster(habmap2,params$HABMAP2,overwrite=TRUE)
  params$HABMAP3 <- "setup_files/habmap3.tif"
  # unlink(params$HABMAP3)
  terra::writeRaster(habmap3,params$HABMAP3,overwrite=TRUE)
  params$zerorast <- "setup_files/zerorast.tif"
  # unlink(params$zerorast)
  terra::writeRaster(zerorast,params$zerorast,overwrite=TRUE)
  params$countrast <- "setup_files/countrast.tif"
  # unlink(params$countrast)
  terra::writeRaster(countrast,params$countrast,overwrite=TRUE)
  params$onerast <- "setup_files/onerast.tif"
  # unlink(params$onerast)
  terra::writeRaster(onerast,params$onerast,overwrite=TRUE)
  
  ## set age/stage structure in advance of setting initial abundance
  
  params$NSTAGES <- 20     # hardcode 20 stages for now
  juvstages <- params$NSTAGES-1
     # approximate stable age distribution
  initdist <- c(.4,seq(0.15,0.02,length=juvstages-1),.7) # fraction of individuals in each age class
  # initdist[20]/sum(initdist[1:19])
  # abundances per cell, by stage 
  INITABUND2 <- initdist/sum(initdist)
  fracadult <- INITABUND2[params$NSTAGES]

  initabund <- initabund_a/fracadult   # convert to total abundance
  
  #  compute maximum/ceiling abundance for "density dependence"
  
  abundceil <-initabund * 2 # KTS: we could change this to ensure that even cells with limited hs have capacity to support a population
  params$ABUNDCEIL <- "setup_files/abundceil.tif"
  # unlink(params$ABUNDCEIL)
  terra::writeRaster(abundceil,params$ABUNDCEIL,overwrite=TRUE)

  #  set initial abundance
  
  tempinit <- lapply(1:(params$NSTAGES),function(t) initabund * INITABUND2[t] ) 
  tempinit <- do.call("c",tempinit)
  params$ST_NAMES <- c(paste0("F_",1:params$NSTAGES))
  names(tempinit) <-  params$ST_NAMES
  tempinit <- round(tempinit)
  # plot(tempinit[[20]])
  params$INITABUND <- "setup_files/initabund.tif"   # starting abundance by stage
  # unlink(params$INITABUND)
  terra::writeRaster(tempinit,params$INITABUND,overwrite=TRUE)
  
  # plot(params$INITABUND$F_1)
  
  # terra::global(terra::app(tempinit,sum),"sum",na.rm=T)

  return(params)
}

# FUNCTION "DoBurnin --------------------------

# run simulation across the landscape under baseline conditions
   #  to ensure that initial abundance can be used as a reference.
   #  returns an updated stable initial abundance raster
# y=10;file=paste0("new_initabund",Sys.Date(),".tif"); useold=T
#  KTS: this now uses the first year of climate data to ensure that differences in habitat suitability
    # across the range is accounted for
# KTS: this now uses a random sample from the first 10 years of climate data to ensure that the burnin
#    accounts for both spatial and temporal variation

doBurnin <- function(params,y,dir="setup_files/initabund", useold=T){
  
  if(!dir.exists(dir)) dir.create(dir)
  ceil_filename <- sprintf("%s/ceiling.tif",dir)
  params$ABUNDCEIL <- ceil_filename
  
  if(useold){
    filelist <- sprintf("%s/initabund_%s.tif",dir,params$clim_scenarios)
  }else{
    # unlink(file)
    filelist = NULL
    thisparams <- sample_params(params,baseline=T)
    thisscen <- params$allscenarios[1,]
    
    # generate a function for generating the matrix
    # thispar = thisparams; nstages=params$NSTAGES; stagenames = params$ST_NAMES;thisdat = params$paraminfo; vitals=F
    thismatfunc <- GenerateMatrixSampler(thispar = thisparams, nstages=params$NSTAGES, 
                                         stagenames = params$ST_NAMES,
                                         thisdat = params$paraminfo,
                                         thisscen = thisscen,vitals=F) 
    
    yrsamps <- sample(1:10,y,replace = T)
    clim=1
    for(clim in 1:length(params$clim_scenarios)){
      this_climmod <- params$clim_scenarios[clim]
      # thisabund <- terra::rast(params$INITABUND)
      thisabund <- raster::brick(params$INITABUND)
            # for debug
      # raster::writeRaster(thisabund,sprintf("sims/clim_%s_burnin_%s.tif",clim,0),overwrite=TRUE) 
      
      # initabund_tot <- terra::app(thisabund,sum)
      initabund_tot <- sum(thisabund)
      
      spateff <- terra::rast(params$zerorast)
      
      yr=3
      for(yr in 2:y){   
        
        yeareff = 0 #rnorm(1)   # account for good years and bad years...
        
        # read brick of climate data for this year
        thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",this_climmod,yrsamps[yr]) ) 
        
        thisyr_brick <- c(thisyr_brick,spateff)
        names(thisyr_brick)[terra::nlyr(thisyr_brick)] <- "spateff"
        
        # totabund <- terra::app(thisabund,sum)  # compute total abundance
        totabund <- sum(thisabund)
        over <- totabund > initabund_tot*1.1   # exceeds limit
        # over_cells <- which(terra::values(over,dataframe=T)[,1])
        over_cells <- which(raster::getValues(over))
        
        haspop <- totabund>0
        # haspop_cells <- which(terra::values(haspop,dataframe=T)[,1])
        haspop_cells <- which(raster::getValues(haspop))
        
        if(length(over_cells)>0){
          which_over <- sapply(over_cells,function(t) which(haspop_cells==t))
        }else{
          which_over <- numeric(0)
        } 
        
        temp <- rast("setup_files/OutputTifs/rcp26_1.tif")
        baseclim <- temp[1]
        temp <- lapply(params$clim_var_names,function(t) baseclim[[t]] <<- 0)
        baseclim$spateff <- 0
        
        badmat <- thismatfunc(baseclim,yeareff=-0.5)*0.75
        
        if(length(haspop_cells)>0){
          this_climext <- thisyr_brick[haspop_cells]
          mat_list <- lapply(1:nrow(this_climext),function(t) thismatfunc(climvars=this_climext[t,],yeareff=yeareff) )
          if(length(which_over)>0) mat_list[which_over] <- lapply(1:length(which_over),function(x) badmat )
        }else{
          mat_list <- NULL
        }
        
        prevabund <- thisabund
        
        # do matrix updating with stochasticity
        #     only apply these functions to the grid cells they correspond to 
        
        tempdf <- prevabund[haspop_cells]
        tofill <- t(sapply(1:nrow(tempdf), function(t) do.stoch.proj(as.numeric(tempdf[t,]),mat_list[[t]],matndx=params$MATRIX_NDX)    ))
        thisabund[haspop_cells] <- tofill
        
        # store data (for debug)
        # raster::writeRaster(thisabund,sprintf("sims/clim_%s_burnin_%s.tif",clim,yr),overwrite=TRUE)
      }  # end year loop
      
      thisinit <- thisabund
      thisfile <- sprintf("%s/initabund_%s.tif",dir,this_climmod)
      filelist <- c(filelist,thisfile)
      terra::writeRaster(thisinit,thisfile,overwrite=TRUE)  # save it...
      
      if(this_climmod=="rcp26"){   # KTS: set maximum abundance based on rcp2.6
        temp2 <- sum(thisinit)
        thisabundceil <- temp2 * 2  # maxabund is 2*initabund...
        # oldabundceil <- terra::rast(params$ABUNDCEIL)
        # terra::global(thisabundceil,"sum",na.rm=T)
        # terra::global(oldabundceil,"sum",na.rm=T)
        # unlink(params$ABUNDCEIL)
        raster::writeRaster(thisabundceil,params$ABUNDCEIL,overwrite=TRUE)  # overwrite previous abundance ceiling
      }
      
    } # end clim loop
  
  } # end 'useold==F' 
    
  params$INITABUND <- filelist   # now has separate initabund for every clim scenario
  # temp2 <- terra::app(thisinit,sum)
  
  return(params)
}


# FUNCTION "MakeEmpirGenerator" ---------------

MakeEmpirGenerator <- function(vec){
  function(rand=runif(1)){
    ndx <- floor(rand*length(vec)+1)
    vec[ndx]
  }
}
# MakeEmpirGenerator(thismcmc[[thisvar]])()

# FUNCTION "MakeUnifGenerator" -----------------

MakeUnifGenerator <- function(min, max){
  function(rand=runif(1)){
    rand*(max-min)+min
  }
}
# MakeUnifGenerator(2,4)()

# a <- MakeUnifGenerator(2,4)

# FUNCTION "MakeNormalGenerator" -----------------

MakeNormalGenerator <- function(mean, sd, constraint="none"){
  if(constraint=="none"){
    ret = function(randn=rnorm(1)){
      (randn*sd)+mean
    }
  }else if(constraint=="c01"){
    ret = function(randn=rnorm(1)){
      max(0.001,min(.999,(randn*sd)+mean))
    }
    
  }
  
  return(ret)
}
# MakeNormalGenerator(0.5,0.5,"c01")()

# FUNCTION "load_paraminfo"  -------------------------
# Load key parameter information and R objects listed in google sheets

load_paraminfo <- function(){
  
  datalist <- list()
  
  # read in CSV file with base parameter values
  thiscsv <- read.csv("setup_files/InitialDistributions_DT.csv")
  thiscsv <- thiscsv[,1:(ncol(thiscsv)-1)]
  thiscsv <- thiscsv[thiscsv$Used=="Y",]
  # head(thiscsv)
  
  thisrobj <- read.csv("setup_files/RObjects_DT.csv")
  thisrobj <- thisrobj[,1:(ncol(thisrobj)-1)]
  thisrobj <- thisrobj[thisrobj$Used=="Y",]
  # thisrobj
  
  load("setup_files/RObjects.RData")    # load R objects
  
  datalist$allparams <- sort(unique(c(thiscsv$Parameter,thisrobj$Parameter)))
  temp <- sapply(1:length(datalist$allparams),function(t)  datalist$paraminfo[[datalist$allparams[t]]] <<- list() )
  
  r=20
  ## loop through params.data
  temp <- sapply(1:nrow(thiscsv),function(r){
    thisparam <- thiscsv$Parameter[r]
    if(thiscsv$Type[r]=="ancillary"){
      datalist$paraminfo[[thisparam]][[thiscsv$Subparameter[r]]] <<- thiscsv$Value[r]
    }else if(thiscsv$Type[r]=="param"){
      datalist$paraminfo[[thisparam]][[thiscsv$Subparameter[r]]]$constraint <<- thiscsv$Constraint[r]
      datalist$paraminfo[[thisparam]]$subparams <<- c(datalist$paraminfo[[thisparam]]$subparams,thiscsv$Subparameter[r])
      if(thiscsv$Distribution[r]=="uniform"){
        datalist$paraminfo[[thisparam]][[thiscsv$Subparameter[r]]]$gen <<- MakeUnifGenerator(thiscsv$Min[r],thiscsv$Max[r])
      } else if(thiscsv$Distribution[r]=="normal"){
        datalist$paraminfo[[thisparam]][[thiscsv$Subparameter[r]]]$gen <<- MakeNormalGenerator(thiscsv$Mean[r],thiscsv$SD[r],thiscsv$Constraint[r])
      }
    }
  }   )
  
  ## loop through model objects
  temp <- sapply(1:nrow(thisrobj),function(r){
    thisparam <- thisrobj$Parameter[r]
    datalist$paraminfo[[thisparam]]$subparams <<- c(datalist$paraminfo[[thisparam]]$subparams,thisrobj$Subparameter[r])
    if(thisrobj$Type[r]=="Posterior"){
      thispost <- thisrobj$Model_out[r]
      thismcmc <- eval(parse(text=thispost))   
      thisvar <- thisrobj$R_obj[r]
      datalist$paraminfo[[thisparam]][[thisrobj$Subparameter[r]]]$post <<- MakeEmpirGenerator(thismcmc[[thisvar]])
    }else if(thisrobj$Type[r]=="Function"){
      thisfunc <- thisrobj$R_obj[r]
      datalist$paraminfo[[thisparam]][[thisrobj$Subparameter[r]]]$func <<- eval(parse(text=thisfunc))
    }else if(thisrobj$Type[r]=="LM_object"){
      datalist$paraminfo[[thisparam]][[thisrobj$Subparameter[r]]]$mpar <<- thisrobj$R_obj[r]
      thisobj <- thisrobj$Model_out[r]
      datalist$paraminfo[[thisparam]][[thisrobj$Subparameter[r]]]$mod <<- eval(parse(text=thisobj))
    } else if(thisrobj$Type[r]=="R_object"){
      thisobj <- thisrobj$R_obj[r]
      datalist$paraminfo[[thisparam]][[thisrobj$Subparameter[r]]]$obj <<- eval(parse(text=thisobj))
    }
  }   )

  
    ## TEST CODE (for debugging)
  # datalist$paraminfo$gamma$dispersal()
  # get("max", env = environment(datalist$paraminfo$gamma$dispersal))
  # datalist$paraminfo$CS$junetmax.mn
  # datalist$paraminfo$MA$ma.intercept.mn()
  # datalist$paraminfo$HS$hs.tempeff.mn()
  # get("mean", env = environment(datalist$paraminfo$HS$hs.tempeff.mn))
  # datalist$paraminfo$NS$nestsucc()
  # get("max", env = environment(datalist$paraminfo$NS$nestsucc))
  # datalist$paraminfo$CS$CS_intercept() 
  # datalist$paraminfo$CS$CS_bodysize()
  # 
  # read in RData file with models for expressing params as a function of climate. 

 
  return(datalist)
}


# FUNCTION "sample_params" ---------------------------

sample_params <- function(params,baseline=F){
  if(!baseline){
    post_draw <- runif(1)
    toret <- list()
    mainparams <- params$allparams
    i=4
    for(i in 1:length(mainparams)){
      thisparam <- mainparams[i]
      toret[[thisparam]] <- list()
      thisinfo <- params$paraminfo[[thisparam]]
      thissubs <- thisinfo$subparams
      j=1
      if(length(thissubs)>0){
        for(j in 1:length(thissubs)){
          thissub <- thissubs[j]
          thisinfo2 <- thisinfo[[thissub]]
          
          if(!is.null(thisinfo2$post)){
            toret[[thisparam]][[thissub]]$this <- thisinfo2$post(post_draw)
          }else if(!is.null(thisinfo2$gen)){
            toret[[thisparam]][[thissub]]$this <- thisinfo2$gen()
          }else if(!is.null(thisinfo2$mod)){
            thismod <- thisinfo2$mod      # need two passes to deal with model objects
            thisvar <- thisinfo2$mpar
            varcovar <- vcov(thismod)
            varnames <- colnames(varcovar)
            # varcovar <- varcovar
            toret[[thisparam]]$mod <- thismod
            toret[[thisparam]]$modpars <- c(toret[[thisparam]]$modpars,thisvar)
            toret[[thisparam]]$subpars <- c(toret[[thisparam]]$subpars,thissub)
            # thissubs     # KTS: hardcoded for now
            
          }else if(!is.null(thisinfo2$func)){
            toret[[thisparam]][[thissub]]$func <- thisinfo2$func
          }
          
        }
      }
      
      i=4
      for(i in 1:length(mainparams)){
        thisparam <- mainparams[i]
        thisinfo <- params$paraminfo[[thisparam]]
        thissubs <- thisinfo$subparams
        if(!is.null(toret[[thisparam]]$modpars)){
          thismod <- toret[[thisparam]]$mod
          thismodpars <- toret[[thisparam]]$modpars
          thissubpars <- toret[[thisparam]]$subpars
          thisvcov <- vcov(thismod)
          thisvcov <- thisvcov[thismodpars,thismodpars,drop=F]
          thismeans <- coef(thismod)[thismodpars]
          thissamp <- rmvnorm(1,thismeans,thisvcov)
          v=1
          for(v in 1:length(thissubpars)){
            thissub <- thissubpars[v]
            toret[[thisparam]][[thissub]]$this <- thissamp[1,thismodpars[v]]
          }
        }
      }
    }
    
  }else{     # for baseline/mean draw (for running burn-in)
    post_draw <- runif(1)
    toret <- list()
    mainparams <- params$allparams
    i=4
    for(i in 1:length(mainparams)){
      thisparam <- mainparams[i]
      toret[[thisparam]] <- list()
      thisinfo <- params$paraminfo[[thisparam]]
      thissubs <- thisinfo$subparams
      j=4
      if(length(thissubs)>0){
        for(j in 1:length(thissubs)){
          thissub <- thissubs[j]
          thisinfo2 <- thisinfo[[thissub]]
          
          if(!is.null(thisinfo2$post)){
            toret[[thisparam]][[thissub]]$this <- mean(replicate(100,thisinfo2$post()))
          }else if(!is.null(thisinfo2$gen)){
            toret[[thisparam]][[thissub]]$this <- mean(replicate(100,thisinfo2$gen()))
          }else if(!is.null(thisinfo2$mod)){
            thismod <- thisinfo2$mod      
            thisvar <- thisinfo2$mpar
            varcovar <- vcov(thismod)
            varnames <- colnames(varcovar)
            toret[[thisparam]]$mod <- thismod
            toret[[thisparam]]$modpars <- c(toret[[thisparam]]$modpars,thisvar)
            toret[[thisparam]]$subpars <- c(toret[[thisparam]]$subpars,thissub)
            thismeans <- coef(thismod)[toret[[thisparam]]$modpars]
            v=1
            for(v in 1:length(toret[[thisparam]]$subpars)){
              thissub <- toret[[thisparam]]$subpars[v]
              toret[[thisparam]][[thissub]]$this <- thismeans[toret[[thisparam]]$modpars[v]]
            }
            
            # thissubs     # KTS: hardcoded for now
            
          }else if(!is.null(thisinfo2$func)){
            toret[[thisparam]][[thissub]]$func <- thisinfo2$func
          }
          
        }
      }

    }
    
  }
  
  return(toret)
}

# sample_params(params)


## FUNCTION "MakeLHSSamples"  (for global sensitivity testing) -----------------------

# Samples from the uniform LHS and translates into desired parameter space
#  Returns a master data frame that will specify params for all runs
#  KTS: not used at the moment

MakeLHSSamples <- function(add=F){
  
  LHSParms <- list()    # initialize the container for parameter bounds
  
  ####  FECUNDITY
  LHSParms <- specifyLHSParam(LHSParms,"FECUNDITY",type="CONT",lb=0,ub=1) #ub=2?
  
  #### SURVIVAL
  LHSParms <- specifyLHSParam(LHSParms,"SURVIVAL",type="CONT",lb=0,ub=1)
  
  ### ROADMORT
  LHSParms <- specifyLHSParam(LHSParms,"ROADMORT",type="CONT",lb=0.01,ub=0.2)       
  
  ### TRAPMORT
  LHSParms <- specifyLHSParam(LHSParms,"TRAPMORT",type="CONT",lb=0.01,ub=0.25)
  
  ### DISPF
  LHSParms <- specifyLHSParam(LHSParms,"DISPF",type="CONT",lb=0.1,ub=0.4) 
  
  ### DISPM
  LHSParms <- specifyLHSParam(LHSParms,"DISPM",type="CONT",lb=0.4,ub=0.8) 
  
  #### FRACINIT
  LHSParms <- specifyLHSParam(LHSParms,"FRACINIT",type="CONT",lb=0.25,ub=0.9)
  
  #### MAXABUND
  LHSParms <- specifyLHSParam(LHSParms,"MAXABUND",type="INT",lb=20,ub=40)
  
  # GENERATE LATIN HYPERCUBE SAMPLE
  
  nVars <- length(names(LHSParms))  
  
  LHS <- lhs::randomLHS(N_LHS_SAMPLES, nVars)   # generate multiple samples from parameter space according to a LHS sampling scheme
  
  temp <- as.data.frame(LHS)
  
  ### translate raw lhs samples into desired parameter space
  colnames(temp) <- names(LHSParms)
  parm=1
  for(parm in 1:nVars){
    if(LHSParms[[parm]]$type=="CONT"){
      temp[,parm] <- LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb)
    }
    if(LHSParms[[parm]]$type=="CAT"){
      temp[,parm] <- ceiling(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
    }
    if(LHSParms[[parm]]$type=="INT"){
      temp[,parm] <- round(LHSParms[[parm]]$lb + LHS[,parm]*(LHSParms[[parm]]$ub-LHSParms[[parm]]$lb))
    }
  }
  
  if(add==FALSE) masterDF <- temp    #  storage container (data frame) to record relevant details for each MP file. Rows:MP file/LHS samples. Cols: relevant variables
  
  if(add==TRUE) masterDF <- rbind(masterDF,temp)
  
  ## name file for LHS parameters 
  write.csv(masterDF,sprintf("masterDF_prelim%s.csv",Sys.Date()),row.names=F)
  
  return(masterDF)
}

## FUNCTION "specifyLHSParam"   (for global sensitivity testing) -----------------------

# Information necessary to translate standard uniform LHS sample into parameters of interest for paleo project 

specifyLHSParam <- function(paramslist,name,type,lb,ub){
  newlist <- paramslist
  eval(parse(text=sprintf("newlist$%s <- list()",name)))
  eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
  eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
  eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
  return(newlist)
}

# Geometric mean function --------------------

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}


# FUNCTION 'do.stoch.proj': project one year with demographic stochasticity  --------------------

# prev=as.numeric(tempdf[2,])   
# mymat = mat_list[[2]]
# matndx = params$MATRIX_NDX
do.stoch.proj <- function(prev,mymat,matndx){
  outabund <- floor(prev) 				
  nstages <- ncol(mymat)
 
  # if(!any(is.na(outabund))){
    # if(any(outabund > 0)){
  outabund[1] <- rpois(1,sum(mymat[1,]*prev[1:nstages])) #sum(rpois(nstages,mymat[1,]*prev[1:nstages]))   
  outabund[2:nstages] <- sapply(2:nstages,function(z) sum(rbinom(length(matndx[[z]]),prev[matndx[[z]]],mymat[z,matndx[[z]]]))   )      #sapply(2:nstages,function(t) sum(rbinom(nstages,prev[1:nstages],mymat[t,])) )

    # }
  # }
  
  return(outabund)
}
# do.stoch.proj(prev,mymat)

# do.stoch.proj(prev,mymat)

# FUNCTION 'GenerateMatrixSampler': Generate function for sampling matrices for a given replicate -----------
#    KTS: in some ways, this is the 'critical habitat breadth' model, and the main simulation model
#        is the 'metapopulation model'

# thispar = thisparams
# nstages=params$NSTAGES
# stagenames = params$ST_NAMES
# thisdat = params$paraminfo
# vitals=F
GenerateMatrixSampler <- function(thispar,nstages,stagenames,thisdat,thisscen,vitals=F){   
  allpars <- names(thispar)
    ## NOTE: to be more efficient we may want to keep track of which params are climate sensitive
  basemat <- matrix(0,nrow=nstages,ncol=nstages)
  rownames(basemat) = stagenames
  colnames(basemat) = stagenames
  
  mat_ndx <- matrix(1:nstages^2,ncol=nstages)
  ndx_subdiags <- diag(mat_ndx[-1,-ncol(mat_ndx)])

  if(vitals) toret=list()

  # generate VB function for this replicate... (to model tortoise growth)
  
  thiskparams = c(
    k.site.prec = thispar$GROW$ANC_growth.k.siteprec$this,
    k.all.meanl = thispar$GROW$ANC_growth.k.int$this,
    k.t.maxspr = thispar$GROW$ANC_growth.k.tmax$this,
    k.pptwinsp = thispar$GROW$ANC_growth.k.precip$this,
    k.tmin = thispar$GROW$ANC_growth.k.tmin$this
  )
  thisaparams = c(
    a.site.prec = thispar$GROW$ANC_growth.a.siteprec$this,
    a.all.meanl = thispar$GROW$ANC_growth.a.int$this,
    a.tmax = thispar$GROW$ANC_growth.a.tmax$this,
    a.ppt = thispar$GROW$ANC_growth.a.precip$this,
    a.tmin = thispar$GROW$ANC_growth.a.tmin$this
  )
  thist0 = thispar$GROW$ANC_growth.t0$this
  
  thisSAM = thispar$MA$size.at.maturity$this
  
  thisPHI_J = max(0.5,min(0.99,thispar$PHI_J$juvsurv$this))    # KTS: should make this climate dependent if adult survival is??
  
  TEST=F
  # TEST=T
  if(TEST){
    temp <- rast(sprintf("setup_files/OutputTifs/%s_1.tif",params$clim_scenarios[1]))
    climvars <- temp[1]
    temp <- lapply(params$clim_var_names,function(t) climvars[[t]] <<- 0)
    climvars$spateff <- 0
    yeareff=0
  }
  
  matrix_sampler <- function(climvars,yeareff=0){
    
    thisclimvars = c(
      t.maxspr = as.numeric(climvars[["feb_may_tmax_"]]/thisdat$GROW$t.maxspr.sd)*thisscen$growth_clim,   # KTS: changed from spr_tmean_
      pptwinsp = as.numeric(climvars[["nov_apr_pr_"]]/thisdat$GROW$pptwinsp.sd)*thisscen$growth_clim,
      tmin = as.numeric(climvars[["ann_tmin_"]]/thisdat$GROW$tmin.sd)*thisscen$growth_clim
    )
    
    #  Use this VB generator in the other functions!  # KTS: move this into matrix sampler to enable site effect to vary by pixel?
    thisVB = thispar$GROW$ANC_VBgen$func(a.params=thisaparams,k.params=thiskparams,climvars=thisclimvars,t0=thist0)
    # thisVB(30)
    
    thisMA = max(9,min(nstages,thispar$GROW$ANC_MAfunc$func(thisVB,thisSAM)))
    if(vitals) toret$MA=thisMA
    
    this_adult.MCL = thispar$GROW$ANC_adultsize$func(thisVB,thisMA)
    if(vitals) toret$adMCL <- this_adult.MCL
    
    # spateff <- rnorm(1)     # TODO double check the scale of the precip effect
    varpart <- sqrt(thispar$PHI_A$PHI_A_siteyear.sd$this^2/2)
    thisbasesurv <- plogis(qlogis(thispar$PHI_A$PHI_A_intercept$this) + varpart*yeareff + varpart*climvars[["spateff"]])   # using PSP model for base survival
    thisPHI_A = max(0.75,min(0.99, plogis(qlogis(thisbasesurv) + thispar$PHI_A$PHI_A_precip$this * (climvars[["sepx2_mar_pr_"]]/thisdat$PHI_A$precip.sd) * thisscen$phi.a_clim +   
                                                thispar$PHI_A$PHI_A_tmax$this * (climvars[["mar_nov_tmax_"]]/thisdat$PHI_A$tmax.sd) * thisscen$phi.a_clim ) ))    # KTS: updated with Steve's model- update param table with sd, make sure to implement standardization here
    if(vitals) toret$PHI_A = thisPHI_A
    
    thisCS = max(0.01,min(12,exp(thispar$CS$CS_intercept$this + thispar$CS$CS_springtmean$this * (climvars[["spr_tmean_"]]/thisdat$CS$mmt.sd) * thisscen$cs_clim +   # KTS: update with my newest version of mitchell model- just need sd
                                   thispar$CS$CS_bodysize$this * ((this_adult.MCL-thisdat$CS$bodysize.mn)/(thisdat$CS$bodysize.sd)) +
                                   thispar$CS$CS_precip$this * (climvars[["mar_mar_pr_"]]/thisdat$CS$precip.sd) * thisscen$cs_clim) ))   # KTS check
    if(vitals) toret$CS = thisCS
    
    thisHS = max(0.001,min(.99, as.numeric(thispar$HS$HS_clim_eff$func(deltaT = climvars[["jun_jul_aug_tmean_"]] * thisscen$hs_clim,   # EH check? [we're good!?]
                                                                       hs.tempeff.draw = thispar$HS$HS_nest_tempeff$this,
                                                                       nest.hs.slope.draw = thispar$HS$nest.hs.slope$this,
                                                                       hs_mean = thispar$HS$hs.mn$this  ) ) ))
    if(vitals) toret$HS = thisHS
    
    thisNS = max(0.001,min(.99, thispar$NS$nestsucc$this))
    
    thisPFmod = thispar$PF$PF_clim_eff$func(deltaTjuly = climvars[["july_tmean_"]] *thisscen$pf_clim, 
                                            t.mean.ref = thisdat$PF$t.mean.ref,
                                            hs.tempeff.draw = thispar$HS$HS_nest_tempeff$this,    # NOTE: using HS draw so that they are the same -look into if there's more elegant way
                                            sex.tmean.draw = thispar$PF$sex.tmean$this)
    thisPF = max(0.001,min(.999, thisPFmod["perc.female"] ))
    if(vitals) toret$PF = thisPF
    
    thisPHI_H = max(0.001,min(.999, thispar$PHI_H$hatchsurv$this))
    
    # TODO: implement the prior reproduction effect
    thisPR = max(0.001,min(.999, plogis(thispar$PR$PR_intercept$this + thispar$PR$PR_priorrep$this*0.9 +
                                          thispar$PR$PR_precip$this * (climvars[["mar_mar_pr_"]]/thisdat$PR$precip.sd) * thisscen$pr_clim  +   
                                          thispar$PR$PR_bodysize$this * ((this_adult.MCL-thisdat$PR$bodysize.mn)/thisdat$PR$bodysize.sd) +
                                          thispar$PR$PR_temprange$this * (climvars[["mar_apr_mtr_"]]/thisdat$PR$mtr.sd) * thisscen$pr_clim  ) ))   
    
    if(vitals) toret$PR <- thisPR
  
    
       # construct matrix
    thismat = basemat
    thismat[1,thisMA:nstages] = thisPR*thisCS*thisNS*thisHS*thisPF*thisPHI_H
    thisMA_J = thisMA
    thismat[ndx_subdiags[1:(thisMA_J-1)]] = thisPHI_J
    thismat[ndx_subdiags[thisMA_J:nstages]] = thisPHI_A
    thismat[nstages,nstages] = thisPHI_A
    
    if(vitals){
      toret$lambda = popbio::lambda(thismat)
      return(toret)
    }else{
      return(thismat)
    }  
  }
  
  return(matrix_sampler)   # matrix(out,ncol=nstages)
}

# FUNCTION 'make_dispmat': define dispersal kernel for the lattice population model  ---------

# dmean <- params$DISP_DIST
# dmax <- params$DISP_MAX
# res <- raster::res(params$HABMAP)[1]
make_dispmat <- function(dmean,dmax,res){
  levs <- floor(dmax/res)
  levs2 <- 2*levs + 1
  out <- matrix(0,nrow=levs2,ncol=levs2,byrow = T)
  dist <- out
  cent <- levs2-levs
  i=2
  for(i in 1:levs){
    dist[cent+i,cent] <- res*i
    dist[cent-i,cent] <- res*i
    dist[cent,cent+i] <- res*i
    dist[cent,cent-i] <- res*i
    if(i==levs){
      for(j in 1:min(i,levs)){
        k=1
        for(k in 1:min(i,levs)){
          dist[cent+j,cent+k] <- sqrt((res*j)^2+(res*k)^2)
          dist[cent+j,cent-k] <- sqrt((res*j)^2+(res*k)^2)
          dist[cent-j,cent+k] <- sqrt((res*j)^2+(res*k)^2)
          dist[cent-j,cent-k] <- sqrt((res*j)^2+(res*k)^2)
        }
      }
    }
  }
  
  lam <- 1/dmean
  out <- exp(-lam*dist)
  out[dist>dmax] <- 0
  out[cent,cent] <- 0
  out <- out/sum(out)
  return(out)
}


# FUNCTION 'dispfn': do dispersal - for use with 'raster::focal' moving window function ----------
#   For debug...
# x=fv[2,];move=getValues(move_f);emptiness=getValues(howempty);ceil=getValues(raster(params$ABUNDCEIL));dispmat=params$dispmat;center=params$cent;habsuit=getValues(raster(params$HABMAP3));zero=getValues(raster(params$zerorast))
dispfn <- function(move=move_f,emptiness,ceil,dispmat=params$dispmat,center=params$cent,habsuit=params$HABMAP3,zero=params$zerorast){
  thisdisp <- zero
  thisempt <- emptiness
  function(x){   # x is window of dispersal indices
    thiscell <- x[center]
    thisval <- move[thiscell]   #[1,1]
    x[center] <- NA
    thisndx <- which(!is.na(x))    # relative available locations within dispersal window
    nonna <- x[thisndx]   # indices of non-na elements within the dispersal regions (focal adds NAs as 'padding' to account for edge issues)
    if(length(nonna)>0){
      hs <- habsuit[nonna]  #[,1]    # absolute suitability of those non-na elements
      ceilabund <- ceil[nonna]  #[,1]
      empwt <- thisempt[nonna]/ceilabund    # [,1]  # degree to which there is 'free space' in these non-na elements
      empwt[empwt<0] <- 0
      empwt[is.nan(empwt)] <- 0
      empwt <- empwt + 0.1    # make the effect of emptiness a little less strong
      # if(sum(dispmat[thisndx]*hs*empwt)>0){
      todisp <- rmultinom(1,thisval,dispmat[thisndx]*hs*empwt)[,1]
      thisdisp[nonna] <<- thisdisp[nonna] + todisp
      thisempt[nonna] <<- thisempt[nonna] - todisp    # update the emptiness raster each time
      # }else{
      # thisdisp[thiscell] <<- thisval
      # }
    }else{
      thisdisp[thiscell] <<- thisval
    }
    
    return(1)
  }
}


# FUNCTION 'FinalizeParams': set final params needed for simulation ---------------

# nyears=2;ddist=1000
FinalizeParams <- function(params,nyears=2,ddist=1000){
  # params <- list()
  params$NYEARS <- nyears
  
  params$MATRIX_NDX <- list()    # index of which matrix elements to consider for each transition (for efficient matrix computation)
  params$MATRIX_NDX[[1]] <- NA
  params$MATRIX_NDX[2:(params$NSTAGES-1)] <- lapply(1:(params$NSTAGES-2), function(t) t  )
  params$MATRIX_NDX[[params$NSTAGES]] <- c(params$NSTAGES-1,params$NSTAGES)
  
  # maxabund <- params$max_density*params$CELL_AREA_KM  # not needed for now?
  
  # params$ABUNDCEIL_MAX <- maxabund     # abundance ceiling per cell in high-quality habitat (not used for now) 
  params$DISP_STAGE <- 20   # note: all adults should disperse- adult stages depend on temp etc... 
  
  # params$DISP_PROB <- dprob   # female dispersal prob  (sampled in params- not needed)    
  params$DISP_DIST <- ddist      # in meters
  params$DISP_MAX <- 8000        # set to only disperse to neighboring cells
  
  ## specify climate scenarios and variable names (move into params?)
  
  params$clim_scenarios = c("rcp26","rcp45","rcp60","rcp85")
  
  params$clim_var_names <- c( 
    "mar_mar_pr_",
    "mar_apr_mtr_",
    "spr_tmean_",
    "nov_apr_pr_" ,
    "ann_tmin_",
    "feb_may_tmax_",  
    # "nov_apr_pr_norm_",
    "july_tmean_",
    "jun_jul_aug_tmean_",
    "sepx2_mar_pr_",
    "mar_nov_tmax_"
    #"spateff"
  )

  # GET READY TO DISPERSE
 
  # hs_raster@crs
  # raster::res(hs_raster)
  
  params$dispmat <- make_dispmat(dmean=params$DISP_DIST,dmax=params$DISP_MAX,res=params$RESOLUTION)
  params$weightsmat <- matrix(1,nrow(params$dispmat),ncol(params$dispmat),byrow = T)
  params$cent <- ceiling(length(params$dispmat)/2)   # identify center cell
  
  countrast <- terra::rast(params$countrast)
  params$fv <- terra::focalValues(x=countrast,w=nrow(params$weightsmat),fill=NA)
  params$fv[is.nan(params$fv)] <- NA
  
  return(params)
} 


# FUNCTION 'doReplicate': run a single complete replicate ------------------
   # KTS: now runs all scenarios as part of each replicate
   # KTS: now reverted to raster (just the main abundance structure, all climate data are in terra)
doReplicate <- function(params,run_name,rep,scen,dispersal=T){
  
  # first, read in the set of parameters randomly drawn from posterior distributions, etc.
  # use this same set of parameters for all climate models and all scenarios!
  
  thisparams <- readRDS(sprintf("%s/replist_%s.rds",run_name,rep))   #params$allreps[[rep]]
  
  spateff <- terra::rast(sprintf("%s/spateff_%s.tif",run_name,rep))  # spatial variation
  # plot(spateff)

  abundceil <- raster::raster(params$ABUNDCEIL)
  
  countrast <- raster::raster(params$countrast)  # maybe not needed?
  
  zerorast <- raster::raster(params$zerorast)
  
  habmap3 <- raster::raster(params$HABMAP3)
  
  thisscen <- params$allscenarios[scen,]   # this scenario specs
  
  # generate a function for generating the matrices
  # thispar = thisparams; nstages=params$NSTAGES; stagenames = params$ST_NAMES;thisdat = params$paraminfo;vitals=F
  thismatfunc <- GenerateMatrixSampler(thispar = thisparams, nstages=params$NSTAGES, 
                                       stagenames = params$ST_NAMES,
                                       thisdat = params$paraminfo,
                                       thisscen = thisscen,
                                       vitals=F) 
  
  this_climmod <- thisscen$clim_scenario 
  clim <- which(params$clim_scenarios==this_climmod)  

    # thisabund <- terra::rast(params$INITABUND)
  thisabund <- raster::brick(params$INITABUND[clim])  # KTS- now initialized differently by climate scen
  names(thisabund) <- params$ST_NAMES
  raster::writeRaster(thisabund,paste0(simspath,"/scenario_",scen,"_rep_",rep,"_year_",0,"_",this_climmod,".tif"),overwrite=TRUE)
  
	#make fresh poparray array if loop is starting
  poparray<-array(0,dim=c(nyears+1,params$NSTAGES))
  abunds <- cellStats(thisabund,sum)
  #abunds <- raster::global(thisabund,"sum",na.rm=T)[,1]
  poparray[1,] <- abunds #store abundances
  
    
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$growth_clim) clim_vars_used <- c(clim_vars_used,"feb_may_tmax_","nov_apr_pr_","ann_tmin_")
  if(thisscen$phi.a_clim) clim_vars_used <- c(clim_vars_used,"sepx2_mar_pr_","mar_nov_tmax_")
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"spr_tmean_","mar_mar_pr_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_aug_tmean_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"july_tmean_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"mar_mar_pr_","mar_apr_mtr_")
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used,"spateff"))
  
  yr=2
  for(yr in 2:(params$NYEARS+1)){
    
    yeareff = rnorm(1)   # account for good years and bad years...
    
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",this_climmod,yr-1) ) 
    thisyr_brick <- c(thisyr_brick,spateff)
    names(thisyr_brick)[terra::nlyr(thisyr_brick)] <- "spateff"
    
    # totabund <- terra::app(thisabund,sum)  # compute total abundance
    totabund <- sum(thisabund)
    
    # do dispersal
    if(dispersal){
      dispstag_f <- paste0("F_",params$DISP_STAGE)
      disp_prob_f <- thisparams$gamma$dispersal$this
      
      # move_f <- terra::app(thisabund[[dispstag_f]], ismov_func, disp_prob_f  )
      move_f <- raster::calc(thisabund[[dispstag_f]],function(x) ifelse(any(!is.na(x)),rbinom(length(x),x,disp_prob_f),NA))
      stay_f <- thisabund[[dispstag_f]]-move_f
      howempty <- abundceil-(totabund-move_f)  # degree to which cells have room for dispersers
      hasmove_cells <- which(raster::getValues(move_f>0))
      if(length(hasmove_cells)>0){
        fv <- params$fv[hasmove_cells,,drop=F]
        # if(is.vector(fv)){fv<-matrix(fv,nrow=1)} #fix the case where there is only a single movecell
        dispfun_f <- dispfn(move=raster::getValues(move_f),emptiness=raster::getValues(howempty),ceil=raster::getValues(abundceil),dispmat=params$dispmat,center=params$cent,habsuit=raster::getValues(habmap3),zero=raster::getValues(zerorast))
        temp <- apply(fv,1,dispfun_f)   # do dispersal for all cells that are dispersal sources  # if(dim(fv)[1]>0){
        thisabund[[dispstag_f]] <- stay_f + raster::setValues(zerorast, get("thisdisp", env = environment(dispfun_f)) ) #thisdisp_f
        # plot(thisabund[[20]])
      }
      
    }
    
    # do matrix simulations
    totabund <- sum(thisabund)
    over <- totabund > abundceil 
    over_cells <- which(raster::getValues(over))
    
    haspop <- totabund>0
    haspop_cells <- which(raster::getValues(haspop))
    
    badmat <- thismatfunc(climvars = c(mar_mar_pr_ = 0,
                                       mar_apr_mtr_ = 0,
                                       spr_tmean_ = 0,
                                       nov_apr_pr_ = 0,
                                       ann_tmin_ = 0,
                                       feb_may_tmax_ = 0,
                                       #nov_apr_pr_norm_ = 0,
                                       july_tmean_ = 0,
                                       jun_jul_aug_tmean_ = 0,
                                       sepx2_mar_pr_ = 0,
                                       mar_nov_tmax_ = 0,
                                       spateff=0),yeareff=-0.5)*0.75
    
    if(length(haspop_cells)>0){
      if(length(over_cells)>0){
        which_over <- sapply(over_cells,function(t) which(haspop_cells==t))
      }else{
        which_over <- numeric(0)
      } 
      this_climext <- thisyr_brick[haspop_cells]
      temp <- sapply(clim_vars_unused,function(t) this_climext[[t]] <<- 0   )
      mat_list <- lapply(1:nrow(this_climext),function(t) thismatfunc(climvars=this_climext[t,],yeareff=yeareff) )
      if(length(which_over)>0) mat_list[which_over] <- lapply(1:length(which_over),function(x) badmat )
    
     prevabund <- thisabund
    
    # KTS: only apply these functions to the grid cells they correspond to
    tempdf <- prevabund[haspop_cells]
       # x=1;prev=as.numeric(tempdf[x,]);mymat=mat_list[[x]];matndx=params$MATRIX_NDX
    tofill <- t(sapply(1:nrow(tempdf), function(x) do.stoch.proj(prev=as.numeric(tempdf[x,]),mymat=mat_list[[x]],matndx=params$MATRIX_NDX) )) 
    # names(tofill) <- params$ST_NAMES
    # temp <- sapply(1:length(haspop_cells), function(x) thisabund[haspop_cells[x]] <<- tofill[x,])  
    thisabund[haspop_cells] <- tofill
   
    #all(thisabund[haspop_cells[10]]==tofill[10,])
    
    }else{
      mat_list <- NULL
    }
    
    # plot(params$countrast)
    # plot(totabund)
    # 
    
   

    
    # store data
    raster::writeRaster(thisabund,paste0(simspath,"/scenario_",scen,"_rep_",rep,"_year_",yr-1,"_",this_climmod,".tif"),overwrite=TRUE)

#get abundances
    abunds <- cellStats(thisabund,sum)
    #abunds <- terra::global(thisabund,"sum",na.rm=T)[,1]
    poparray[yr,] <- abunds #store abundances
    
    if(yr==(params$NYEARS+1)){
      thisfile <- sprintf("%s/poparray2_scen_%s_%s_rep_%s.rds",poppath,scen,this_climmod,rep)
      saveRDS(poparray,thisfile)
      }  
  }
  
  return(NULL)
}


