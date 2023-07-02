
# GOPHER TORTOISE SIMULATION FUNCTIONS ----------

# note: many functions are shared with analogous "DT_functions_*.R" source. If key 
#   updates are made, the corresponding update should also be made to the DT script
#   and vice versa!


# FUNCTION Explore_CHB  ----------------------------------
#    this script returns a data frame useful for visualizing critical habitat breadth
#    for a given scenario and replicate

# params=params;haspop=haspop;thispar=thispar;scen=scen;rep=rep;thisyears=1:94;thisclim=thisclim;sample=250
Explore_CHB <- function(params,haspop,thispar,scen,rep,thisyears=1:94,thisclim,sample=250){ 
  
  allpars <- names(thispar)
  thisscen <- params$allscenarios[scen,]
  
  thismatfunc <- GenerateMatrixSampler(thispar = thispar, nstages=params$NSTAGES, 
                                       stagenames = params$ST_NAMES,
                                       thisdat = params$paraminfo,
                                       thisscen = thisscen,vitals=T) 
  
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"june_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"apr_may_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_diff_","jun_jul_pr_diff_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"jun_diff_","jun_pr_diff_")
  if(thisscen$ma_clim) clim_vars_used <- c(clim_vars_used,"mn_ann_temp_")
  #cat(sprintf("%s_scen_%s_line_%s",Sys.time(),scenario,742),file="errors/errorfile.txt",append=T,sep="\n")
  
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used))
  
  burnfolder <- sprintf("%s/burnbricks",run_name) 
  
  toret=NULL
  
  yr=1
  for(yr in thisyears){ 
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",thisclim,yr) ) 
    baseT <- terra::rast(params$baseT)
    YSB <- terra::rast(sprintf("%s/YSB_%s_r%s_y%s.tif",burnfolder,thisclim,rep,yr))
    
    thisyr_brick <- c(thisyr_brick,baseT,YSB)
    names(thisyr_brick)[(terra::nlyr(thisyr_brick)-1):terra::nlyr(thisyr_brick)] <- c("baseT","YSB")
    
    thissamp <- sample(haspop,sample,replace=T)
    
    toret2 <- thisyr_brick[thissamp]
    temp <- sapply(clim_vars_unused,function(t) toret2[[t]] <<- 0   )
    
    VR_list <- lapply(1:nrow(toret2),function(t) as.data.frame(thismatfunc(climvars=toret2[t,])) )
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


ForecastVitalRates <- function(params,thispar,scen,rep,thisyears=1:94,thispixel,thisclim){ 
  
  allpars <- names(thispar)
  thisscen <- params$allscenarios[scen,]
  
  burnfolder <- sprintf("%s/burnbricks",run_name) 
  
  thismatfunc <- GenerateMatrixSampler(thispar = thispar, nstages=params$NSTAGES, 
                                       stagenames = params$ST_NAMES,
                                       thisdat = params$paraminfo,
                                       thisscen = thisscen,vitals=T) 
  
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"june_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"apr_may_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_diff_","jun_jul_pr_diff_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"jun_diff_","jun_pr_diff_")
  if(thisscen$ma_clim) clim_vars_used <- c(clim_vars_used,"mn_ann_temp_")
  #cat(sprintf("%s_scen_%s_line_%s",Sys.time(),scenario,742),file="errors/errorfile.txt",append=T,sep="\n")
  
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used))
  
  toret=NULL
  
  yr=1
  for(yr in thisyears){ 
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",thisclim,yr) ) 
    baseT <- terra::rast(params$baseT)
    YSB <- terra::rast(sprintf("%s/YSB_%s_r%s_y%s.tif",burnfolder,thisclim,rep,yr))
    
    thisyr_brick <- c(thisyr_brick,baseT,YSB)
    names(thisyr_brick)[(terra::nlyr(thisyr_brick)-1):terra::nlyr(thisyr_brick)] <- c("baseT","YSB")
    
    
    climvars <- thisyr_brick[thispixel]
    #  tempdf <- thisyr_brick[haspop]
    temp <- sapply(clim_vars_unused,function(t) climvars[[t]] <<- 0   )
    
    toret2 <- list()
    toret2$year = yr
    
    toret2 <- c(toret2,thismatfunc(climvars))
    toret <- rbind(toret,as.data.frame(toret2))
  }
  # toret
  
  return(toret)   # matrix(out,ncol=nstages)
} 


# FUNCTION "prepWorkspace" --------------------------
#    prepares the workspace with key data for GT simulation
# nreps=4; nyears=94
prepWorkspace <- function(nreps,nyears){
  # Set up workspace 
  
  params <- load_paraminfo()    # load parameter information and key R objects from file (based on google sheets info)
  # save(params,file=sprintf("%s/params.RData",run_name))
  
  # Prepare spatial data layers needed for sim, and store in 'params'
  params <- prepspatial(params,fn_abund="setup_files/baseabund.tif",
                        fn_hsi="setup_files/basehsi.tif")    # prepare spatial data
  
  # set number of reps
  params$nreps <- nreps    # should represent the total number of replicates...
  
  # Finalize Parameter list 
  params <- FinalizeParams(params,nyears=nyears,ddist=1000)
  
  # Make Scenarios 
  params$allscenarios <- makescenlist(params)   # KTS: add climate in with these scenarios
  params$nscenarios <- nrow(params$allscenarios)   # 3584 scenarios (64 scens, 4 climates)
  
  # Final preparation of climate data 
  # Break up climate data into one-year chunks (for memory management in cluster)
  # only need to do this once- can take some time...
  
  temp <- prepclimdat(params,climyears=nyears,folder="setup_files/OutputTifs",DO=F)
  
  # Do burnin simulation (if needed) -----------------------------
  # y=5;dir="setup_files/initabund_2023_05_26";useold=F    # KTS: burnin not needed for GT, not sure this makes sense for steeply declining pop!
  # params <- doBurnin(params,y=50,dir="setup_files/initabund_2023_05_26",useold=T)   # NOT IMPLEMENTED YET- RETURN TO THIS LATER!
  
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
# climyears=94;folder="setup_files/OutputTifs"; DO=T
prepclimdat <- function(params,climyears=94,folder="setup_files/OutputTifs", DO=F){
  if(DO){
    clim=1
    for(clim in 1:length(params$clim_scenarios)){
      
      this_climmod = params$clim_scenarios[clim]
      
      climvar_list <- list()
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
#  KTS: this function now stores samples as RDS files for better parallelization

makereplist <- function(params,run_name){
  
  if(!file.exists(sprintf("%s/replist_1.rds",run_name))){
    
    replist <- list()
    
    # r=1
    for(r in 1:params$nreps){
      # sample a set of parameters randomly drawn from posterior distributions, etc.
      # use this same set of parameters for all climate models and all scenarios!
      replist[[r]] <- sample_params(params,baseline=F)
      saveRDS(replist[[r]],file=sprintf("%s/replist_%s.rds",run_name,r))
      temp <- make_burnbrick(params,"noCC",r,run_name)
      c=1
      for(c in 1:length(params$clim_scenarios)){
        this_climmod <- params$clim_scenarios[c]
        YSB_brick <- make_burnbrick(params,this_climmod,r,run_name) 
      }
    }
    
  }
  
  # saveRDS(replist,file="replist.rds")
  # return(replist)
}

# FUNCTION "makescenlist" ------------------
# this function makes the list of scenarios for appending to the params object

makescenlist <- function(params){
  
  scenlist <- list()
scenlist$burn_clim <- c(T,F)    # is burn frequency dependent on climate
scenlist$phi.a_burn <- c(T,F)   # does adult survival depend on burn frequency? 
scenlist$phi.j_burn <- c(T,F)   # does juv survival depend on burn frequency?
#scenlist$cs_dep <- c(T,F)       # does climate effect on clutch size disappear at warmer latitudes TURNED OFF
#scenlist$pr_dep <- c(T,F)       # does climate effect on prob rep disappear at warmer latitudes TURNED OFF
scenlist$cs_dep <- c(F)       # does climate effect on clutch size disappear at warmer latitudes
scenlist$pr_dep <- c(F)
scenlist$cs_clim <- c(T,F)      # does clutch size change with climate change?
scenlist$pr_clim <- c(T,F)      # does rep prob depend on climate?
scenlist$hs_clim <- c(T,F)      # does hatch succ depend on climate
scenlist$pf_clim <- c(T,F)      # does sex ratio depend on climate
scenlist$ma_clim <- c(T,F)      # does age at maturity depend on temperature? 
scenlist$phi.j_short <- c(T,F)  # does juv surv jump to adult surv at half of MA?
scenlist$phi_burn_half <-c(T,F) # should we use 1/2 the slope for ysb effect on adult and juvenile mortality?
scenlist$clim_scenario <- params$clim_scenarios

allscenarios <- do.call(expand.grid,scenlist)    # 8192 scenarios to start off... 
#nrow(allscenarios)
# remove any scenarios with no climate effect but with latitude dependency  
#allscenarios <- allscenarios[!((allscenarios$cs_dep&!allscenarios$cs_clim)|(allscenarios$pr_dep&!allscenarios$pr_clim)),]

# remove any scenarios with no burn effect but an effect of burn on climate (6144 scenarios left)
allscenarios <- allscenarios[!(allscenarios$burn_clim&!(allscenarios$phi.j_burn)),]
#nrow(allscenarios)

# subset(allscenarios,burn_clim==T & phi.j_burn==F) #okay it works
# remove scenarios with juv ysb eff but no adult ysb eff? (4096 scenarios left)
allscenarios <- allscenarios[!(!allscenarios$phi.a_burn&allscenarios$phi.j_burn),]
#nrow(allscenarios)
# subset(allscenarios,phi.a_burn==F & phi.j_burn==T) #okay it works...we have only adult, or both adult and juv

#remove scenarios with half burn effect but no adult or juvenile burn effect (3584 scenarios left)
allscenarios <- allscenarios[!(allscenarios$phi_burn_half&!allscenarios$phi.j_burn&!allscenarios$phi.a_burn),]
#nrow(allscenarios)

# subset(allscenarios,phi_burn_half==T & phi.j_burn==F & phi.a_burn==F) #okay it works...we have only adult, or both adult and juv

# subset(allscenarios,phi_burn_half==T & phi.j_burn==F & phi.a_burn==F) #okay it works...we have only adult, or both adult and juv
allscenarios$index<-1:nrow(allscenarios)

#which scenarios do we want to pull out to make maps to examine various burning scenario effects
burn_scen_tifs<-subset(allscenarios,cs_dep==F & pr_dep==F & cs_clim==T & pr_clim==T & hs_clim==T & pf_clim==T & ma_clim==T & phi.j_short==T ) 
#these are in sims_1 and sims_11

no_clim_tifs<-subset(allscenarios,cs_dep==F & pr_dep==F & cs_clim==F & pr_clim==F & hs_clim==F & pf_clim==F & ma_clim==F & phi.j_short == T ) 
no_clim_tifs$group<-floor(no_clim_tifs$index/42)+1

write.csv(burn_scen_tifs,sprintf("%s/burn_scen_tifs.csv",run_name),row.names=F)
write.csv(no_clim_tifs,sprintf("%s/no_clim_tifs.csv",run_name),row.names=F)

allscenarios
  
  # NSCENARIOS <- nrow(allscenarios)
  return(allscenarios)
}


# FUNCTION "prepspatial" ------------------------
# this function reads in desert tortoise initial habitat suitability map, spreads the total 
#  abundance across all cells and writes a habitat map to file as geotiff
#  this function also ensures that key spatial data are bundled in 'params'

# fn_abund="setup_files/baseabund.tif";fn_hsi="setup_files/basehsi.tif"
prepspatial <- function(params,fn_abund="setup_files/baseabund.tif",fn_hsi="setup_files/basehsi.tif"){

# habitat suitability-- maybe not needed??
  hsimap <-  terra::rast(fn_hsi) 
  hsimap <-  hsimap/(terra::global(hsimap,"max",na.rm=T)["layer","max"])   # Note: base.habmap map comes from the data bundle
  # should be max of 1. Min should be close to zero
  
# abundance  
  base.abund <- terra::rast(fn_abund)
  base.abund <- ifel(base.abund<=2,0,base.abund)
  base.abund <- base.abund/2  #50:50 sex ratio
  
  habmap <- base.abund/terra::global(base.abund,"max",na.rm=T)["layer","max"]
  
  params$CRS <- terra::crs(habmap,proj=T)
  params$RESOLUTION <- terra::res(habmap)[1]     # in m per side
  params$CELL_AREA_KM <- (params$RESOLUTION/1000)^2   # in km2
  
  initabund <- base.abund
  
  mat <- terra::rast("setup_files/mat.tif")  # do we need this?
  
  #### other global params
  
  params$max_density <- 1/0.01    #0.00404686   # convert to per km (FOLT: 2 female tortoises per acre) (I made it half that to allow for heterogeneity within 25km2) (note: it looks to be ha not acres)
  
  # Burning related objects  
  
  stacknames <- list.files("setup_files/OutputTifs")[grepl("burn.yr.change",list.files("setup_files/OutputTifs"))]
  stacknames <- stacknames[grepl(".tif",stacknames)]
  stacknames_full <- paste0("setup_files/OutputTifs/",stacknames)
  clim_scenarios <- sapply(strsplit(stacknames,"_"),function(t) t[4])
  params$clim_scenarios <- gsub(".tif","",clim_scenarios[grepl("tif",clim_scenarios)])
  
  names(stacknames_full) <- params$clim_scenarios
  params$burnprob_bricks <- stacknames_full # now packaged with climate vars
  
  # write habmap to file and record the name
  params$HABMAP <- "setup_files/habmap.tif"
  terra::writeRaster(habmap,params$HABMAP,overwrite=TRUE)   # write map to working directory
  
  # derived spatial params objects
  habmap2 <- terra::ifel(is.na(habmap),0,habmap)   # same, but with zeros instead of NAs
  habmap3 <- terra::ifel(habmap2==0,0.0001,habmap2)    # needed for dispersal algorithm
  habmap3[habmap3==0] <- 0.0001
  temp <- terra::rast(habmap)
  zerorast <- terra::init(temp,0) # all zeros
  countrast <- terra::init(temp,fun="cell")  # raster value is cell index
  countrast <- terra::ifel(is.na(habmap),NA,countrast)  # make non-habitat cells NA to ensure no dispersal to those cells
  # plot(params$countrast)
  onerast <- terra::init(temp,1)
  
  # write to file and store filenames
  
  params$HABMAP2 <- "setup_files/habmap2.tif"
  terra::writeRaster(habmap2,params$HABMAP2,overwrite=TRUE)
  params$HABMAP3 <- "setup_files/habmap3.tif"
  terra::writeRaster(habmap3,params$HABMAP3,overwrite=TRUE)
  params$zerorast <- "setup_files/zerorast.tif"
  terra::writeRaster(zerorast,params$zerorast,overwrite=TRUE)
  params$countrast <- "setup_files/countrast.tif"
  terra::writeRaster(countrast,params$countrast,overwrite=TRUE)
  params$onerast <- "setup_files/onerast.tif"
  terra::writeRaster(onerast,params$onerast,overwrite=TRUE)
  
  ## set age/stage structure in advance of setting initial abundance
  
  params$NSTAGES <- 20     # hardcode 20 stages for now
  juvstages <- params$NSTAGES-1
     # approximate stable age distribution
  initdist <- c(seq(0.125,0.05,length=juvstages),5) # fraction of individuals in each age class
  # initdist[20]/sum(initdist[1:19])
  # abundances per cell, by stage 
  INITABUND2 <- initdist/sum(initdist)
  fracadult <- INITABUND2[params$NSTAGES]

  params$ABUNDCEIL_MAX <- params$max_density*params$CELL_AREA_KM  # abundance ceiling per cell in high-quality habitat 
 
  #  compute maximum/ceiling abundance for "density dependence"
  
  abundceil <-initabund * 2  # assume all cells could support twice as many as curren
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

  params$baseT <- "setup_files/mat.tif"
  
  return(params)
}

# FUNCTION "DoBurnin --------------------------

# run simulation across the landscape under baseline conditions
   #  to ensure that initial abundance can be used as a reference.
   #  returns an updated stable initial abundance raster

doBurnin <- function(params,y,dir="setup_files/initabund_2023_05_26", useold=T){
  
  if(!dir.exists(dir)) dir.create(dir)
  ceil_filename <- sprintf("%s/ceiling.tif",dir)
  params$ABUNDCEIL <- ceil_filename
  
  if(useold){
    filelist <- sprintf("%s/initabund_%s.tif",dir,params$clim_scenarios)
  }else{
    filelist = NULL
    thisparams <- sample_params(params,baseline=T)
    
    thisscen <- params$allscenarios[1,]
    
    # generate a function for generating the matrix
    # thispar = thisparams; nstages=params$NSTAGES; stagenames = params$ST_NAMES;thisdat = params$paraminfo;vitals=T
    thismatfunc <- GenerateMatrixSampler(thispar = thisparams, nstages=params$NSTAGES, 
                                         stagenames = params$ST_NAMES,
                                         thisdat = params$paraminfo,
                                         thisscen = thisscen,vitals=F) 
    
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
      
      yr=2
      for(yr in 2:y){   
        
        yeareff = rnorm(1)   # account for good years and bad years...
        
        # read brick of climate data for this year
        thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",this_climmod,sample(1:10,1)) ) 
        
        thisyr_brick <- c(thisyr_brick,spateff)
        names(thisyr_brick)[terra::nlyr(thisyr_brick)] <- "spateff"
        
        # totabund <- terra::app(thisabund,sum)  # compute total abundance
        totabund <- sum(thisabund)
        over <- totabund > initabund_tot*1.3   # exceeds limit
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
  thiscsv <- read.csv("setup_files/InitialDistributions_GT.csv")
  thiscsv <- thiscsv[,1:(ncol(thiscsv)-1)]
  thiscsv <- thiscsv[thiscsv$Used=="Y",]
  # head(thiscsv)
  
  thisrobj <- read.csv("setup_files/RObjects_GT.csv")
  thisrobj <- thisrobj[,1:(ncol(thisrobj)-1)]
  thisrobj <- thisrobj[thisrobj$Used=="Y",]
  # thisrobj
  
  load("setup_files/RObjects.RData")    # load R objects
  
  datalist$allparams <- sort(unique(c(thiscsv$Parameter,thisrobj$Parameter)))
  temp <- sapply(1:length(datalist$allparams),function(t)  datalist$paraminfo[[datalist$allparams[t]]] <<- list() )
  
  r=1
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


# FUNCTION "MakeLHSSamples"  (for global sensitivity testing) -----------------------

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

# FUNCTION "specifyLHSParam"   (for global sensitivity testing) -----------------------

# Information necessary to translate standard uniform LHS sample into parameters of interest for paleo project 

specifyLHSParam <- function(paramslist,name,type,lb,ub){
  newlist <- paramslist
  eval(parse(text=sprintf("newlist$%s <- list()",name)))
  eval(parse(text=sprintf("newlist$%s$type <- \"%s\"",name,type)))
  eval(parse(text=sprintf("newlist$%s$lb <- %s",name,lb)))
  eval(parse(text=sprintf("newlist$%s$ub <- %s",name,ub))) 	
  return(newlist)
}

# FUNCTION Geometric mean function --------------------

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

# FUNCTION 'make_burnbrick'----------------------------- 
#   generate a brick of YSB for all years for a given scenario/replicate
# this_climmod="CES45"; rep=1
make_burnbrick <- function(params, this_climmod, rep, run_name){
  # set initial years since burn (YSB)
  
  # initialize years since burn
  hab <- terra::rast(params$HABMAP)
  YSB <- terra::ifel(is.na(hab),NA,0)
  tofill_ndx <- which(!is.na(terra::values(YSB)[,1]))
  tofill <- sample(0:3,length(tofill_ndx),replace=T)
  YSB[tofill_ndx] <- tofill
    # plot(YSB)
  
  # get burn probs if needed
  if(this_climmod%in%params$clim_scenarios){
    burnprob_brick <- rast(params$burnprob_bricks[this_climmod])
    # plot(burnprob_brick$index_2)
  }
  thisfolder <- sprintf("%s/burnbricks",run_name) 
  if(!file.exists(thisfolder)) dir.create(thisfolder)
  
  # terra::writeRaster(YSB,filename =sprintf("%s/YSB_%s_r%s_y%s.tif",thisfolder,this_climmod,rep,0) ,overwrite=T)
  
  yr=2
  for(yr in 2:(params$NYEARS+1)){
    # update YSB raster
    if(this_climmod%in%params$clim_scenarios){
      these_burnprobs <- terra::values(burnprob_brick[[yr-1]])[tofill_ndx]*params$BURNPROB
    }else{
      these_burnprobs <- params$BURNPROB    # alternative if burn-climate link is OFF
    }
    
    didburn <- rbinom(length(tofill_ndx),1,these_burnprobs)
    these_YSB <- terra::values(YSB)[tofill_ndx,1]
    these_YSB[as.logical(didburn)] <- 0
    these_YSB[!as.logical(didburn)] <- these_YSB[!as.logical(didburn)] + 1
    these_YSB <- ifelse(these_YSB>4,4,these_YSB)
    YSB[tofill_ndx] <- these_YSB
    terra::writeRaster(YSB,filename =sprintf("%s/YSB_%s_r%s_y%s.tif",thisfolder,this_climmod,rep,yr-1) ,overwrite=T)
  }
  # toret <- raster::brick(YSB_brick)
  # names(toret) <- paste0("year",0:(params$NYEARS))
  # plot(toret$year0)
  return(1)
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
  
  TEST=F
  # TEST=T
  if(TEST){
    temp <- rast(sprintf("setup_files/OutputTifs/%s_1.tif",params$clim_scenarios[1]))
    climvars <- temp[1]
    temp <- lapply(params$clim_var_names,function(t) climvars[[t]] <<- 0)
    climvars$baseT = 20
    climvars$YSB = 1
  }
  
  thisNS = max(0.001,min(.99,thispar$NS$nestsucc$this))
  thisPHI_H = max(0.05,min(1,thispar$PHI_H$hatchsurv$this)) 
  
  matrix_sampler <- function(climvars){
    
    thisMA = max(9,min(nstages,round(thispar$MA$ma.intercept.mn$this + thispar$MA$ma.tempeff.mn$this * (climvars[["baseT"]] + (climvars[["mn_ann_temp_"]]*thisscen$ma_clim) ) ) ))
    if(vitals) toret$MA=thisMA
    
    thisPHI_J = max(0.5,min(0.99,thispar$PHI_J$juvsurv$this + thispar$PHI_A$phi.a.burneff.mn$this * min(4,climvars[["YSB"]]) * max(0,(thisscen$phi.j_burn-0.5*thisscen$phi_burn_half)) ) ) 
    if(vitals) toret$PHI_J=thisPHI_J
    
    thisPHI_A = max(0.75,min(0.99, thispar$PHI_A$adultsurv.mn$this + thispar$PHI_A$phi.a.burneff.mn$this * min(4,climvars[["YSB"]]) * max(0,(thisscen$phi.a_burn-0.5*thisscen$phi_burn_half)) )  )
    if(vitals) toret$PHI_A = thisPHI_A
    
    thisCSmean = max(0.01,thispar$CS$cs.intercept.mn$this + thispar$CS$cs.tempeff.mn$this * climvars[["baseT"]] )  # set baseline no climate change
    thisCS = max(0.01,thisCSmean + thispar$CS$CS_intercept$this + 
                   (thispar$CS$CS_junetmax$this*
                      (climvars[["june_"]])/thisdat$CS$junetmax.sd) * thisscen$cs_clim)  # climvars[["thisatten_cs"]]
    if(vitals) toret$CS = thisCS
    
    thisHSmod = thisdat$HS$HS_clim_eff$func(deltaT=climvars[["jun_jul_diff_"]]*thisscen$hs_clim,
                                            deltaP=climvars[["jun_jul_pr_diff_"]]*thisscen$hs_clim,
                                            hs.tempeff.draw=thispar$HS$HS_nest_tempeff$this,
                                            hs.precipeff.draw = thispar$HS$HS_nest_precipeff$this,
                                            means.nest=thisdat$HS$HS_meansnest$obj)
    thisHS =  thispar$HS$hs.mn$this*thisHSmod
    if(vitals) toret$HS = thisHS
    
    # if(vitals) toret$NS = thisNS
    # if(vitals) toret$PHI_H = thisPHI_H
    
    thisPFmod = thisdat$PF$PF_clim_eff$func(deltaTjune = climvars[["jun_diff_"]]*thisscen$pf_clim,
                                            deltaPjune = climvars[["jun_pr_diff_"]]*thisscen$pf_clim,
                                            hs.tempeff.draw=thispar$PF$PF_nest_tempeff$this,
                                            hs.precipeff.draw = thispar$PF$PF_nest_precipeff$this,
                                            t.mean.ref = thisdat$PF$t.mean.ref,
                                            p.mean.ref = thisdat$PF$p.mean.ref,
                                            sex.tmean.pmean.top = thisdat$PF$PF_clim_eff_model$obj)
    
    thisPF = thisPFmod["perc.fem"] 
    if(vitals) toret$PF = thisPF
    
    thisPR = plogis(thispar$PR$PR_intercept$this + 
                      (thispar$PR$PR_aprmaytmax$this*
                         (climvars[["apr_may_"]])/thisdat$PR$aprmaytmax.sd)*thisscen$pr_clim   # climvars[["thisatten_pr"]]  
                    +1.8) # correction for timing of nest discovery   
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
FinalizeParams <- function(params,nyears=2,ddist=1000,burnprob=0.25){
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
  
  params$clim_var_names <- c( 
    "apr_may_",
    "jun_diff_",
    "jun_jul_diff_",
    "jun_jul_pr_diff_" ,
    "jun_pr_diff_",
    "june_",  
    "mn_ann_temp_"
    # "burn_yr_change_"   # note: this doesn't actually belong here- only used to create "burnbricks"
  )
  
  params$BURNPROB = burnprob

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

  abundceil <- raster::raster(params$ABUNDCEIL)
  countrast <- raster::raster(params$countrast)  # maybe not needed?
  zerorast <- raster::raster(params$zerorast)
  habmap3 <- raster::raster(params$HABMAP3)
  
  thisscen <- params$allscenarios[scen,]   # this scenario specs
  
  burnfolder <- sprintf("%s/burnbricks",run_name) 
  
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
  thisabund <- raster::brick(params$INITABUND)  # KTS- if no burnin, not initialized differently by climate scen
  names(thisabund) <- params$ST_NAMES
  if(!bigfile){
  raster::writeRaster(thisabund,paste0(simspath_tmp,"/scenario_",scen,"_rep_",rep,"_year_",0,"_",this_climmod,".tif"),overwrite=TRUE)
  }
  
  if(bigfile){
    abunds_out<-thisabund
  names(abunds_out)<-paste("Y_0.",names(abunds_out),sep="")
  }
  #make fresh poparray array if loop is starting
  poparray<-array(0,dim=c(nyears+1,params$NSTAGES))
  abunds <- cellStats(thisabund,sum)
  #abunds <- raster::global(thisabund,"sum",na.rm=T)[,1]
  poparray[1,] <- abunds #store abundances
  
    
  # determine which climate vars are used
  clim_vars_used <- c()
  if(thisscen$cs_clim) clim_vars_used <- c(clim_vars_used,"june_")
  if(thisscen$pr_clim) clim_vars_used <- c(clim_vars_used,"apr_may_")
  if(thisscen$hs_clim) clim_vars_used <- c(clim_vars_used,"jun_jul_diff_","jun_jul_pr_diff_")
  if(thisscen$pf_clim) clim_vars_used <- c(clim_vars_used,"jun_diff_","jun_pr_diff_")
  if(thisscen$ma_clim) clim_vars_used <- c(clim_vars_used,"mn_ann_temp_")
  #cat(sprintf("%s_scen_%s_line_%s",Sys.time(),scenario,742),file="errors/errorfile.txt",append=T,sep="\n")
  
  clim_vars_used <- sort(unique(clim_vars_used))
  clim_vars_unused <- setdiff(params$clim_var_names,c(clim_vars_used))
  
  # all_var_names <- c(clim_vars_used,clim_vars_unused,"baseT","YSB")
  
  temp <- rast(sprintf("setup_files/OutputTifs/%s_1.tif",params$clim_scenarios[1]))
  climvars2 <- temp[1]
  temp <- lapply(params$clim_var_names,function(t) climvars2[[t]] <<- 0)
  climvars2$baseT = 20
  climvars2$YSB = 3
  
  yr=2
  for(yr in 2:(params$NYEARS+1)){
    
    # read brick of climate data for this year
    thisyr_brick <- terra::rast(sprintf("setup_files/OutputTifs/%s_%s.tif",this_climmod,yr-1) ) 
    
    baseT <- terra::rast(params$baseT)
    YSB <- terra::rast(sprintf("%s/YSB_%s_r%s_y%s.tif",burnfolder,this_climmod,rep,yr-1))
    
    thisyr_brick <- c(thisyr_brick,baseT,YSB)
    names(thisyr_brick)[(terra::nlyr(thisyr_brick)-1):terra::nlyr(thisyr_brick)] <- c("baseT","YSB")
    
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
    
    badmat <- thismatfunc(climvars = climvars2)*0.75
    
    if(length(haspop_cells)>0){
      if(length(over_cells)>0){
        which_over <- sapply(over_cells,function(t) which(haspop_cells==t))
      }else{
        which_over <- numeric(0)
      } 
      this_climext <- thisyr_brick[haspop_cells]
      temp <- sapply(clim_vars_unused,function(t) this_climext[[t]] <<- 0   )
      mat_list <- lapply(1:nrow(this_climext),function(t) thismatfunc(climvars=this_climext[t,]) )
      if(length(which_over)>0) mat_list[which_over] <- lapply(1:length(which_over),function(x) badmat )
      
      prevabund <- thisabund
    
    # KTS: only apply these functions to the grid cells they correspond to
    tempdf <- prevabund[haspop_cells]
       # x=1;prev=as.numeric(tempdf[x,]);mymat=mat_list[[x]];matndx=params$MATRIX_NDX
    tofill <- t(sapply(1:nrow(tempdf), function(x) do.stoch.proj(prev=as.numeric(tempdf[x,]),mymat=mat_list[[x]],matndx=params$MATRIX_NDX) )) 
    # names(tofill) <- params$ST_NAMES
    # temp <- sapply(1:length(haspop_cells), function(x) thisabund[haspop_cells[x]] <<- tofill[x,])  
    thisabund[haspop_cells] <- tofill
    
    }else{
      mat_list <- NULL
    }
    
    # plot(params$countrast)
    # plot(totabund)
    # 
    

   
    #all(thisabund[haspop_cells[10]]==tofill[10,])

    
    # store data
    if(!bigfile){
    raster::writeRaster(thisabund,paste0(simspath_tmp,"/scenario_",scen,"_rep_",rep,"_year_",yr-1,"_",this_climmod,".tif"),overwrite=TRUE)
    }
    #get abundances
    abunds <- cellStats(thisabund,sum)
    #abunds <- terra::global(thisabund,"sum",na.rm=T)[,1]
    poparray[yr,] <- abunds #store abundances
    #names(thisabund)<-paste(yr-1,"_0-",names(thisabund),sep="")
    if(bigfile){
    abunds_out<-stack(abunds_out,thisabund)
    nnames<-length(names(abunds_out))
    names(abunds_out)[(nnames-19):nnames]<-paste("Y_",yr-1,".",names(thisabund),sep="")
    }
    
    if(yr==(params$NYEARS+1)){
      thisfile <- sprintf("%s/poparray2_scen_%s_%s_rep_%s.rds",poppath,scen,this_climmod,rep)
      saveRDS(poparray,thisfile)
      if(bigfile){
      raster::writeRaster(abunds_out,paste0(simspath_tmp,"/scenario_",scen,"_rep_",rep,"_all_",this_climmod,".tif"),overwrite=TRUE)
      write.csv(data.frame(names=names(abunds_out)),paste0(simspath_tmp,"/scenario_",scen,"_rep_",rep,"_all_",this_climmod,"_names.csv"))
      }
    }
    
  
    }
  
  return(NULL)
}

90*94
