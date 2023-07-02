
# SCRIPT: RUN DT SIMULATIONS--------------------------------

# Notes/Decisions ---------------------------------

  # now moves climate scenario back out of the doReplicate function- climate scenario is just part of the scenarios list
  # now models all scenarios for every replicate parameter set
  # cut off all populations east of the colorado
  # assume Allison and Mcluckie number is just adults- new juveniles are now added to the init abund
  # KL added toggle so we don't need multiple versions to run on the cluster vs testing... 
    

# TODO -------------------------------

# add Ken's density info...
# add Steve Hromada survival model when Elizabeth compiles the layers... 

# I don't think we were modeling the impact of climate change appropriately. What does a 'no climate change'
#   scenario really represent. It's not necessarily 'turning off' the effect of climate drivers, since that
#   basically eliminates spatial variation. Probably makes more sense to sample randomly from
#   among the first 10 years of a climate simulation- i.e., using the burn-in as a baseline... 

# Questions ----------------------------


# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Set cluster 'toggle switch'  --------------------------------------------------------

cluster=T

# Set the run name and output paths --------------------------------------------------------

#if running on personal computer
if(!cluster){
  run_name<-"6_3_2023_test3" #set the run name
  arrayID<-1 #only one array if not on cluster
  if(!dir.exists(run_name)){
    dir.create(run_name) #make an output directory
  }#else{
  #   stop("you've already run a model with this run_name. Either delete the previous model or make a new run_name\n")
  # } 
  #make the necessary paths
  simspath<-sprintf("%s/sims_%s",run_name,arrayID)
  if(!dir.exists(simspath)) dir.create(simspath)
  scenpath<-sprintf("%s/scendata_%s",run_name,arrayID)
  if(!dir.exists(scenpath)) dir.create(scenpath)
  dumppath<-sprintf("%s/dumps_%s",run_name,arrayID)
  if(!dir.exists(dumppath)) dir.create(dumppath)
  poppath<-sprintf("%s/poparray",run_name)
  if(!dir.exists(poppath)) dir.create(poppath)
}

#if running on cluster
if(cluster){
  #get the array id from the batch job
  args <- commandArgs(trailingOnly = TRUE)
  arrayID<- as.numeric(args[1]) #array ID
  run_name<- args[2] #run_name
  print(paste("arrayID is:",arrayID))
  print(paste("run_name is:",run_name))


  #use globalscratch for file management
  globalscratch=T
  
      #folder for output files from this run
  if(globalscratch) {if(!dir.exists("/globalscratch/kjloope/")) dir.create("/globalscratch/kjloope/")} 
  if(globalscratch) {if(!dir.exists(sprintf("/globalscratch/kjloope/%s",run_name))) dir.create(sprintf("/globalscratch/kjloope/%s",run_name))}   
  
  #make one in project folder
  if(!dir.exists(run_name)) dir.create(run_name)
  
  #make a folder for each set of scenarios (ie one folder per .sh submission file)
  if(globalscratch) {simspath<-sprintf("/globalscratch/kjloope/%s/sims_%s",run_name,arrayID)} else {
    simspath<-sprintf("%s/sims_%s",run_name,arrayID)}
  dir.create(simspath)
  
  #folder for tryCatchLog dumps
  dumppath<- "/projects/birdnet/PVA_DT/Rcode_DT/dumps"
  dir.create(dumppath)
  
  #folder for output poparray data
  if(globalscratch) {poppath<-sprintf("/globalscratch/kjloope/%s/poparray",run_name)} else {
    poppath<-sprintf("%s/poparray",run_name)}
  dir.create(poppath)
  
}

# Load Functions --------------------------------------------------------

source('DT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the workspace for simulations -----------------------------

#generate params and make replist if it's the first job, otherwise wait for the first job to generate params and load it
if(arrayID==1){
  params <- prepWorkspace(nreps=120,nyears = 89)

  # Set up replicates --------------------------
  #   make the spatial random effect raster and save the params for each replicate
  #    NOTE: only run this once (prior to running simulations) and don't run again
  #          otherwise it will overwrite all the parameter information for each replicate...

  temp <- makereplist(params,params$nreps,run_name)
  save(params,file=sprintf("%s/params.RData",run_name))

} else{
  #wait for params to appear from first ArrayID run
  while (!file.exists(sprintf("%s/params.RData",run_name))) {
    Sys.sleep(10)
    print("I'm waiting for the params file to appear")
  } 
  load(sprintf("%s/params.RData",run_name)) #load the R environment with the params
  print("the new params file appeared and i've now loaded it'")
  
}

# Set up data structure to define simulation runs (to help with parallelization- can be broken up into convenient chunks)
allruns <- expand.grid(rep=1:params$nreps,scen=1:params$nscenarios)
write.csv(allruns, file=sprintf("%s/allruns.csv",run_name),row.names = F)

#define start and stop values for ndx in this run
if(cluster){
  runs_per_job<-ceiling(nrow(allruns)/24)
} else{
  allruns <- allruns[1:4,]    # for debugging
  write.csv(allruns, file=sprintf("%s/allruns.csv",run_name),row.names = F)
  runs_per_job<-nrow(allruns)
}

ndx_start<-(arrayID-1)*(runs_per_job)+1
ndx_end<-min(arrayID*runs_per_job, nrow(allruns))

# Prepare for harvesting results from TIF files ----------------------

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios

print("done with setup, now run replicates!")

# Run replicates ------------------------        
ncores=if(cluster){120} else{2}
myCluster <- makeCluster(spec=ncores,type="PSOCK")    # PREPARE FOR RUNNING PARALLEL REPLICATES
registerDoParallel(myCluster)

temp <- foreach(ndx=ndx_start:ndx_end,.packages = c('terra','raster','tryCatchLog'))%dopar% {
  
    #  params=params;run_name=run_name;rep=1;scen=1;dispersal=T;ndx=1
  
  TRYCATCH=F
  if(TRYCATCH){  
    tryCatchLog(doReplicate(params=params,run_name,
                rep=allruns$rep[ndx],
                scen=allruns$scen[ndx],
                dispersal=T),
                write.error.dump.file = TRUE,
                write.error.dump.folder=dumppath)
  }else{
    doReplicate(params=params,run_name,
                    rep=allruns$rep[ndx],
                    scen=allruns$scen[ndx],
                    dispersal=T)
  }
   
}

print("Done running scenarios...now making poparray files'")
#when it's done, make the poparray file
temp <- foreach(ndx=ndx_start:ndx_end,.packages = c('terra','raster','tryCatchLog'))%dopar% {
  
  thisscen<-allruns$scen[ndx]
  thisrep <- allruns$rep[ndx]
  thisclim <- params$allscenarios$clim_scenario[thisscen]
  
  poparray <- array(0,dim=c(nyears+1,params$NSTAGES))
  
  yr=1
  for(yr in 1:(nyears+1)){
    
    thismap <- terra::rast(sprintf("%s/scenario_%s_rep_%s_year_%s_%s.tif",simspath,thisscen,thisrep,yr-1,thisclim)) 
    
    abunds <- terra::global(thismap,"sum",na.rm=T)[,1]
    
    poparray[yr,] <- abunds
    
  }  # end year loop
  
  thisfile <- sprintf("%s/poparray_scen_%s_%s_rep_%s.rds",poppath,thisscen,thisclim,thisrep)
  saveRDS(poparray,thisfile)
}  # end foreach

stopCluster(myCluster)



#check if all of the poparray files are finished  if so, run the make data frame script
if(length(list.files(poppath, pattern = 'rds'))==nrow(allruns)){source("DT_Make_Repdf.R")}

### debugging

# load("/Users/kevinloope/Dropbox/GopherTortoise/Analysis/PVA_DT/dumps_full/dump_2023-05-12_at_00-03-59.761_PID_84043.rda")    # load the dump into the global environment
# load("/Users/kevinloope/Dropbox/GopherTortoise/Analysis/PVA_DT/dumps_full/dump_2023-05-12_at_00-04-09.896_PID_33966.rda")
# debugger(last.dump)     





## End script



















