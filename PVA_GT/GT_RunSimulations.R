
# SCRIPT: RUN GT SIMULATIONS--------------------------------

# Notes/Decisions ---------------------------------
    

# TODO -------------------------------


# Questions ----------------------------


# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Set cluster and bigfile and storage 'toggle switches'  --------------------------------------------------------

cluster=T
bigfile=F
#use temp storage for faster writing
#tempstor<-"tmpdir"
#tempstor<-"tmpfs"
tempstor<-"none"

#put stuff on global scratch?
globalscratch=T

#prep new params and burn maps?  otherwise just use ones in folder
prepnew=F

ncores=if(cluster){120} else{2}

# Set the run name and output paths --------------------------------------------------------

#if running on personal computer
if(!cluster){
  run_name<-"06_02_23" #set the run name
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

#get tmpdir and tmpfs
  tmpdir<- Sys.getenv("TMPDIR")
  tmpfs<-Sys.getenv("TMPFS")
  
      #folders for output files from this run
   if(globalscratch) {if(!dir.exists(sprintf("/globalscratch/kjloope/%s",run_name))) dir.create(sprintf("/globalscratch/kjloope/%s",run_name))} else {  
  if(!dir.exists(run_name)) dir.create(run_name)} 
  if(tempstor=="tmpdir") {if(!dir.exists(sprintf("%s/%s",tmpdir,run_name))) dir.create(sprintf("%s/%s",tmpdir,run_name))}
  if(tempstor=="tmpfs") {if(!dir.exists(sprintf("%s/%s",tmpfs,run_name))) dir.create(sprintf("%s/%s",tmpfs,run_name))}
    
  
  #make a folder for each set of scenarios (ie one folder per .sh submission file)
  if(globalscratch) {simspath<-sprintf("/globalscratch/kjloope/%s/sims_%s",run_name,arrayID)} else {
    simspath<-sprintf("%s/sims_%s",run_name,arrayID)}
  dir.create(simspath)
  
  #create a temporary storage for writing files before compressing and transferring to long-term storage
  if(tempstor=="tmpdir") {simspath_tmp<-sprintf("%s/%s/sims_%s",tmpdir,run_name,arrayID)}
  if(tempstor=="tmpfs") {simspath_tmp<-sprintf("%s/%s/sims_%s",tmpfs,run_name,arrayID)}
  if(tempstor=="none") {simspath_tmp<-simspath}#if not using tempstor, just make simspath_tmp the same as simspath
  if(!dir.exists(simspath_tmp)) {dir.create(simspath_tmp)} 
  
  #folder for tryCatchLog dumps
  dumppath<- sprintf("%s/dumps_%s",run_name,arrayID)
  dir.create(dumppath)
  
  #folder for output poparray data
  if(globalscratch) {poppath<-sprintf("/globalscratch/kjloope/%s/poparray",run_name)} else {
    poppath<-sprintf("%s/poparray",run_name)}
  dir.create(poppath)
  
}

# Load Functions --------------------------------------------------------

source('GT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the workspace for simulations -----------------------------
#note: KJL parallelized this to replace the makereplist function since it was taking ~5 hrs to generate replists for 120 reps

#generate params and make replist if it's the first job, otherwise wait for the first job to generate params and load it
if(arrayID==1 & prepnew==T){
  params <- prepWorkspace(nreps=112,nyears = 94)
  # Set up replicates --------------------------
  #   make the 'years since burn' rasters and save the specific params for each replicate
  #    NOTE: only run this once (prior to running simulations) and don't run again
  #          otherwise it will overwrite all the parameter information for each replicate...

 # temp <- makereplist_par(params,run_name,ncores=ncores)   # takes a while- also makes the YSB rasters... could parallelize this easily if we want to... 
  replist <- list()
  
  myCluster <- makeCluster(spec=ncores,type="PSOCK")    # PREPARE FOR RUNNING PARALLEL
  registerDoParallel(myCluster)
  
  #r=1
  temp <- foreach(r=1:params$nreps,.packages = c('terra','raster','tryCatchLog','mvtnorm'))%dopar% {
    
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
  stopCluster(myCluster)
  
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
  runs_per_job<-ceiling(nrow(allruns)/240)
} else{
  allruns <- allruns[1:4,]    # for debugging
  write.csv(allruns, file=sprintf("%s/allruns.csv",run_name),row.names = F)
  runs_per_job<-nrow(allruns)
}

ndx_start<-(arrayID-1)*(runs_per_job)+1
ndx_end<-min(arrayID*runs_per_job, nrow(allruns))

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios

print("done with setup, now run replicates!")

# Run replicates ------------------------        

myCluster <- makeCluster(spec=ncores,type="PSOCK")    # PREPARE FOR RUNNING PARALLEL REPLICATES
registerDoParallel(myCluster)

temp <- foreach(ndx=ndx_start:ndx_end,.packages = c('terra','raster','tryCatchLog'))%dopar% {
  
    #  params=params;run_name=run_name;rep=1;scen=1;dispersal=T;ndx=1
  
  TRYCATCH=T
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
files<-list.files(simspath_tmp)
runscens<-as.numeric(unique(sapply(strsplit(files,"_"), `[`, 2)))
mywd<-getwd()


print("Done running scenarios...now compress and copy the files out of tempstor")
temp <- foreach(myscen=runscens)%dopar% {
setwd(simspath_tmp)
thesefiles<-files[grep(paste0("scenario_",myscen,"_"),files)]
tar(sprintf("%s/sims_%s_scen_%s.tar.bz2",simspath,arrayID,myscen), files = thesefiles,
    compression = "bzip2")
file.remove(thesefiles) #now delete the untarred individual files
}

stopCluster(myCluster)
setwd(mywd)
#check if all of the poparray files are finished  if so, run the make data frame script
print(paste("number of poparray files is:",length(list.files(poppath, pattern = 'rds'))))
print(paste("number of rows of allruns is:",nrow(allruns)))

if(length(list.files(poppath, pattern = 'rds'))==nrow(allruns)){source("GT_Make_Repdf.R")}

### debugging

# load("/Users/kevinloope/Dropbox/GopherTortoise/Analysis/PVA/Rcode/troubleshooting/dump_2023-06-01_at_13-23-43.279_PID_243911.rda")    # load the dump into the global environment
# load("/Users/kevinloope/Dropbox/GopherTortoise/Analysis/PVA_DT/dumps_full/dump_2023-05-12_at_00-04-09.896_PID_33966.rda")
# debugger(last.dump)     





## End script



















