# SCRIPT: PREPARE ANIMATIONS FOR GT MODEL   --------------------------------
#    this script prepares images for animation

# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Load Functions --------------------------------------------------------

source('GT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the Workspace for simulations -----------------------------
   # make sure reps and years match the actual simuations run

run_name<-"5_27_2023_test" #set the run name
arrayID=1
simspath<-sprintf("%s/sims_%s",run_name,arrayID)
animpath<-sprintf("%s/animate",simspath)
if(!dir.exists(animpath)) dir.create(animpath)
params <- prepWorkspace(nreps=4,nyears = 94)

# Read in list of scens and reps -------------------------------

allruns <- read.csv(file=sprintf("%s/allruns.csv",run_name))

# Prepare for storing plots for animations  ---------------------------

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios

# Define what scen/rep/clim we want to run animation for ----------------

thisrep <- 1
thisscen <- 1

thisscen_name <- "All climate effects" ##paste0("scen_",thisscen)

thisclim <- params$allscenarios$clim_scenario[thisscen]

onlyAD <- F

# Run through years and store frames for animation -----------------------

yr=1
for(yr in 1:(nyears+1)){
  thismap <- terra::rast(sprintf("%s/scenario_%s_rep_%s_year_%s_%s.tif",simspath,thisscen,thisrep,yr-1,thisclim))
      
  if(onlyAD){
    totabund <- thismap[[params$NSTAGES]]
  }else{
    totabund <- terra::app(thismap,sum)
  }

  totabund <- ifel(totabund<20,0,totabund)
  tiff(sprintf("%s/scen_%s_rep_%s_year_%s_clim_%s.tif",animpath,thisscen,thisrep,yr,thisclim),
      width=5,height=5,res=150,units="in")
  
    plot(totabund,main=sprintf("rep %s, year %s, clim %s", thisrep, yr-1, thisclim),range=c(0,100) )
  dev.off()
        
      
} 



# End script  --------------------------
    
  
 


