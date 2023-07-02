# SCRIPT: Plot trajectories for selected simulation scenarios   --------------------------------

#    this script makes 'spaghetti' plots for visualizing total abundance over time

# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Load Functions --------------------------------------------------------

source('DT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the Workspace for simulations -----------------------------
   # make sure reps and years match the actual simuations run

run_name<-"DT_06_04_23B" #set the run name
poppath<-sprintf("%s/poparray",run_name)
params <- prepWorkspace(nreps=120,nyears = 89)

# Read in list of scens and reps -------------------------------

allruns <- read.csv(file=sprintf("%s/allruns.csv",run_name))

# Prepare for making trajectory plots  ---------------------------

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios
nstages <- params$NSTAGES


# Determine which scenario to plot -------------------------

thisscen <- 1

thisclim <- params$allscenarios$clim_scenario[thisscen]

# Grab the relevant 'poparray' objects and combine them ---------------------

poparray <- array(0,dim=c(nreps,(nyears+1),nstages))

r=1
for(r in 1:nreps){
  poparray[r,,] <- readRDS(sprintf("%s/poparray_scen_%s_%s_rep_%s.rds",poppath,thisscen,thisclim,r))
}

##Make plot -------------------

maxabund <- max(sapply(1:nreps,function(t) max(apply(poparray[t,,],1,sum)/1000) ))*1.2


plotdir <- sprintf("%s/plots",run_name)
if(!dir.exists(plotdir)) dir.create(plotdir)

pdf(sprintf("%s/NvsTime_scen%s_%s.pdf",plotdir,thisscen,thisclim), 5,4 )

#svg(sprintf("%s/NvsTime_scen%s_%s.svg",plotdir,thisscen,thisclim), 5,4 )
plot(1,1,pch="",xlim=c(1,nyears+1),
     ylim=c(0,maxabund),
     xaxt="n",xlab="Years",ylab="Abundance (thousands)", main=paste0("Scenario ",thisscen," ",thisclim))

temp <- lapply(1:nreps,function(t) lines(1:(nyears+1),apply(poparray[t,,],1,sum)/1000,col=gray(0.6)) ) 

lines(1:(nyears+1),apply(sapply(1:nreps,function(t) apply(poparray[t,,],1,sum)),1,median)/1000,lwd=2)

axis(1,at=seq(1,nyears+1,length=5),labels = round(seq(2010,2010+nyears,length=5)))
dev.off()

  
 


