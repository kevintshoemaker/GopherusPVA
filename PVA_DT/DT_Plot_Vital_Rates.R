# SCRIPT: Plot vital rates across specified times/places and climate spaces    --------------------------------

# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Load Functions --------------------------------------------------------

source('DT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the Workspace for simulations -----------------------------
   # make sure reps and years match the actual simuations run

run_name<-"DT_06_04_23B" #set the run name
params <- prepWorkspace(nreps=120,nyears = 89)

# Read in list of scens and reps -------------------------------

allruns <- read.csv(file=sprintf("%s/allruns.csv",run_name))

# Prepare workspace  ---------------------------

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios
nstages <- params$NSTAGES


# Determine which scenario to plot -------------------------

scen <- 1
rep <- 1
thisclim <- params$allscenarios$clim_scenario[scen]

habmap <- rast(params$HABMAP)
haspop <- which(terra::values(habmap)>0.1)
thispixel <- haspop[1]

thispar <- readRDS(sprintf("%s/replist_%s.rds",run_name,rep)) 

stagenames <- params$ST_NAMES
thisdat = params$paraminfo

# params,thispar,scen,rep,thisyears=1:89,thispixel,thisclim
toplot <- ForecastVitalRates(params,thispar,scen,rep,1:89,thispixel,thisclim)

#  plot all relationships
temp <- sapply(setdiff(names(toplot),c("year","adMCL")),function(t) plot(toplot$year,toplot[[t]],type="l",lwd=2,xlab="year",ylab=t,main=sprintf("pix %s, rep %s, scen %s",thispixel,rep,scen)) )


#  make specific plots for visualization
plotdir <- sprintf("%s/plots",run_name)
if(!dir.exists(plotdir)) dir.create(plotdir)

pdf(sprintf("%s/vr_over_time.pdf",plotdir),4.5,6)
#svg(sprintf("%s/vr_over_time.svg",plotdir),4.5,6)
par(mfrow=c(3,2))
par(mai=c(0.6,0.6,0.1,0.1))
# plot(toplot$year,zoo::rollmean( toplot$MA,10,fill=NA),type="l",lwd=2,xlab="year",ylab="Age at maturity")
plot(toplot$year,toplot$MA,type="l",lwd=2,xlab="year",ylab="Age mat.")
plot(toplot$year,toplot$CS,type="l",lwd=2,xlab="year",ylab="Rep. output")
plot(toplot$year,toplot$HS,type="l",lwd=2,xlab="year",ylab="Hatch success")
plot(toplot$year,toplot$PF,type="l",lwd=2,xlab="year",ylab="Frac. fem.")
plot(toplot$year,toplot$PR,type="l",lwd=2,xlab="year",ylab="Prob of rep.")
plot(toplot$year,toplot$lambda,type="l",xlab="year",ylab="Lambda",col="red",lwd=3)
dev.off()

# End script ---------------






