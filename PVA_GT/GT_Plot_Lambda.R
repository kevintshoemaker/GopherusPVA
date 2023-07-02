# SCRIPT: Plot lambda across climate space    --------------------------------

# Clear workspace --------------------------------------------------------

rm(list=ls())   # clear workspace

# Load Functions --------------------------------------------------------

source('GT_functions.R')   # load functions from file

# Load Packages -----------------------------------------------------------

suppressWarnings(loadPackages())

# Set up the Workspace for simulations -----------------------------
   # make sure reps and years match the actual simuations run

run_name<-"GT_06_04_23" #set the run name
params <- prepWorkspace(nreps=120,nyears = 94)

# Read in list of scens and reps -------------------------------

allruns <- read.csv(file=sprintf("%s/allruns.csv",run_name))

# Prepare workspace  ---------------------------

nreps <- max(allruns$rep)
nyears <- params$NYEARS
nclims <- length(params$clim_scenarios)

allcmods <- params$clim_scenarios
allscens <- params$allscenarios
nstages <- params$NSTAGES


# Determine which scenario/clim to plot -------------------------

scen <- 1
rep <- 2
thisclim <- params$allscenarios$clim_scenario[scen]

habmap <- rast(params$HABMAP)
haspop <- which(terra::values(habmap)>0.1)
thispixel <- haspop[1]

thispar <- readRDS(sprintf("%s/replist_%s.rds",run_name,rep)) 

stagenames <- params$ST_NAMES
thisdat = params$paraminfo

# params=params;haspop=haspop;thispar=thispar;scen=scen;rep=rep;thisyears=1:89;thisclim=thisclim;sample=250
thisdf <- Explore_CHB(params,haspop,thispar,scen,rep,thisyears=1:94,thisclim,sample=250)

loadings <- thisdf[[2]]
thisdf <- thisdf[[1]]

## visualize results as heatmap or contour plot

factors_toplot <- c(1,2)
facnames <- paste0("Factor",factors_toplot)

keep <- c(facnames,"lambda")

df <- thisdf[,keep]

bins1 <- seq(min(df[,1]),max(df[,1]),length=25)
binsize1 <- bins1[20]-bins1[19]
bins2 <- seq(min(df[,2]),max(df[,2]),length=25)
binsize2 <- bins2[20]-bins2[19]

df_p <- expand.grid(bins1,bins2)
df_p$lambda <- NA

r=100
for(r in 1:nrow(df_p)){
  this1 <- df_p[r,1]
  this2 <- df_p[r,2]
  ndx1 <- which(near(df[,1],this1,binsize1)) 
  ndx2 <- which(near(df[,2],this2,binsize1))
  ndx <- intersect(ndx1,ndx2)
  if(length(ndx)>0){
    df_p$lambda[r] <- mean(df$lambda[ndx],na.rm=T)
  }
}

df_p$lambda2 <- cut(df_p$lambda,breaks=c(-Inf,0.99,0.995,1.005,1.01,Inf),labels=c("strong neg","weak neg","stable","weak pos","strong pos"))

library(ggplot2)

plotdir <- sprintf("%s/plots",run_name)
if(!dir.exists(plotdir)) dir.create(plotdir)

pdf(sprintf("%s/CHBfig1.pdf",plotdir),5.5,4)
#svg(sprintf("%s/CHBfig1.svg",plotdir),5.5,4)
thisplot <- ggplot(df_p, aes(Var1, Var2, fill= lambda2)) +
  geom_tile() +
  #scale_fill_gradient(low="red", high="green")   # +
  scale_fill_brewer(palette = "RdYlGn") +
  labs(x="Factor 1 (temperature)",y="Factor 2 (precip)",fill="Lambda",title=sprintf("DT, rep %s, %s",rep, thisclim))
  
print(thisplot)
dev.off()





# End script ---------------






