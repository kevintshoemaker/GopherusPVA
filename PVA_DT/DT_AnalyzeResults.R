# ANALYZE SIMULATION RESULTS- DESERT TORTOISE
#  Use machine learning to analyze replicate-level results

# Clear workspace ----------------------

rm(list=ls())   # clear workspace

   # NOTE: don't need to load the DT_functions script here

# Load packages -----------------


# if(!require("tidyverse")) install.packages("tidyverse")
# if(!require("caret")) install.packages("caret")
# if(!require("ranger")) install.packages("ranger")
# if(!require("ggplot2")) install.packages("ggplot2")
# if(!require("dplyr")) install.packages("dplyr")
# if(!require("party")) install.packages("party")
# if(!require("faux")) install.packages("faux")
# if(!require("doParallel")) install.packages("doParallel")

library(tidyverse)
library(caret)
library(ranger)
library(ggplot2)
library(dplyr)
library(party)
library(faux)
library(doParallel)

# myCluster <- makeCluster(spec=120,type="PSOCK")    # PREPARE FOR RUNNING PARALLEL REPLICATES
# registerDoParallel(myCluster)

# library(randomForest)

# Load functions -------------------

sigfig <- function(vec, n=3){ 
  ### function to round values to N significant digits
  # input:   vec       vector of numeric
  #          n         integer is the required sigfig  
  # output:  outvec    vector of numeric rounded to N sigfig
  
  formatC(signif(vec,digits=n), digits=n,format="fg", flag="#") 
  
}      # end of function   sigfig

   # compare key vital rates between original distribution and distribution for sustainable sims
orig_vs_sustainable <- function(obj,var){
  # obj <- sust_list
  # var <- "PHI_A_base"
  thisobj <- obj[[var]]
  thisdf <- data.frame(
    data = c(thisobj$orig,thisobj$sust),
    cat = c(rep("orig",times=length(thisobj$orig)),rep("sust",times=length(thisobj$sust)))
  )
  a<- ggplot(thisdf) +
    geom_density(aes(x=data,fill=cat),alpha=0.5) +
    labs(title=var)
  print(a)
}

   # compare scenario frequency between original and sustainable 
orig_vs_sustainable2 <- function(obj,var){
  # obj <- sust_list2
  # var <- "burn_clim"
  thisobj <- obj[[var]]
  thisdf <- data.frame(
    data = c(thisobj$orig,thisobj$sust),
    cat = c(rep("orig",times=length(thisobj$orig)),rep("sust",times=length(thisobj$sust)))
  )
  
  thisdf2 <- thisdf %>% 
    group_by(cat) %>% 
    summarize(frac = mean(data))
  
  a<- ggplot(thisdf2) +
    geom_col(aes(y=frac,x=cat),alpha=0.5) +
    labs(title=var) +
    scale_y_continuous(limits=c(0,1))
  print(a)
}

# Set results directory -------------------
KJL<-T
cluster<-F
if(KJL){
  setwd('/Users/kevinloope/Dropbox/GopherTortoise/Analysis/PVA_DT/Results/DT_03_31_2023A')
} 
if(cluster){
  setwd('/projects/birdnet/PVA/results/01_19_2023')
}

setwd('C:\\Users\\Kevin\\Dropbox\\GopherTortoise\\Analysis\\PVA_DT\\Results\\DT_03_21_2023_results_summary')

# Load data ----------------------

filetag <- "allreps_results"
allfilenames <- list.files()[grepl(filetag,list.files())]

datalist <- lapply(allfilenames,function(t) read.csv(t,row.names=1)  )

mydf <- do.call(rbind,datalist)

# Prep data ---------------------

# names(mydf)
nsims <- nrow(mydf)  # 3072 model runs

   # note: we never integrated climate scenarios into the results- need to re-integrate 
all_climscens <- c("rcp26","rcp45","rcp60","rcp85")
nclimscens <- length(all_climscens)

all_scens <- sort(unique(mydf$scenario))
nscenarios <- length(all_scens)   # 64

all_vitalrates <- c("GRa","GRk","PHI_A","PHI_J","PHI_H","MA","CS","PR","PF","HS","NS","dispersal")
nvitalrates <- length(all_vitalrates)

names(all_vitalrates) <-  c("Asympt. body size","Somatic growth rate","Adult Survival", "Juvenile Survival",
                      "Hatchling Survival", "Age at Maturity",
                      "Clutch Size","Reproductive fraction",
                      "Sex Ratio (frac female)","Hatching Success",
                      "Nest Success","Dispersal rate")

mydf$climscen <- rep(all_climscens,times=nscenarios)

rm(filetag,allfilenames,datalist)

# make sure that names of vital rates always match

colnames <- names(mydf)
colnames <- gsub("phi.a_pr","phi.a_clim",colnames)
colnames <- gsub("phi.a_","PHI_A_",colnames)
#colnames <- gsub("phi.j_","PHI_J_",colnames)
colnames <- gsub("cs_","CS_",colnames)
colnames <- gsub("pr_","PR_",colnames)
colnames <- gsub("hs_","HS_",colnames)
colnames <- gsub("pf_","PF_",colnames)
colnames <- gsub("ma_","MA_",colnames)
colnames <- gsub("growth_","GR_",colnames)
colnames <- gsub("SAM","GR_SAM",colnames)
colnames <- gsub("nestT","HS_nestT",colnames)
colnames

names(mydf) <- colnames

# NAMING VARIABLES ----------------

names(mydf)           # scenario (T/F) variables 
scen.vars=c(   #"climscen",        # climate scenario
               "GR_clim",         # is growth dependent on climate       
               "PHI_A_clim",      # is adult survival dependent on burn freq 
                #"PHI_J_burn",      # is juv survival dependent on burn freq
               "CS_clim",          # is clutch size dependent on latitude?
               "PR_clim",          # is reproduction prob dependent on latitude?
               "HS_clim",         # is hatching success dependent on climate?   
               "PF_clim"         # is sex ratio dependent on climate?
)
names(mydf) 
key_params <- c(
               "GRk_int",         # somatic growth rate, intercept
               "GRk_siteprec",    # growth rate, site effect
               "GRk_tmax",        # growth rate, temp effect
               "GRk_p",           # growth rate, precip effect
               "GRk_tmin",        # growth rate, tmin effect
               "GRa_int",         # asymptotic body size, intercept
               "GRa_siteprec",    # body size, site effect
               "GRa_tmax",        # body size, temp effect
               "GRa_p",           # body size, precip effect
               "GRa_tmin",        # asymptotic body size, tmin effect
               #"GRt0",            # body size at age 0 (used in VB growth model)
               "GR_SAM",          # size at maturity
               "PHI_A_base",      # baseline adult survival
               "PHI_A_p",         # effect of precip on adult survival
               "PHI_A_syeff",     # magnitude of site/year random effect on adult survival
               "PHI_H_base",      # baseline hatchling survival
               "PHI_J_base",      # baseline juv survival
               "CS_base",         # baseline clutch size
               "CS_t",            # effect of temperature on clutch size
               "CS_bs",           # effect of body size on clutch size
               "CS_p",            # effect of precip on clutch size
               "HS_base",         # baseline hatching success
               "HS_t",            # effect of temperature on hatching success
               #"HS_nestT",        # effect of air temperature on nest temp
               "NS_base",         # baseline nest success
               "PF_t",            # effect of temperature on sex ratio
               "PR_base",         # baseline repro fraction
               "PR_t",             # effect of temp on prob repro
               "PR_p",             # effect of precip on prob repro
               "PR_bs",            # effect of body size on prob repro
               "dispersal"       # dispersal rate
)

# note: ensure climate effects are 0 when not in effect
key_params
mydf$GRk_p <- mydf$GRk_p*mydf$GR_clim
mydf$GRk_tmax <- mydf$GRk_tmax*mydf$GR_clim
mydf$GRk_tmin <-  mydf$GRk_tmin*mydf$GR_clim
mydf$GRa_p <- mydf$GRa_p*mydf$GR_clim
mydf$GRa_tmax <- mydf$GRa_tmax*mydf$GR_clim
mydf$GRa_tmin <-  mydf$GRa_tmin*mydf$GR_clim
mydf$PHI_A_p <- mydf$PHI_A_p*mydf$PHI_A_clim
mydf$CS_t <- mydf$CS_t *mydf$CS_clim
mydf$CS_p <- mydf$CS_p*mydf$CS_clim
mydf$PR_t <- mydf$PR_t *mydf$PR_clim
mydf$PR_p <- mydf$PR_p*mydf$PR_clim
mydf$HS_t <- mydf$HS_t*mydf$HS_clim  # account for nest to air T relationship?
mydf$PF_t <- mydf$PF_t*mydf$PF_clim


# ancillary variables- need to keep, but not used directly
ancil.vars <- c(
  "scenario",
  "rep",
  "climscen"
)

# Define response variables

names(mydf) 
respvars <- c(    # ignore elasticity for now
  "lambda",          # lambda for this scenario at baseline climate conditions
  "final_abund",     # final abundance
  "delta_abund",     # change in abundance across simulation
  "frac_change",     # final abund as fraction of initial abund
  "final_adult",     # final adult abundance
  "final_hatch",     # final hatchling abundance
  "final_agerat"    # age ratio at end of simulation
)


# assess params for lambda >= 1 -----------------------

# Determine the sets of vital rates for which lambda was greater than 1- compare with 
   # tested ranges of parameters. 

sustainable <- subset(mydf,final_abund>=75000)
nrow(sustainable)   # 5679 runs with final abund >= 50000

nrow(sustainable)/nrow(mydf)  #58% of simulations

# read in original ranges for all variables!

# setwd("C:/Users/Kevin/Dropbox/GopherTortoise/Analysis/PVA/Rcode")
# source('GT_functions_Nov2022_v3.R')   # load functions from file
# datalist <- load_data()    # load all data needed for running models
# generators <- datalist$paraminfo
# rm(datalist)
# generators$gamma$dispersal$gen()

key_params

sust_list <- list()
p=1
for(p in 1:length(key_params)){
  thisparam <- key_params[p]
  sust_list[[thisparam]] <- list()
  sust_list[[thisparam]]$orig <- mydf[[thisparam]]  # replicate(1000,generators$PHI_A$adultsurv.mn$gen())
  sust_list[[thisparam]]$sust <- sustainable[[thisparam]]
}

# sust_list$CS_base <- list()
# csfunc <- function() max(0.01,thispar$CS$cs.intercept.mn$this + thispar$CS$cs.tempeff.mn$this * base_climvars["baseT"] + thispar$CS$CS_intercept$this )
# sust_list$CS_base$orig <- replicate(1000,generators$CS$CS_intercept$post())

v=1
for(v in 1:length(key_params)){
  orig_vs_sustainable(sust_list,key_params[v])
}


scen.vars

sust_list2 <- list()
p=2
for(p in 1:length(scen.vars)){
  thisparam <- scen.vars[p]
  sust_list2[[thisparam]] <- list()
  sust_list2[[thisparam]]$orig <- mydf[[thisparam]]  # replicate(1000,generators$PHI_A$adultsurv.mn$gen())
  sust_list2[[thisparam]]$sust <- sustainable[[thisparam]]
}

v=1
for(v in 1:length(scen.vars)){
  orig_vs_sustainable2(sust_list2,scen.vars[v])
}

# notes:

# look at climate scenarios

orig <- mydf$climscen
sust <- sustainable$climscen

thisdf <- data.frame(
  data = c(orig,sust),
  cat = c(rep("orig",times=length(orig)),rep("sust",times=length(sust)))
)
thisdf2 <- thisdf %>% 
  group_by(cat) %>% 
  summarize(rcp26 = mean(data=="rcp26"),
            rcp45 = mean(data=="rcp45"),
            rcp60 = mean(data=="rcp60"),
            rcp85 = mean(data=="rcp85")
  )

thisdf2
library(tidyr)
thisdf3 <- thisdf2 %>% pivot_longer(cols=!cat,names_to="climscen",values_to="frac") %>% filter(cat=="sust")

svg("plots/climscen_bars.svg",4,4)
a<- ggplot(thisdf3) +
  geom_col(aes(y=frac,x=climscen),alpha=0.5) +
  scale_y_continuous(limits=c(0,.5)) +
  labs(x="Climate Scenario",y="Fraction sustaintable")
print(a)
dev.off()




# Random forest and CART --------------------------

#### Define our formula (response ~ predictors)

response <- respvars[2]  # final abundance may be the best response var...

formula_scen <- as.formula(paste(response,"~",paste(scen.vars,collapse="+")))
formula_keyvars <- as.formula(paste(response,"~",paste(key_params,collapse="+")))
allvars <- c(scen.vars,key_params,"climscen")
formula_allvars <- as.formula(paste(response,"~",paste(allvars,collapse="+")))

# Conditional inference trees -------------------------

all_climscens  # pick a climate scenario...  
mydf_byclim <- list()
temp <- lapply(all_climscens,function(t) mydf_byclim[[t]] <<- subset(mydf,climscen==t))

scen_tree <- ctree(formula=formula_scen, data=mydf_byclim$rcp60, controls = ctree_control(mincriterion = 0.85,
                                                                      maxdepth = 5))
# graphics.off()

svg("plots/scen_tree_rerun.svg",9,5)
plot(scen_tree)
dev.off()


keyvars_tree <- ctree(formula=formula_keyvars, data=mydf_byclim$rcp60, controls = ctree_control(mincriterion = 0.85,
                                                                                          maxdepth = 5))
pdf("plots/keyvars_tree_rerun.pdf",40,10)
plot(keyvars_tree)
dev.off()

mydf$climscen <- as.factor(mydf$climscen)
allvars_tree <- ctree(formula=formula_allvars, data=mydf, controls = ctree_control(mincriterion = 0.85,
                                                                                                maxdepth = 5))
pdf("plots/allvars_tree_rerun.pdf",40,10)
plot(allvars_tree)
dev.off()


# RFE -------------------------------

# KTS: doesn't seem to have worked?

# Add random variables 
library(faux)

response

# mydf
mydf <- mydf %>% 
  mutate(
    rand1 = as.factor(sample(sample(letters[1:3], nrow(mydf), replace = TRUE))),
    rand2 = rnorm(nrow(mydf), mean = 10, sd = 2),
    corr1 = rnorm_pre(mydf[,response], mu = 10, sd = 2, r = 0.2, empirical = TRUE),
    corr2 = rnorm_pre(mydf[,response], mu = 10, sd = 2, r = 0.5, empirical = TRUE)
  )

library(caret)
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "cv", # repeated cv
                      # repeats = 5, # number of repeats
                      number = 5) # number of folds


# Features
x <- mydf %>%
  select(all_of(setdiff(allvars,c(response,"rand1","rand2","corr1","corr2")))) %>%
  as.data.frame()

# Target variable
y <- mydf[[response]]

# run the RFE algorithm- takes a VERY LONG time- four days at least. 
# result_rfe1 <- rfe(x = x, 
#                    y = y, 
#                     sizes = c(1:15),
#                     rfeControl = control)

# Print the results
# result_rfe1
# 

# plot(result_rfe1)   # KTS: this looks very strange..

# Print the selected features
# predictors(result_rfe1)

# result_rfe1
# save.image(file='rfe_environment.RData')

# ggplot(data = result_rfe1, metric = "Rsquared") + theme_bw()
# no2t <- ifelse(no2,"NO2","noNO2")
# svg(sprintf("plots/rfe_%s_%s.svg",response,no2t),3,4)
# thisplot <- ggplot(data = result_rfe1, metric = "RMSE") + theme_bw() + 
#   scale_x_continuous(limits=c(0,10),breaks = c(0, 2, 4, 6, 8, 10) )
# print(thisplot)
# dev.off()


# # select top variables
# set.size = 5 
# lm.vars <- result_rfe1$variables 
# lm.set <- lm.vars[lm.vars$Variables==set.size,  ] 
# lm.set <- aggregate(lm.set[, c("Overall")], list(lm.set$var), mean)
# lm.order <- order(lm.set[, c("x")], decreasing = TRUE)[1:set.size]
# tops <- lm.set[lm.order, ][["Group.1"]]
# 
# save.image(file='rfe_environment.RData')

# Random forest analysis --------------------

#  first make the formulas
scen.vars_plus <- c(scen.vars,c("rand1","rand2","corr1","corr2"))
formula_scen <- as.formula(paste(response, "~", paste(scen.vars_plus,collapse="+")))

key_params_plus <- c(key_params,c("rand1","rand2","corr1","corr2"))
formula_keyvars <- as.formula(paste(response, "~", paste(key_params_plus,collapse="+")))

allvars_plus <- c(allvars,c("rand1","rand2","corr1","corr2"))
formula_allvars <- as.formula(paste(response, "~", paste(allvars_plus,collapse="+")))
formula_allvars2 <- as.formula(paste(response, "~", paste(allvars,collapse="+")))

# fit model with just the scenario variables
rfmod_scen <- ranger(formula_scen,data=mydf,importance = 'permutation')

varimp <- importance(rfmod_scen)

options(scipen=20)

varimp <- varimp[order(varimp,decreasing=T)]

svg("plots/varimp_scen.svg",4,6)
  par(mai=c(1,2,0.1,0.1))
  par(las=2)
  plot(length(varimp):1~varimp,yaxt="n",xlab="importance",ylab="",
       xlim=c(0,max(varimp)+0.1),xaxt="n")
  axis(2,at=length(varimp):1,labels = names(varimp)  )
dev.off()

# fit model with the quantitative variables
rfmod_keyvars <- ranger(formula_keyvars,data=mydf,importance = 'permutation')

varimp <- importance(rfmod_keyvars)

options(scipen=20)

varimp <- varimp[order(varimp,decreasing=T)]

svg("plots/varimp_keyvars.svg",4,7)
  par(mai=c(1,2,0.1,0.1))
  par(las=2)
  plot(length(varimp):1~varimp,yaxt="n",xlab="importance",ylab="",
       xlim=c(0,max(varimp)+0.1),xaxt="n")
  axis(2,at=length(varimp):1,labels = names(varimp)  )
dev.off()

# fit rf model with all variables
rfmod_allvars <- ranger(formula_allvars2,data=mydf,importance = 'permutation')

varimp <- importance(rfmod_allvars)

options(scipen=20)

varimp <- varimp[order(varimp,decreasing=T)]

svg("plots/varimp_allvars.svg",4,6)
  par(mai=c(1,2,0.1,0.1))
  par(las=2)
  plot(20:1~varimp[1:20],yaxt="n",xlab="importance",ylab="",
       xlim=c(0,max(varimp)+0.1),xaxt="n")
  axis(2,at=20:1,labels = names(varimp[1:20])  )
dev.off()

varimp_allvars <- varimp 


# RF effects plots

# bestvars <- predictors(result_rfe1)

toplot2 <- names(varimp)[!names(varimp)%in%c("rand1","rand2","corr1","corr2")]

toplot <- toplot2[1:10]

p=1
for(p in 1:length(toplot)){
  thisvar <- toplot[p]
  svg(sprintf("plots/%s_effectsRF_%s.svg",p,thisvar),4,4)
  # par(mfrow=c(2,2))
  if(is.numeric(mydf[[thisvar]])){
    nd <- data.frame(x=seq(min(mydf[[thisvar]]),max(mydf[[thisvar]]),length=50))
  }else{
    nd <- data.frame(x = levels(as.factor(mydf[[thisvar]])))
  }
  
  names(nd) <- thisvar
  
  if(class(mydf[[thisvar]])=="logical") nd[[thisvar]] <- as.logical(nd[[thisvar]])
  
  othervars <- setdiff(allvars_plus,thisvar)
  
  othervars_numeric <- othervars[sapply(mydf[,othervars],class  ) %in% c("numeric","integer")]
  othervars_logical <- othervars[sapply(mydf[,othervars],class  ) %in% c("logical")]
  othervars_cat <- othervars[sapply(mydf[,othervars],class  ) %in% c("factor","character")]
  
  temp <- sapply(othervars_numeric,function(t) nd[[t]] <<- mean(mydf[[t]])  )
  temp <- sapply(othervars_logical,function(t) nd[[t]] <<- FALSE  )
  temp <- sapply(othervars_cat,function(t) nd[[t]] <<- factor(names(which.max(table(mydf[[t]]))),levels=levels(as.factor(mydf[[t]]))) )
  
  pred = predict(rfmod_allvars,data=nd,type="response")$predictions
  
  if(is.numeric(mydf[[thisvar]])){
    par(mai=c(1,1,.1,.1))
    plot(pred/1000~nd[,1],type="l",xlab=thisvar,main="",ylab=response,xaxt="n")
    rug(mydf[[thisvar]][sample(1:nrow(mydf),100,replace=T)])
    axis(1,at=seq(min(mydf[[thisvar]]),max(mydf[[thisvar]]),length=5),
         labels=sigfig(seq(min(mydf[[thisvar]]),max(mydf[[thisvar]]),length=5)) )
  }else{
    par(mai=c(1,1,.5,.1))
    plot(pred/1000~as.factor(nd[,1]),ylim=c(0,max(pred/1000)+2),
         xlab="",main=thisvar,ylab=response)
  }
  dev.off()
}



# interactions -----------------------

# just test interactions among top vars?

# select top variables
set.size = 15
tops <- toplot2[1:set.size]

intdf <- expand.grid(tops,tops,stringsAsFactors = F)
intdf <- intdf[intdf$Var1!=intdf$Var2,]
intdf2 <- NULL
i=1
for(i in 1:nrow(intdf)){
  if(!paste0(intdf$Var2[i],intdf$Var1[i])%in%paste0(intdf2$Var1,intdf2$Var2)){
    intdf2 <- rbind(intdf2,intdf[i,])
  }
}
intdf2$IntStrength <- NA
intdf2$IntStrength2 <- NA

i=1
for(i in 1:nrow(intdf2)){
  v1 <- intdf2$Var1[i]
  v2 <- intdf2$Var2[i]
  if(class(mydf[[v1]])%in%c("numeric","integer")){
    x1 <- as.numeric(quantile(mydf[[v1]],seq(0,1,length=10)))
  }else if(class(mydf[[v1]])%in%c("logical")){
    x1 <- c(T,F)
  }else if(class(mydf[[v1]])%in%c("factor")){
    x1 <- factor(levels(mydf[[v1]]),levels=levels(mydf[[v1]]))
  } 
  
  if(class(mydf[[v2]])%in%c("numeric","integer")){
    x2 <- as.numeric(quantile(mydf[[v2]],seq(0,1,length=10)))
  }else if(class(mydf[[v2]])%in%c("logical")){
    x2 <- c(T,F)
  }else if(class(mydf[[v2]])%in%c("factor")){
    x2 <- factor(levels(mydf[[v2]]),levels=levels(mydf[[v2]]))
  } 
  nd <- expand.grid(x1,x2)
  names(nd) <- c(v1,v2)
  
  othervars <- setdiff(allvars_plus,c(v1,v2))
  
  othervars_numeric <- othervars[sapply(mydf[,othervars],class  ) %in% c("numeric","integer")]
  othervars_logical <- othervars[sapply(mydf[,othervars],class  ) %in% c("logical")]
  othervars_cat <- othervars[sapply(mydf[,othervars],class  ) %in% c("factor","character")]
  
  temp <- sapply(othervars_numeric,function(t) nd[[t]] <<- mean(mydf[[t]])  )
  temp <- sapply(othervars_logical,function(t) nd[[t]] <<- FALSE  )
  temp <- sapply(othervars_cat,function(t) nd[[t]] <<- factor(names(which.max(table(mydf[[t]]))),levels=levels(as.factor(mydf[[t]]))) )
  
  nd$pred = predict(rfmod_allvars,data=nd,type="response")$predictions
  tempform <- as.formula(sprintf("pred ~ as.factor(%s)+as.factor(%s)",v1,v2))
  addmod <- lm(tempform,data=nd)
  nd$pred2 = predict(addmod)
  intdf2$IntStrength[i] <- sqrt(mean((nd$pred-nd$pred2)^2))
  intdf2$IntStrength2[i] <- intdf2$IntStrength[i] / sd(nd$pred)
}

topint <- intdf2[which.max(intdf2$IntStrength2),1:2]  
# cor(mydf2[,toplot])  # double check if top variables are highly correlated



# visualize top interactions --------------

   # note: not yet set up for non-quantitative variables.

intdf2 <- intdf2[order(intdf2$IntStrength2,decreasing=T),]

toplot <- intdf2[1:5,]

p=1
for(p in 1:nrow(toplot)){
  thisint <- unlist(toplot[p,1:2])
  v1 <- thisint[1]
  v2 <- thisint[2]
  svg(sprintf("plots/%s_interaction_%s_%s.svg",p,v1,v2),5,3)
  if(class(mydf[[v1]])%in%c("numeric","integer")){
    x1 <- seq(min(mydf[[v1]]),max(mydf[[v1]]),length=15)
  }else if(class(mydf[[v1]])%in%c("logical")){
    x1 <- c(T,F)
  }else if(class(mydf[[v1]])%in%c("factor")){
    x1 <- factor(levels(mydf[[v1]]),levels=levels(mydf[[v1]]))
  } 
  if(class(mydf[[v2]])%in%c("numeric","integer")){
    x2 <- seq(min(mydf[[v2]]),max(mydf[[v2]]),length=15)
  }else if(class(mydf[[v2]])%in%c("logical")){
    x2 <- c(T,F)
  }else if(class(mydf[[v2]])%in%c("factor")){
    x2 <- factor(levels(mydf[[v2]]),levels=levels(mydf[[v2]]))
  } 
  nd <- expand.grid(x1,x2)
  names(nd) <- c(v1,v2)
  
  othervars <- setdiff(allvars_plus,c(v1,v2))
  
  othervars_numeric <- othervars[sapply(mydf[,othervars],class  ) %in% c("numeric","integer")]
  othervars_logical <- othervars[sapply(mydf[,othervars],class  ) %in% c("logical")]
  othervars_cat <- othervars[sapply(mydf[,othervars],class  ) %in% c("factor","character")]
  
  temp <- sapply(othervars_numeric,function(t) nd[[t]] <<- mean(mydf[[t]])  )
  temp <- sapply(othervars_logical,function(t) nd[[t]] <<- FALSE  )
  temp <- sapply(othervars_cat,function(t) nd[[t]] <<- factor(names(which.max(table(mydf[[t]]))),levels=levels(as.factor(mydf[[t]]))) )
  
  nd[[response]] = predict(rfmod_allvars,data=nd,type="response")$predictions

  # nd
  if(all(c(class(mydf[[v1]]),class(mydf[[v2]]))%in%c("numeric","integer") ) ){
    # Heatmap
    thisplot <- ggplot(nd, aes(.data[[v1]], .data[[v2]], fill= .data[[response]])) +
      geom_tile() +
      scale_fill_gradient(low="white", high="blue")   # +
    # scale_x_continuous( #limits=c(max(0,min(x)),max(x)),
    #   breaks = c(min(x),seq(min(x),max(x),length=3),max(x)),
    #   labels = sigfig(pmax(0,c(min(x),seq(min(x),max(x),length=3),max(x)) ),2) ) +
    # scale_y_continuous( #limits=c(max(0,min(y)),max(y)),
    #   breaks = c(min(y),seq(min(y),max(y),length=3),max(y)),
    #   labels = sigfig(pmax(0,c(min(y),seq(min(y),max(y),length=3),max(y)) ),2) )
    print(thisplot)
  }else{
    which.fac <- which(c(class(mydf[[v1]]),class(mydf[[v2]]))%in%c(c("logical","factor","character")) )
    facvar <- c(v1,v2)[which.fac]
    which.nonfac <- setdiff(c(1,2),which.fac)
    nonfac <- c(v1,v2)[which.nonfac]
    thisplot <- ggplot(nd,aes(.data[[nonfac]], .data[[response]])) +
      geom_point(aes(col=.data[[facvar]])) +
      geom_smooth(aes(col=.data[[facvar]]))
    # plot(nd[[v1]])
    print(thisplot)
  }   
  dev.off()

}


write.csv(intdf2,file="plots/interaction_table.csv",row.names=F)

