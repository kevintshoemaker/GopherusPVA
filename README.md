
## Population viability analysis for Gopherus tortoises (G. polyphemus and G. agassizii)
Kevin T. Shoemaker, Kevin J. Loope, Elizabeth A. Hunter

Updated: 6/28/2023

### Overview
This repository contains code to run simulation models projecting population demographic responses to climate change under a variety of climate change scenarios and vital rate sensitivities to climate.  Each species is modeled separately using a shared framework: populations are simulated on range-wide 5x5km grids, divided into 20 age classes with spatially variable, environmentally-determined age at maturity based on empirical data.  Models are initialized with a burn-in step that accounts for unstable starting age-distributions, and models are run through the year 2099.  Vital rates, including parameters representing climate sensitivity, are drawn from empirically-derived distributions, with unique replicate draws of all parameters simultaneously analyzed across all scenarios.  Desert tortoise simulations include 256 scenarios representing various "on-off" combinations of climate sensitivities, as well as 4 climate change projection models.  Gopher tortoise simulations include 3584 scenarios representing various combinations of climate sensitivities and 8 climate change projection models.  For more details, see final report for SERDP project RC18-1103.  

Models are implemented in Program R.  Each species has 4 associated R scripts, an optional bash script for running on a SLURM cluster, and one associated setup_files folder:

**RunSimulations.R** - Open this file to run simulations.  

Running this script locally will start a model run on your computer.  To run on a computing cluster, generate and submit associated submission files that call this script using the accompanying .sh file (see below).

Before starting a run, open the code and specify:

* paths for output files
* whether the model will be run locally or on a SLURM-based computing cluster
* the number of replicates per scenario
* the number of cores to parallelize across

This script outputs one multi-band tif file for each year of each replicate of each scenario, with layers corresponding to age classes 1-20.  It also outputs a folder of population array rds files that contain range-wide summaries of each replicate in tabular form.  

**Functions.R** - This script contains all functions necessary to run the simulations.  This file is sourced by the other scripts and does not need to be loaded manually.

**Make_Repdf.R** - This script is run after the simulations have completed, and processes the output poparray files (containing the age structured results), computing summary population statistics and combining them with drawn model parameters for each replicate.  The resulting dataframe can be used for random forest analyses.  This is automatically run for the DT model, and can be automatically run or manually run for the GT model.

**Analyze_Results.R** - This file runs random forest analyses on replicate-level outputs from the simulation, identifying important baseline and climate-sensitivity parameters across the entire range.  

**Plot_Trajectories.R** - This script will plot the total abundance through time for each replicate for a specified scenario.  

**Plot_Vital_Rates.R** - This script will plot the vital rates (Age at maturity, Reproductive output, Hatching success, Fraction female, Probability of reproduction, and Lambda) for a given scenario, replicate and pixel.

**Prep_Animations.R** - This script will prepare an animation of range-wide population change through time for a given scenario and rep.

**Plot_Lambda.R** - This script will generate a bivariate heat map showing the effects of the two largest climate variables on lambda.  


### Notes on run times: 

Given the large number of scenarios, the gopher tortoise model, as currently specified, takes a very long time to run.  Our final run with 120 reps per scenario took ~480 compute hours on 120 cores when run on the VT ARC computing cluster.  To allow runs to meet the queueing requirements of the cluster, the total number of scenarios x replicates was divided into 240 chunks that were submitted as separate jobs (corresponding to arrayIDs of 1-240 in the RunSimulations.R script), each of which was run on a 120-core node (parallelized with doParallel).  Each job took ~2 hours to complete.  To reduce run times, reduce the number of scenarios by editing lines 288-302 in Functions.R, or reduce the number of replicates by editing line 112 in RunSimulations.R.  

The desert tortoise model was similarly broken into 24 jobs and submitted to VT ARC's queuing system.  Each job took ~1.75 hours of compute time on a single node, for a total of 42 compute hours on 120 cores.    









