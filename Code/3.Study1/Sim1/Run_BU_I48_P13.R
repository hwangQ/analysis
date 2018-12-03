###############################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 11/26/2018
#Title: Bottom-up assembly code 
###############################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(cacIRT)
library(mstR)
library(tidyverse)

##----------------------------------------------------------------------------
# Read source code
src.files <- list.files("R/", pattern="*.R")
purrr::walk(src.files, function(x) source(file.path("R", x)))

##----------------------------------------------------------------------------
# set conditions
test.length <- 48
n.stage <- 2
pform <- "13"
theta <- seq(-4, 4, .1)
RDP <- list(NA, c(-0.44, 0.44), c(-0.44, 0.44))
size.module <- rep(test.length/n.stage, n.stage)
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.1, mip.gap=c(0.1, 0.05))
route.map <- read.csv(paste0("Input/Routing_map_P", pform, ".csv"), header=FALSE) # read a routing map
cont_require <- read.csv(paste0("Input/Content_BU_I", test.length, "_P", pform,".csv"))
constraints <- list(route.map=route.map, size.module=size.module, content=cont_require[, 2:(1 +n.stage)])

##----------------------------------------------------------------------------
# read item bank and content category information
bank_info <- readRDS("Input/item_pool_1.rds")
bank_info <- data.frame(ID=1:nrow(bank_info), bank_info)
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS")

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- select(bank_info, ID, CATEGORY, MODEL, CLASS, starts_with("PARAM"))
item.pool <- bank_info

##----------------------------------------------------------------------------
# create a list of target module information functions (MIFs)
theta.stg1 <- list(c(-2, median(c(-2.0, -0.44)), median(c(-0.44, 0.44)), median(c(0.44, 2.0)), 2))
theta.stg2 <- list(seq(-2.0, -0.44, length.out=5), 
                   seq(-0.44, 0.44, length.out=5),
                   seq(0.44, 2.0, length.out=5))
theta.stg3 <- theta.stg2
target.theta <- list(stage1=theta.stg1, stage2=theta.stg2, stage3=theta.stg3)
targetMIF <- info_mif(item.pool, target.theta, route.map, test.length, D)

##----------------------------------------------------------------------------
# Automated test assembly using the bottom-up approach 
mstBU <- ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), 
                   D=1.702, divide.D=FALSE, lp.control=lp.control)

# save results
dir_out <- file.path("Output/Study1/Sim1", paste0("I", test.length, "_", pform, "MST"))
saveRDS(mstBU, file=file.path(dir_out, "mstBU.rds"))

# read results
# mstBU <- readRDS(file.path(dir_out, "mstBU.rds"))

##----------------------------------------------------------------------------
# Analytical evaluation 1
# read the assembled test forms
ata_forms <- mstBU$prm.df

# set condtions
range.theta <- c(-5, 5)
interval <- c(-2.5, 2.5)

# read a routing map
panel.info <- panel_info(route.map)
config.info <- panel.info$config
n.module <- panel.info$n.module
pathway <- panel.info$pathway

# estimate the observed equated scores across all (sub) pathways
module.form <- ata_forms
eos_list <- est_eos(module.form, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)

# find points where two adjacent TIFs intersect across all stages
cut.score <- cutoff(module.form, route.map, RDP, D, range.theta, interval)

# conditional raw score distribution of each module given an ability value
cond_dist <- ability_dist_sep(module.form, theta, route.map, D)

# conditional ability distribution of total-test ability estimates given a true ability
x <- cond_dist
eos.path <- eos_list$eos_path
n.path <- eos_list$n_path
joint.dist <- ability_dist(x, eos.path, n.path, cut.score)

# calculate a mean and sd of conditional ability distribution at each ability value
if(n.stage == 2) {
  cond_moments <- sapply(joint.dist$joint.dist$stage2, cal_moments, node=eos_list$eos_path$stage2)
}
if(n.stage == 3) {
  cond_moments <- sapply(joint.dist$joint.dist$stage3, cal_moments, node=eos_list$eos_path$stage3)
}
bias_anal1 <- cond_moments[1, ] - theta
csee_anal1 <- sqrt(cond_moments[2, ])

# save results
saveRDS(bias_anal1, file=file.path(dir_out, "bias_anal1.rds"))
saveRDS(csee_anal1, file=file.path(dir_out, "csee_anal1.rds"))

# read results
# bias_anal1 <- readRDS(file.path(dir_out, "bias_anal1.rds"))
# csee_anal1 <- readRDS(file.path(dir_out, "csee_anal1.rds"))

##----------------------------------------------------------------------------
# Analytical evaluation 2: using Park et al. (2017)
# find points where adjacent two primary routes intercet  
if(n.stage == 2) {
  x <- eos_list$df_path$stage2[c(1, 2, 3)]
} 
if(n.stage == 3){
  x <- eos_list$df_path$stage3[c(1, 4, 7)]
}
cuts <- cross_info(x, RDP=c(-.44, .44), range.theta=range.theta, D, interval=c(-2, 2))$cut.score

# compute MST test-level information using the simple moving average method
info_MST <- purrr::map_dbl(theta, mov_ave, x=x, cut.score=cuts, D=D)

# compute the conditional SEEs  
csee_anal2 <- 1/sqrt(info_MST)

# save results
saveRDS(csee_anal2, file=file.path(dir_out, "csee_anal2.rds"))

# read results
# csee_anal2 <- readRDS(file.path(dir_out, "csee_anal2.rds"))

##----------------------------------------------------------------------------
# Simulation 1: Equated NC scoring method
source.files <- file.path(getwd(), "R", src.files)
nrep <- c(100, 1000, 5000)

# implement runMST_rep
simEOS <- map(nrep, function(x) runMST_rep(x, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE))

# evaluation
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval_EOS <- map(simEOS, evalMST, weights=w)

# save results
saveRDS(simEOS, file=file.path(dir_out, "simEOS.rds"))
saveRDS(myeval_EOS, file=file.path(dir_out, "myeval_EOS.rds"))

# read results
# simEOS <- readRDS(file.path(dir_out, "simEOS.rds"))
# myeval_EOS <- readRDS(file.path(dir_out, "myeval_EOS.rds"))

bias_simEOS <- map(myeval_EOS, function(x) x$cond.eval[, 1])
csee_simEOS <- map(myeval_EOS, function(x) sqrt(x$cond.eval[, 2]))

##----------------------------------------------------------------------------
# Simulation 2: MLE scoring
module.mat <- mstBU$itemvar.df[[1]]
route.scoring <- "EAP"
final.scoring <- "ML"
range <- range.theta 
cut.mat <- cutoff_mat(cut.score, n.stage, n.module)
source.files <- file.path(getwd(), "R", src.files)

# implement runMST_MLErep
simMLE <- map(nrep, function(x) runMST_MLErep(x, theta, item.pool, module.mat, ini.module=1, route.map, cut.mat, D, route.scoring=route.scoring, 
                                              final.scoring=final.scoring, range=range, full.out=TRUE, iseed=0L))

# evaluation
myeval_MLE <- map(simMLE, evalMST, weights=w)

# save results
saveRDS(simMLE, file=file.path(dir_out, "simMLE.rds"))
saveRDS(myeval_MLE, file=file.path(dir_out, "myeval_MLE.rds"))

# read results
# simMLE <- readRDS(file.path(dir_out, "simMLE.rds"))
# myeval_MLE <- readRDS(file.path(dir_out, "myeval_MLE.rds"))

bias_simMLE <- map(myeval_MLE, function(x) x$cond.eval[, 1])
csee_simMLE <- map(myeval_MLE, function(x) sqrt(x$cond.eval[, 2]))





