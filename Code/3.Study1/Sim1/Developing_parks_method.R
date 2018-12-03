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
test.length <- 24
n.stage <- 3
pform <- "133"
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

##----------------------------------------------------------------------------
# Analytical evaluation
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
bias.anal <- cond_moments[1, ] - theta
csem.anal <- sqrt(cond_moments[2, ])

##----------------------------------------------------------------------------
# Analytical evaluation using Park et al. (2017)
# find points where adjacent two primary routes intercet  
if(n.stage == 2) {
  x <- eos_list$df_path$stage2[c(1, 4, 7)]
} 
if(n.stage == 3){
  x <- eos_list$df_path$stage3[c(1, 2, 3)]
}
cuts <- cross_info(x, RDP=c(-.44, .44), range.theta=range.theta, D, interval=c(-2, 2))$cut.score

# compute MST test-level information using the simple moving average method
info_MST <- purrr::map_dbl(theta, mov_ave, x=x, cut.score=cuts, D=D)

# compute the conditional SEEs  
csem.anal2 <- 1/sqrt(info_MST)

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
bias.simEOS <- map(myeval_EOS, function(x) x$cond.eval[, 1])
csem.simEOS <- map(myeval_EOS, function(x) sqrt(x$cond.eval[, 2]))

##----------------------------------------------------------------------------
# Simulation 2: MLE scoring
module.mat <- mstBU$itemvar.df[[1]]
route.scoring="EAP"
final.scoring="ML"
range=c(-5, 5)
source.files <- file.path(getwd(), "R", src.files)

nrep <- 1000

# implement runMST_MLErep
simMLE <- map(nrep, function(x) runMST_MLErep(x, theta, item.pool, module.mat, ini.module=1, route.map, D, route.scoring=route.scoring, 
                                              final.scoring=final.scoring, range=range, full.out=TRUE, iseed=0L))

# evaluation
myeval_MLE <- map(simMLE, evalMST, weights=w)
bias.simMLE <- map(myeval_MLE, function(x) x$cond.eval[, 1])
csem.simMLE <- map(myeval_MLE, function(x) sqrt(x$cond.eval[, 2]))

##----------------------------------------------------------------------------
# Classification accuracy
quadrature <- list(theta, w)
cacIRT::Rud.D(cutscore=1.0, quadrature, sem=csem.anal)

cacIRT::Rud.D(cutscore=1.0, quadrature, sem=csem.simEOS[[1]])

cacIRT::Rud.D(cutscore=1.0, quadrature, sem=csem.simMLE[[1]])

##----------------------------------------------------------------------------
# plot bias
plot(bias.simMLE[[1]] ~ theta, type='l', col=1, lwd=2, 
     ylab="BIAS",
     xlab="Proficiency")
lines(bias.anal ~ theta, col=2, lwd=2)
lines(bias.simEOS[[1]] ~ theta, col=3, lwd=3)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)

# plot csee
plot(csem.simMLE[[1]] ~ theta, type='l', col=1, lwd=2, 
     ylab="CSEM",
     xlab="Proficiency")
lines(csem.anal ~ theta, col=2, lwd=2)
lines(csem.simEOS[[1]] ~ theta, col=3, lwd=3)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)


