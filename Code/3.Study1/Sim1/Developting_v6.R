######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 11/05/2018
#Title: Bottom-up assembly code 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(cacIRT)
library(mstR)
library(tidyverse)

##########################################################################
# Read source code
src.files <- list.files("R/", pattern="*.R")
purrr::walk(src.files, function(x) source(file.path("R", x)))

# set directory
# mydir <- "E:/Dissertation/Analysis/Temp"
# setwd(mydir)

# set conditions
test.length <- 32

# read item bank and content category information
bank_info <- readRDS("Input/item_pool_1.rds")
bank_info <- data.frame(ID=1:nrow(bank_info), bank_info)
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS")
cont_require <- read.csv(paste0("Input/Content_BU_I", test.length, ".csv"))

# transform a content variable to numeric variable
# x <- bank_info[, 5]
# x <- as.numeric(x)
# bank_info$CLASS_RE <- x

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- select(bank_info, ID, CATEGORY, MODEL, CLASS, starts_with("PARAM"))

###########################################################################
# set condtions for ATA
theta <- seq(-8, 8, .1)
RDP <- list(NA, c(-0.44, 0.44), c(-0.44, 0.44))
item.pool <- bank_info
size.module <- c(8, 8, 8)
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))
route.map <- read.csv("Input/Routing_map_1-3.csv", header=FALSE) # read a routing map
constraints <- list(route.map=route.map, size.module=size.module, content=cont_require[, 2:4])

# create a list of target module information functions (MIFs)
theta.stg1 <- list(c(-2, median(c(-2.0, -0.44)), median(c(-0.44, 0.44)), median(c(0.44, 2.0)), 2))
theta.stg2 <- list(seq(-2.0, -0.44, length.out=5), 
                   seq(-0.44, 0.44, length.out=5),
                   seq(0.44, 2.0, length.out=5))
theta.stg3 <- theta.stg2
target.theta <- list(stage1=theta.stg1, stage2=theta.stg2, stage3=theta.stg3)
targetMIF <- info_mif(item.pool, target.theta, route.map, test.length)

# example 
mstBU <- ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
         lp.control=list(timeout=30, epsint=0.2, mip.gap=c(0.1, 0.05)))


################################################################################
# read the assembled test forms
ata_forms <- mstBU$prm.df

# set condtions
range.theta <- c(-6, 6)
interval <- c(-2.5, 2.5)

# read a routing map
panel.info <- panel_info(route.map)
config.info <- panel.info$config
n.module <- panel.info$n.module
pathway <- panel.info$pathway


# estimate the observed equated scores across all (sub) pathways
module.form <- ata_forms
eos_list <- est_eos(module.form, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)

eos_list$df_path

test_info <- test.info(eos_list$df_path$stage3[[7]], theta, D)$testInfo
plot(test_info ~ theta, type='l')




a <- eos_list$df_path$stage3[[1]]$PARAM.1
b <- eos_list$df_path$stage3[[1]]$PARAM.2
g <- eos_list$df_path$stage3[[1]]$PARAM.3
icc <- t(pl.fn(theta, a, b, g, D))
tcc <- colSums(icc)
plot(tcc ~ theta, type='l', col=1, lwd=2)

for(i in 2:7) {
a <- eos_list$df_path$stage3[[i]]$PARAM.1
b <- eos_list$df_path$stage3[[i]]$PARAM.2
g <- eos_list$df_path$stage3[[i]]$PARAM.3
icc <- t(pl.fn(theta, a, b, g, D))
tcc <- colSums(icc)
lines(tcc ~ theta, type='l', col=i, lwd=2)
}
legend("topleft", legend=paste0("Path.", 1:7), col=1:7, lty=1, lwd=2)

res <- simMLE$responses
sum.list <- lapply(res, rowSums)
sum.4 <-  unlist(lapply(sum.list, '[[', 81))

# find points where two adjacent TIFs intersect across all stages
cut.score <- cutoff(module.form, route.map, RDP, D, range.theta, interval)

################################################################################
# Simulation 
source.files <- paste0(src.dir, src.files)
nrep <- 300

# implement runMST_rep
simTCC <- runMST_rep(nrep, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE)

# evaluation
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval <- evalMST(simTCC, w)
bias.sim <- myeval$cond.eval[, 1]
csem.sim <- sqrt(myeval$cond.eval[, 2])

################################################################################
# Analytical evaluation
# conditional raw score distribution of each module given an ability value
cond_dist <- ability_dist_sep(module.form, theta, route.map, D)

# conditional ability distribution of total-test ability estimates given a true ability
x <- cond_dist
eos.path <- eos_list$eos_path
n.path <- eos_list$n_path
joint.dist <- ability_dist(x, eos.path, n.path, cut.score)

# calculate a mean and sd of conditional ability distribution at each ability value
cond_moments <- sapply(joint.dist$joint.dist$stage3, cal_moments, node=eos_list$eos_path$stage3)
bias.anal <- cond_moments[1, ] - theta
csem.anal <- sqrt(cond_moments[2, ])

# plot bias
plot(bias.sim ~ theta, type='l', col=1, lwd=2, 
     ylab="BIAS",
     xlab="Proficiency")
lines(bias.anal ~ theta, col=2, lwd=2)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)

# plot csem
plot(csem.anal ~ theta, type='l', col=1, lwd=2, 
     ylab="CSEM",
     xlab="Proficiency")
lines(csem.anal ~ theta, col=2, lwd=2)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)

################################################################################
# Classification accuracy

quadrature <- list(theta, w)
cacIRT::Rud.D(cutscore=1.0, quadrature, sem=csem.anal)


################################################################################
# Simulation using MLE scoring
# a function to administer MST for individual examinee
module.mat <- mstBU$itemvar.df[[1]]
route.scoring="ML"
final.scoring="ML"
range=c(-6, 6)
source.files <- paste0(src.dir, src.files)
theta <- seq(-4, 4, .1)
nrep <- 300

simMLE <- runMST_MLErep(nrep, theta, item.pool, module.mat, ini.module=1, route.map, D, route.scoring=route.scoring, 
                        final.scoring=final.scoring, range=range, full.out=TRUE, iseed=0L) 

# evaluation
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval_ML <- evalMST(simMLE, w)
bias.simML <- myeval_ML$cond.eval[, 1]
csem.simML <- sqrt(myeval_ML$cond.eval[, 2])
cacIRT::Rud.D(cutscore=1.0, quadrature, sem=csem.simML)

# plot bias
plot(bias.simML ~ theta, type='l', col=1, lwd=2, 
     ylab="BIAS",
     xlab="Proficiency")
lines(bias.anal ~ theta, col=2, lwd=2)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)

# plot
plot(csem.simML ~ theta, type='l', col=1, lwd=2, 
     ylab="CSEM",
     xlab="Proficiency")
lines(csem.anal ~ theta, col=2, lwd=2)
legend("topright", legend=c("Simulation", "Analytic"), col=1:2, lwd=2)

lines(csem.sim ~ theta, col=3, lwd=2)

