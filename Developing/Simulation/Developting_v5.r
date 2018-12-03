######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 09/28/2018
#Title: Bottom-up assembly code 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)

##########################################################################
# Read source code
src.dir <- "E:/Dissertation/Analysis/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/Dissertation/Analysis/Temp"
setwd(mydir)

# read item bank and content category information
bank_info <- read.csv("itempar_bank_2PL.csv")
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS")
cont_require <- read.csv("Content_Requirement.csv")

# transform a content variable to numeric variable
x <- bank_info[, 5]
x <- as.numeric(x)
bank_info$CLASS_RE <- x

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- bank_info[, c(1, 7, 8, 5, 6, 2, 3, 4)]

###########################################################################
# set condtions for ATA
theta <- seq(-4, 4, .1)
RDP <- list(NA, c(-0.44, 0.44), c(-0.44, 0.44))
item.pool <- bank_info
test.length <- 45
size.module <- c(15, 15, 15)
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))
route.map <- read.csv("Routing_map.csv", header=FALSE) # read a routing map
constraints <- list(route.map=route.map, size.module=size.module, content=cont_require[, 2:4])

# create a list of target module information functions (MIFs)
theta.stg1 <- list(c(median(c(-2.0, -0.44)), median(c(-0.44, 0.44)), median(c(0.44, 2.0))))
theta.stg2 <- list(seq(-2.0, -0.44, length.out=5), 
                   seq(-0.44, 0.44, length.out=5),
                   seq(0.44, 2.0, length.out=5))
theta.stg3 <- theta.stg2
target.theta <- list(stage1=theta.stg1, stage2=theta.stg2, stage3=theta.stg3)
targetMIF <- info_mif(target.theta, route.map)

# example 
mstBU <- ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
         lp.control=list(timeout=30, epsint=0.2, mip.gap=c(0.1, 0.05)))


################################################################################
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

################################################################################
# Simulation 
source.files <- paste0(src.dir, src.files)
nrep <- 1000

# implement runMST_rep
x <- runMST_rep(nrep, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE)

# evaluation
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval <- evalMST(x, w)
bias.sim <- myeval$cond.eval[, 1]
csem.sim <- sqrt(myeval$cond.eval[, 2])

################################################################################
# Analytical evaluation
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
plot(csem.sim ~ theta, type='l', col=1, lwd=2, 
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
route.scoring="WL"
final.scoring="ML"
range=c(-5, 5)
source.files <- paste0(src.dir, src.files)
theta <- seq(-4, 4, .1)
nrep <- 1000

simMLE <- runMST_MLErep(nrep, theta, item.pool, module.mat, ini.module=1, route.map, D, route.scoring=route.scoring, 
                        final.scoring=final.scoring, range=range, full.out=FALSE, iseed=0L) 

# evaluation
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval_ML <- evalMST(simMLE, w)
bias.simML <- myeval_ML$cond.eval[, 1]
csem.simML <- sqrt(myeval_ML$cond.eval[, 2])

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

