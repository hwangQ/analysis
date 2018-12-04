###############################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 12/02/2018
#Title: Bottom-up assembly code (Simulation 2)
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
n.stage <- 2
pform <- "13"
cut_thread <- c(0, 0.524)
theta <- seq(-4, 4, .1)
theta_1000 <- readRDS(file.path("Input", "theta_1000.rds"))
theta_5000 <- readRDS(file.path("Input", "theta_5000.rds"))
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
# read results
dir_out1 <- file.path("Output/Study1/Sim1", paste0("I", test.length, "_", pform, "MST"))
dir_out2 <- file.path("Output/Study1/Sim2", paste0("I", test.length, "_", pform, "MST"))
mstBU <- readRDS(file.path(dir_out1, "mstBU.rds"))

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

# read results
bias_anal1 <- readRDS(file.path(dir_out1, "bias_anal1.rds"))
csee_anal1 <- readRDS(file.path(dir_out1, "csee_anal1.rds"))

##----------------------------------------------------------------------------
# Analytical evaluation 2: using Park et al. (2017)
# read results
csee_anal2 <- readRDS(file.path(dir_out1, "csee_anal2.rds"))

##----------------------------------------------------------------------------
# Simulation 1: Equated NC scoring method
source.files <- file.path(getwd(), "R", src.files)
nrep <- 100

# implement runMST_rep
simEOS_1000 <- runMST_rep(nrep, theta_1000, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE)
simEOS_5000 <- runMST_rep(nrep, theta_5000, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE)

# save results
saveRDS(simEOS_1000, file=file.path(dir_out2, "simEOS_1000.rds"))
saveRDS(simEOS_5000, file=file.path(dir_out2, "simEOS_5000.rds"))

# read results
# simEOS_1000 <- readRDS(file.path(dir_out2, "simEOS_1000.rds"))
# simEOS_5000 <- readRDS(file.path(dir_out2, "simEOS_5000.rds"))

# classification accuracy
ccrEOS_1000 <- ccr(x=simEOS_1000, cut_thread)
ccrEOS_5000 <- ccr(x=simEOS_5000, cut_thread)

# write out CCR table
write.csv(ccrEOS_1000, file.path(dir_out2, "ccrEOS_1000.csv"))
write.csv(ccrEOS_5000, file.path(dir_out2, "ccrEOS_5000.csv"))

# remove simulation results
rm(simEOS_1000)
rm(simEOS_5000)

##----------------------------------------------------------------------------
# Simulation 2: MLE scoring
module.mat <- mstBU$itemvar.df[[1]]
route.scoring <- "EAP"
final.scoring <- "ML"
range <- range.theta 
cut.mat <- cutoff_mat(cut.score, n.stage, n.module)
source.files <- file.path(getwd(), "R", src.files)

# implement runMST_MLErep
simMLE_1000 <- runMST_MLErep(nrep, theta=theta_1000, item.pool, module.mat, ini.module=1, route.map, cut.mat, 
                             D, route.scoring=route.scoring, final.scoring=final.scoring, range=range, full.out=FALSE, iseed=0L, n.core=2)
simMLE_5000 <- runMST_MLErep(nrep, theta=theta_5000, item.pool, module.mat, ini.module=1, route.map, cut.mat, 
                             D, route.scoring=route.scoring, final.scoring=final.scoring, range=range, full.out=FALSE, iseed=0L, n.core=2)

# save results
saveRDS(simMLE_1000, file=file.path(dir_out2, "simMLE_1000.rds"))
saveRDS(simMLE_5000, file=file.path(dir_out2, "simMLE_5000.rds"))

# read results
simMLE_1000 <- readRDS(file.path(dir_out2, "simMLE_1000.rds"))
simMLE_5000 <- readRDS(file.path(dir_out2, "simMLE_5000.rds"))

# classification accuracy
ccrMLE_1000 <- ccr(x=simMLE_1000, cut_thread)
ccrMLE_5000 <- ccr(x=simMLE_5000, cut_thread)

# write out CCR table
write.csv(ccrMLE_1000, file.path(dir_out2, "ccrMLE_1000.csv"))
write.csv(ccrMLE_5000, file.path(dir_out2, "ccrMLE_5000.csv"))

# remove simulation results
rm(simMLE_1000)
rm(simMLE_5000)

##----------------------------------------------------------------------------
# Classification accuracy for analytical method
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
quadrature <- list(theta, w)

CCR_anal1 <- matrix(NA, length(cut_thread), 2)
CCR_anal2 <- matrix(NA, length(cut_thread), 2)
colnames(CCR_anal1) <- colnames(CCR_anal2) <- c("CCR", "TER")
rownames(CCR_anal1) <- rownames(CCR_anal2) <- paste0("CUT (", cut_thread, ")")
for(i in 1:length(cut_thread)) {
  ccrs <- cacIRT::Rud.D(cutscore=cut_thread[i], quadrature, sem=csee_anal1)$Marginal[1]
  CCR_anal1[i, ] <- c(ccrs, 1-ccrs)
  
  ccrs <- cacIRT::Rud.D(cutscore=cut_thread[i], quadrature, sem=csee_anal2)$Marginal[1]
  CCR_anal2[i, ] <- c(ccrs, 1-ccrs)
}

# save results
write.csv(CCR_anal1, file.path(dir_out2, "CCR_analEOS.csv"))
write.csv(CCR_anal2, file.path(dir_out2, "CCR_analMLE.csv"))

##----------------------------------------------------------------------------
# Organize the results
ccrEOS <- cbind(CCR_anal1, ccrEOS_1000[, c(1, 4)], ccrEOS_5000[, c(1, 4)])
ccrMLE <- cbind(CCR_anal2, ccrMLE_1000[, c(1, 4)], ccrMLE_5000[, c(1, 4)])
ccr_table <- rbind(ccrEOS, ccrMLE)
colnames(ccr_table) <- c("CCR (Anal)", "TER (Anal)",
                         "CCR (Sim1000)", "TER (Sim1000)", 
                         "CCR (Sim5000)", "TER (Sim5000)")
rownames(ccr_table) <- c("CUT.EOS (0)", "CUT.EOS (0.524)",
                         "CUT.MLE (0)", "CUT.MLE (0.524)")

# save results
write.csv(ccr_table, file.path(dir_out2, "CCR_table.csv"))



