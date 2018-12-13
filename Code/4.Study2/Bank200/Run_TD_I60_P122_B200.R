######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 12/07/2018
#Title: ATA using a top-down approah 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(cacIRT)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

##----------------------------------------------------------------------------
# Read source code
src.files <- list.files("R/", pattern="*.R")
purrr::walk(src.files, function(x) source(file.path("R", x)))

##----------------------------------------------------------------------------
# set the item bank size
bank_size <- 200

# read item bank and content category information
bank_info <- readRDS(paste0("Input/item_pool_2_B", bank_size, ".rds"))
bank_info <- data.frame(ID=1:nrow(bank_info), bank_info)
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS1", "CLASS2")

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- dplyr::select(bank_info, ID, CATEGORY, MODEL, CLASS1, CLASS2, dplyr::starts_with("PARAM"))
item.pool <- bank_info

##----------------------------------------------------------------------------
# set conditions
n.core <- 3 # the number of cores to be used for parallel analysis
test.length <- 60
n.stage <- 3
pform <- "122"
theta <- seq(-4, 4, .1)
D <- 1.0
divide.D <- FALSE
with.end <- FALSE
equal.info <- TRUE
lp.control <- list(timeout=180, epsint=0.2, mip.gap=c(0.1, 0.05))
content1 <- read.csv(paste0("Input/Content_TD_Class1.csv"))
content2 <- read.csv(paste0("Input/Content_TD_Class2.csv"))
content <- list(content1[, 2], content2[, 2])
minmod.p <- 0.2
if(pform == "133") path.group <- list(g1=1, g2=c(2, 3), g3=4, g4=c(5, 6), g5=7) # 1-3-3 MST
if(pform == "13") path.group <- list(g1=1, g2=2, g3=3) # 1-3 MST
if(pform == "122") path.group <- list(g1=1, g2=c(2, 3), g3=4) # 1-2-2 MST
route.map <- read.csv(paste0("Input/Routing_map_P", pform, ".csv"), header=FALSE) # read a routing map
if(pform %in% c("13", "133")) RDP_mat <- data.matrix(expand.grid(seq(-0.8, -0.1, by=0.02), seq(0.1, 0.8, by=0.02)))
if(pform == "122") RDP_mat <- cbind(seq(-0.7, 0.7, by=0.02))

constraints <- vector('list', nrow(RDP_mat)) # a list of constraints
for(i in 1:nrow(RDP_mat)) {
  RDP <- RDP_mat[i, ]
  post <- c(-2.0, RDP, 2.0)
  constraints[[i]] <- list(route.map=route.map, post=post, path.group=path.group, 
                           test.length=test.length, content=content, minmod.p=minmod.p,
                           with.end=with.end, equal.info=equal.info)
}

##----------------------------------------------------------------------------
# Automated test assembly using a top-down approach
# set the number of cpu cores
left.core <- 4 - n.core
numCores <- parallel::detectCores() - left.core

# create a parallel processesing cluster
cl = parallel::makeCluster(numCores, type="PSOCK")

# load packages
parallel::clusterEvalQ(cl, library(lpSolveAPI))
parallel::clusterEvalQ(cl, library(reshape2))

# load some specific variable names into processing cluster
parallel::clusterExport(cl, c("item.pool", "constraints", "theta", "D", "divide.D", "lp.control"))

# load source files used for the replications into processing cluster
source.files <- file.path(getwd(), "R", src.files)
lapply(1:length(source.files), function(i) parallel::clusterCall(cl, fun=source, source.files[i]))

# run ata_mstTD for all combinations of RDPs
f <- function(i) ata_mstTD(item.pool, constraints[[i]], theta, D, divide.D, lp.control)
mstTD <- pbapply::pblapply(X=1:length(constraints), FUN=f, cl=cl) # to see the progress bar

# finish
parallel::stopCluster(cl)

# save results
dir_out <- file.path(paste0("Output/Study2/Bank", bank_size), paste0("I", test.length, "_", pform, "MST"))
saveRDS(mstTD, file=file.path(dir_out, "mstTD.rds"))

# read results
# mstTD <- readRDS(file.path(dir_out, "mstTD.rds"))

# check the number of cases that the ATA process failed
failed_num <- sum(purrr::map_lgl(mstTD, is.null))
failed_rate <- failed_num / length(mstTD)
failed <- c(failed_num, failed_rate) %>% 
  setNames(nm = c("N", "(%)"))

# save the failed case and rates
write.csv(failed, file.path(dir_out, "failed.csv"))

##----------------------------------------------------------------------------
# Analytical evaluation
# set condtions
range.theta <- c(-5, 5)
interval <- c(-2.5, 2.5)

# create an empty containers
cond_moments <- vector('list', length(mstTD))
mrel <- rep(NA, length(mstTD))
ave_csee <- rep(NA, length(mstTD))
max_csee <- rep(NA, length(mstTD))

# evaluation
# initialize the progress bar
pb <- utils::winProgressBar(title="Evaluation Progress: 0%", label="No replication is done", min=0, max=100, initial=0)
for(i in 1:length(mstTD)) {
  
  if(is.null(mstTD[[i]])) { 
    
    next
    
  } else {
    
    # read the assembled test forms
    ata_forms <- mstTD[[i]]$prm.df
    RDP <- RDP_mat[i, ]
    RDPList <- list(NA, RDP, RDP)
    
    # read a routing map
    panel.info <- panel_info(route.map)
    config.info <- panel.info$config
    n.module <- panel.info$n.module
    pathway <- panel.info$pathway
    
    # estimate the observed equated scores across all (sub) pathways
    eos_list <- est_eos(ata_forms, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)
    
    # find points where two adjacent TIFs intersect across all stages
    cut.score <- cutoff(ata_forms, route.map, RDP=RDPList, D, range.theta, interval)
    
    # conditional raw score distribution of each module given an ability value
    cond_dist <- ability_dist_sep(ata_forms, theta, route.map, D)
    
    # conditional ability distribution of total-test ability estimates given a true ability
    x <- cond_dist
    eos.path <- eos_list$eos_path
    n.path <- eos_list$n_path
    # joint.dist <- ability_dist(x, eos.path, n.path, cut.score=RDPList)
    joint.dist <- ability_dist(x, eos.path, n.path, cut.score=cut.score)
    
    # calculate a mean and sd of conditional ability distribution at each ability value
    if(n.stage == 2) {
      cond_moments[[i]] <- sapply(joint.dist$joint.dist$stage2, cal_moments, node=eos_list$eos_path$stage2)
    }
    if(n.stage == 3) {
      cond_moments[[i]] <- sapply(joint.dist$joint.dist$stage3, cal_moments, node=eos_list$eos_path$stage3)
    }
    
    ##----------------------------------------------------------------------------
    # Compute three objective functions
    # generate weights
    w <- dnorm(theta, 0, 1)
    w <- w/sum(w)
    
    # compute objective function values
    var.cond <- cond_moments[[i]][2, ]
    mrel[i] <- objfn(obj="mrel", var.cond, var.pop=1, w=w, range=c(-2, 2))
    ave_csee[i] <- objfn(obj="ave.se", var.cond, var.pop=1, w=w, range=c(-2, 2))
    max_csee[i] <- objfn(obj="max.se", var.cond, var.pop=1, w=w, range=c(-2, 2))
    
  }
  
  # modify the progress bar (using setWinProgressBar function)
  info <- sprintf("Evaluation Progress: %d%%", round((i/length(mstTD))*100))
  utils::setWinProgressBar(pb, i/length(mstTD)*100, title=info, label=paste0("Evaluation of the assembled MST ", i, " is done"))
  
}

###### Closing the Progress Bar
close(pb)

# combine all objective function results
obj_res <- data.frame(mrel=mrel, ave_csee=ave_csee, max_csee=max_csee)

# save results
saveRDS(cond_moments, file=file.path(dir_out, "cond_moments.rds"))
saveRDS(obj_res, file=file.path(dir_out, "obj_res.rds"))

# read results
# cond_moments <- readRDS(file.path(dir_out, "cond_moments.rds"))
# obj_res <- readRDS(file.path(dir_out, "obj_res.rds"))

##----------------------------------------------------------------------------
# summarize the results of objective functions in a specified order
obj_df <- summary_obj(obj_res, order=1:length(mstTD))

# save the summary of the results
write.csv(obj_df, file.path(dir_out, "obj_df.csv"))



