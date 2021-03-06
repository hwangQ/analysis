######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 12/21/2018
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
bank_size <- 400

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
pform <- "133"
theta <- seq(-4, 4, .1)
D <- 1.0
divide.D <- FALSE
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
                           equal.info=equal.info)
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
route.method <- "DPI"

# evaluation
# set the number of cpu cores
left.core <- 4 - n.core
numCores <- parallel::detectCores() - left.core

# create a parallel processesing cluster
cl = parallel::makeCluster(numCores, type="PSOCK")

# load packages
parallel::clusterEvalQ(cl, library(lpSolveAPI))
parallel::clusterEvalQ(cl, library(reshape2))

# load some specific variable names into processing cluster
parallel::clusterExport(cl, c("route.method", "theta", "range.theta", "interval", "D"))

# load source files used for the replications into processing cluster
source.files <- file.path(getwd(), "R", src.files)
lapply(1:length(source.files), function(i) parallel::clusterCall(cl, fun=source, source.files[i]))

# run anal_evalMST for all assembled MSTs
eval_mstTD <- pbapply::pblapply(X=mstTD, FUN=anal_evalMST, route.method=route.method, theta=theta, 
                                range.theta=range.theta, interval=interval, D=D, cl=cl) 

##----------------------------------------------------------------------------
# extract measurement precision and objective function values
cond_moments <- purrr::map(eval_mstTD, .f=function(x) x$cond_moments)
mrel <- purrr::map_dbl(eval_mstTD, .f=function(x) x$object.fn$mrel)
ave_csee <- purrr::map_dbl(eval_mstTD, .f=function(x) x$object.fn$ave_csee)
max_csee <- purrr::map_dbl(eval_mstTD, .f=function(x) x$object.fn$max_csee)

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
objRDP_df <- summary_obj(obj_res, order=1:length(mstTD), showRDP=TRUE, RDP_mat=RDP_mat)

# save the summary of the results
write.csv(obj_df, file.path(dir_out, "obj_df.csv"))
write.csv(objRDP_df, file.path(dir_out, "objRDP_df.csv"))

# check the partion of contents
content_df <- content_table(x=mstTD, which.mst=obj_df$loc.mre[c(1:4)], mod.name=c("1M", "2E", "2M", "2H", "3E", "3M", "3H")) 
write.csv(content_df, file.path(dir_out, "content_df.csv"))

##----------------------------------------------------------------------------
# plot test information functions for all routes for each assembled MST
plot(x = mstTD[[1170]], ylab.text = "Information", range.theta, D=D,
     legend.text=c("1M-2E-3E", "1M-2E-3M", "1M-2M-3E", "1M-2M-3M", "1M-2M-3H", "1M-2H-3M", "1M-2H-3H"))

# plot test information functions for all routes for multiple assembled MSTs
plot(x=mstTD, which.mst=obj_df$loc.mre[c(1:4)], range.theta, D, 
     ylab.text = "Information",
     legend.text=c("1M-2E-3E", "1M-2E-3M", "1M-2M-3E", "1M-2M-3M", "1M-2M-3H", "1M-2H-3M", "1M-2H-3H"), 
     legend.size=15, legend.position="right", strip.size=12,
     layout.col=2)

# plot test information functions across stages for multiple assembled MSTs
plot_mif(x=mstTD, which.mst=obj_df$loc.mre[c(1:4)], 
         ylab.text = "Information",
         line.size=1.5, 
         legend.text=c("1M", "2E", "2M", "2H", "3E", "3M", "3H"), 
         legend.size=15, legend.position="right", strip.size=12) 

# plot CSEEs
plot_csee(cond_moments, which.mst=obj_df$loc.mre[c(1:4)], RDP_mat, 
          ylab.text = "Standard Error",
          lab.size=15, axis.size=15, line.size=1.5, legend.size=15, legend.position="right")



