###############################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 12/13/2018
#Title: Plot CSEE and Bias for Study 1 
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

# set conditions to draw plots
test.length <- c("24 items", "24 items", "48 items", "48 items")
pform <- c("1-3 MST", "1-3-3 MST", "1-3 MST", "1-3-3 MST")
theta <- seq(-4, 4, .1)
dir_fd <- c("I24_13MST", "I24_133MST", "I48_13MST", "I48_133MST")
dir_out1 <- file.path("Output/Study1/Sim1", dir_fd)
dir_out2 <- file.path("Output/Study1/Sim2", dir_fd)

## -----------------------------------------------------------------
# 1. CSEE plot
# prepare data to draw CSEE plots
csee_R100 <- vector('list', length(dir_fd))
csee_R1000 <- vector('list', length(dir_fd))
csee_R5000 <- vector('list', length(dir_fd))
for(i in 1:length(dir_fd)) {

  csee_anal1 <- readRDS(file.path(dir_out1[i], "csee_anal1.rds"))
  csee_anal2 <- readRDS(file.path(dir_out1[i], "csee_anal2.rds"))
  myeval_EOS <- readRDS(file.path(dir_out1[i], "myeval_EOS.rds"))
  myeval_MLE <- readRDS(file.path(dir_out1[i], "myeval_MLE.rds"))
  csee_simEOS <- map(myeval_EOS, function(x) sqrt(x$cond.eval[, 2]))
  csee_simMLE <- map(myeval_MLE, function(x) sqrt(x$cond.eval[, 2]))
  
  csee_R100[[i]] <- 
    data.frame(AnalENC=csee_anal1, AnalMLE=csee_anal2, 
               SimENC=csee_simEOS[[1]], SimMLE=csee_simMLE[[1]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="CSEE") %>% 
    mutate(Panel=pform[i], Length=test.length[i])

  csee_R1000[[i]] <- 
    data.frame(AnalENC=csee_anal1, AnalMLE=csee_anal2, 
               SimENC=csee_simEOS[[2]], SimMLE=csee_simMLE[[2]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="CSEE") %>% 
    mutate(Panel=pform[i], Length=test.length[i])
  
  csee_R5000[[i]] <- 
    data.frame(AnalENC=csee_anal1, AnalMLE=csee_anal2, 
               SimENC=csee_simEOS[[3]], SimMLE=csee_simMLE[[3]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="CSEE") %>% 
    mutate(Panel=pform[i], Length=test.length[i])
  
}

# transform the list to the data.frame
csee_R100 <- do.call('rbind', csee_R100) 
csee_R1000 <- do.call('rbind', csee_R1000)
csee_R5000 <- do.call('rbind', csee_R5000) 

# draw CSEE plots
# (1) R = 100
plot_study1(x = csee_R100, y.var = "CSEE", ylim = c(0, 1.6))

# (2) R = 1000
plot_study1(x = csee_R1000, y.var = "CSEE", ylim = c(0, 1.6))

# (3) R = 5000
plot_study1(x = csee_R5000, y.var = "CSEE", ylim = c(0, 1.6))

## -----------------------------------------------------------------
# 2. Conditional Bias plot
# prepare data to draw conditional bias plots
bias_R100 <- vector('list', length(dir_fd))
bias_R1000 <- vector('list', length(dir_fd))
bias_R5000 <- vector('list', length(dir_fd))
for(i in 1:length(dir_fd)) {
  
  bias_anal1 <- readRDS(file.path(dir_out1[i], "bias_anal1.rds"))
  myeval_EOS <- readRDS(file.path(dir_out1[i], "myeval_EOS.rds"))
  myeval_MLE <- readRDS(file.path(dir_out1[i], "myeval_MLE.rds"))
  bias_simEOS <- map(myeval_EOS, function(x) x$cond.eval[, 1])
  bias_simMLE <- map(myeval_MLE, function(x) x$cond.eval[, 1])
  
  bias_R100[[i]] <- 
    data.frame(AnalENC=bias_anal1, SimENC=bias_simEOS[[1]], 
               SimMLE=bias_simMLE[[1]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="Bias") %>% 
    mutate(Panel=pform[i], Length=test.length[i])
  
  bias_R1000[[i]] <- 
    data.frame(AnalENC=bias_anal1, SimENC=bias_simEOS[[2]], 
               SimMLE=bias_simMLE[[2]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="Bias") %>% 
    mutate(Panel=pform[i], Length=test.length[i])
  
  bias_R5000[[i]] <- 
    data.frame(AnalENC=bias_anal1, SimENC=bias_simEOS[[3]], 
               SimMLE=bias_simMLE[[3]], theta=theta) %>% 
    reshape2::melt(variable.name="Method", id.vars="theta", value.name="Bias") %>% 
    mutate(Panel=pform[i], Length=test.length[i])
  
}

# transform the list to the data.frame
bias_R100 <- do.call('rbind', bias_R100) 
bias_R1000 <- do.call('rbind', bias_R1000)
bias_R5000 <- do.call('rbind', bias_R5000) 

# draw the conditional bias plots
# (1) R = 100
plot_study1(x = bias_R100, y.var = "Bias", ylim = c(-1.6, 2.0))

# (2) R = 1000
plot_study1(x = bias_R1000, y.var = "Bias", ylim = c(-1.6, 2.0))

# (3) R = 5000
plot_study1(x = bias_R5000, y.var = "Bias", ylim = c(-1.6, 2.0))



