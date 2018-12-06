######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 12/05/2018
#Title: Find the optimal values of target maximums for 1-3-3 MST design structure 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(cacIRT)
library(tidyverse)

##----------------------------------------------------------------------------
# Read source code
src.files <- list.files("R/", pattern="*.R")
purrr::walk(src.files, function(x) source(file.path("R", x)))

##----------------------------------------------------------------------------
# set conditions
test.length <- 32
n.stage <- 3
pform <- "122"
theta <- seq(-4, 4, .1)
# RDP <- c(-0.44, 0.44)
RDP <- 0
post <- c(-2.0, RDP, 2.0) 
RDPList <- list(NA, RDP, RDP)
D <- 1
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.1, mip.gap=c(0.1, 0.05))
route.map <- read.csv(paste0("Input/Routing_map_P", pform, ".csv"), header=FALSE) # read a routing map
content1 <- read.csv(paste0("Input/Content_TD_Class1.csv"))
content2 <- read.csv(paste0("Input/Content_TD_Class2.csv"))
content <- list(content1[, 2], content2[, 2])
# path.group <- list(g1=1, g2=c(2, 3), g3=4, g4=c(5, 6), g5=7) # 1-3-3 MST
# path.group <- list(g1=1, g2=2, g3=3) # 1-3 MST
path.group <- list(g1=1, g2=c(2, 3), g3=4) # 1-2-2 MST
minmod.p <- 0.3 
constraints <- list(route.map=route.map, post=post, path.group=path.group, test.length=test.length, content=content, minmod.p=minmod.p)

##----------------------------------------------------------------------------
# read item bank and content category information
bank_info <- readRDS("Input/item_pool_2.rds")
bank_info <- data.frame(ID=1:nrow(bank_info), bank_info)
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS1", "CLASS2")

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- dplyr::select(bank_info, ID, CATEGORY, MODEL, CLASS1, CLASS2, dplyr::starts_with("PARAM"))
item.pool <- bank_info

##----------------------------------------------------------------------------
# a function for automated test assembly of MST
ata_mstTD(item.pool, constraints, theta, D=D, divide.D=FALSE, lp.control=lp.control)




ata_mstTD <- function(item.pool, constraints, theta, D=1.702, divide.D=FALSE, 
                      lp.control=list(timeout=60, epsint=0.1, mip.gap=c(0.1, 0.05))) {
  
  ##------------------------------------------
  # set conditions
  cats <- item.pool$CATEGORY
  model <- item.pool$MODEL
  item.id <- item.pool$ID
  df_dc <- item.pool[item.pool$CATEGORY == 2, ]
  a <- df_dc$PARAM.1
  b <- df_dc$PARAM.2
  g <- df_dc$PARAM.3
  prm_dc <- list(a, b, g)
  cont1 <- item.pool$CLASS1
  cont2 <- item.pool$CLASS2
  
  C1 <- length(constraints$content[[1]]) # number of the first content categories
  C2 <- length(constraints$content[[2]]) # number of the second content categories
  I <- length(cont1) # total number of items in a pool
  panel.info <- panel_info(constraints$route.map)
  config.info <- panel.info$config
  pathway <- panel.info$pathway
  nstg <- panel.info$n.stage
  nmod <- panel.info$n.module
  n.pathway <- nrow(pathway)
  path.group <- constraints$path.group
  n.group <- length(path.group)
  M <-  I * sum(nmod) + 1 # location of an real number of decision variable for an objective function
  F <- sum(nmod) # total number of modules in a test
  test.length <- constraints$test.length
  cont_require <- constraints$content
  Nc_path1 <- cont_require[[1]] * test.length
  Nc_path2 <- cont_require[[2]] * test.length
  post <- constraints$post
  RDP <- post[c(-1, -length(post))] 
  n.rdp <- (n.stage - 1) * length(RDP) # the number of RPD points where two adjacent modules intersect
  minmod.p <- constraints$minmod.p
  
  ##------------------------------------------
  # warning messages
  if(test.length != sum(Nc_path1)) stop("A test length should be the same with the sum of first content categoires for a pathway")
  if(test.length != sum(Nc_path2)) stop("A test length should be the same with the sum of second content categoires for a pathway")
  
  ##------------------------------------------
  # create a data frame for item parameters in a bank
  a <- prm_dc[[1]]
  b <- prm_dc[[2]]
  g <- prm_dc[[3]]
  if(divide.D) {
    prm_dc_list <- list(a=a/D, b=b, g=g)
  } else {
    prm_dc_list <- list(a=a, b=b, g=g)
  }
  df_bank <- shape_df(par.dc=prm_dc_list, item.id=item.id, cats=cats, model=model)
  
  # create a matrix of item information for each subpopulation group of route
  theta_list <- subpop(post, n.stage=nstg)
  info_list <- vector('list', n.group)
  for(i in 1:n.group) {
    info_list[[i]] <- test.info(x=df_bank, theta=theta_list[[i]], D=D)$itemInfo
  }
  
  ##------------------------------------------
  # assign item to each category of content
  Vc1 <- list()
  for(i in 1:C1) {
    Vc1[[i]] <- c(1:I)[cont1 == i]
  }
  Vc2 <- list()
  for(i in 1:C2) {
    Vc2[[i]] <- c(1:I)[cont2 == i]
  }
  
  ##--------------------------------------------------------------------
  # simultaneous assembly
  # make a model
  sim_mod <- make.lp(nrow=0, ncol=M+n.rdp)
  
  # objective function
  set.objfn(lprec=sim_mod, obj=1, indices=M)
  
  # set control parameters: maximization problem; integer tolerance is et to 0.1; 
  # absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
  lp.control(lprec=sim_mod, sense="max", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)
  
  # constraint: type of decision variable
  set.type(lprec=sim_mod, columns=1:(M-1), type="binary")
  set.type(lprec=sim_mod, columns=M:(M+n.rdp), type="real")
  set.bounds(lprec=sim_mod, lower=rep(0, M-1), upper=rep(1, M-1), columns=1:(M-1))
  
  # constraint: test length
  for(w in 1:n.pathway) {
    temp <- cbind((pathway[w, ] * I - I + 1), (pathway[w, ] * I))
    indices <- NULL
    for(s in 1:n.stage) {
      indices <- c(indices, temp[s, 1]:temp[s, 2])
    }
    add.constraint(lprec=sim_mod, xt=rep(1, I * n.stage), type="=", rhs=test.length, indices=indices)
  }
  
  # constraint: first content category
  for(w in 1:n.pathway) {
    for(i in 1:length(Nc_path1)) {
      temp <- pathway[w, ]*I - I
      indices <- NULL
      for(s in 1:n.stage) {
        indices <- c(indices, temp[s] + Vc1[[i]])
      }
      add.constraint(lprec=sim_mod, xt=rep(1, length(Vc1[[i]]) * n.stage), type=">=", rhs=Nc_path1[i], indices=indices)
    }
  }
  
  # constraint: second content category
  for(w in 1:n.pathway) {
    for(i in 1:length(Nc_path2)) {
      temp <- pathway[w, ]*I - I
      indices <- NULL
      for(s in 1:n.stage) {
        indices <- c(indices, temp[s] + Vc2[[i]])
      }
      add.constraint(lprec=sim_mod, xt=rep(1, length(Vc2[[i]]) * n.stage), type=">=", rhs=Nc_path2[i], indices=indices)
    }
  }
  
  # constraint: no overlap between stages
  for(w in 1:n.pathway) {
    for(i in 1:I) {
      indices <- (pathway[w, ]*I-I) + i
      add.constraint(lprec=sim_mod, xt=rep(1, nstg), type="<=", rhs=1, indices=indices)
    }
  }
  
  # constraint: two adjacent modules have the same module information at RDP
  info.RDP <- test.info(x=df_bank, theta=RDP, D=D)$itemInfo
  i <- 1
  for(s in 2:n.stage) {
    for(m in 1:(nmod[s]-1)) {
      index.1 <- (unique(pathway[, s])[m] * I - I + 1):(unique(pathway[, s])[m] * I)
      index.2 <- (unique(pathway[, s])[m+1] * I - I + 1):(unique(pathway[, s])[m+1] * I)
      indices.1 <- c(index.1, M + i)
      indices.2 <- c(index.2, M + i)
      add.constraint(lprec=sim_mod, xt=c(info.RDP[, m], -1), type="=", rhs=0, indices=indices.1)
      add.constraint(lprec=sim_mod, xt=c(info.RDP[, m], -1), type="=", rhs=0, indices=indices.2)
      i <- i + 1
    }
  }
  
  # constraint: minimum module length
  for(f in 1:F) {
    add.constraint(lprec=sim_mod, xt=rep(1, I), type=">=", rhs=(test.length * minmod.p), indices=(I*(f-1)+1):(I*f))
  }
  
  # constraint: target information (relative target)
  for(w in 1:n.pathway) {
    for(g in 1:n.group) {
      if(w %in% path.group[[g]]) {
        for(i in 1:length(theta_list[[g]])) {
          temp <- cbind((pathway[w, ] * I - I + 1), (pathway[w, ] * I))
          indices <- NULL
          for(s in 1:n.stage) {
            indices <- c(indices, temp[s, 1]:temp[s, 2])
          }
          indices <- c(indices, M)
          add.constraint(lprec=sim_mod, xt=c(rep(info_list[[g]][, i], n.stage), -1), type=">=", rhs=0, indices=indices)
        }
      }
    }
  }
  
  # solve teh model
  res_flag_sim <- solve(sim_mod)
  
  # the integer value containing the status code, for example, 0: "optimal solution found"
  eval.model <- res_flag_sim
  if(eval.model == 0) {
    solution <- "optimal"
    print("Optimal soluation is found")
  }
  if(eval.model == 1) {
    solution <- "suboptimal"
    print("The model is sub-optimal")
  }
  if(eval.model > 1) stop(paste0("A status code for this model is ", eval.model))
  
  # retrieve the values of the decision variables 
  sim_opt <- get.variables(sim_mod)
  
  ##------------------------------------------
  # assemble the modules based on the optimized results of decision variables
  sim_opt_list <- vector('list', F)
  for(f in 1:F) {
    sim_opt_list[[f]] <- sim_opt[(f*I-I+1):(f*I)]
  }
  select_mat <- do.call('cbind', sim_opt_list)
  select_df <- cbind(select_mat, rowSums(select_mat))
  select_df <- rbind(select_df, colSums(select_df))
  rownames(select_df) <- c(paste0("item.", 1:I), "Sum")
  colnames(select_df) <- c(paste0("module.", 1:F), "Sum")
  
  df_sim_list <- vector('list', F)
  names(df_sim_list) <- paste0("module.", 1:F)
  class1_list <- vector('list', F)
  class2_list <- vector('list', F)
  for(f in 1:F) {
    df_sim_list[[f]] <- df_bank[sim_opt_list[[f]] == 1, ]
    x1 <- item.pool$CLASS1[sim_opt_list[[f]] == 1]
    x2 <- item.pool$CLASS2[sim_opt_list[[f]] == 1]
    df_sim_list[[f]]$CLASS1 <- x1
    df_sim_list[[f]]$CLASS2 <- x2
  }
  
  ##------------------------------------------
  # test information for each module
  info_mod <- lapply(1:F, function(f) test.info(x=df_sim_list[[f]], theta=theta, D=D)$testInfo)
  names(info_mod) <- paste0("module.", 1:F)
  
  ## return results 
  rr <- list(var.opt=sim_opt, itemvar.list=sim_opt_list, itemvar.df=list(select_mat, select_df), 
             prm.df=df_sim_list, info.mod=info_mod, theta=theta, solution=solution, 
             metainfo=constraints)
  
  rr
  
}
  
