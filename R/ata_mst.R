# a function for ATA with top-down approach
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
  n.rdp <- (nstg - 1) * length(RDP) # the number of RPD points where two adjacent modules intersect
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
    for(s in 1:nstg) {
      indices <- c(indices, temp[s, 1]:temp[s, 2])
    }
    add.constraint(lprec=sim_mod, xt=rep(1, I * nstg), type="=", rhs=test.length, indices=indices)
  }
  
  # constraint: first content category
  for(w in 1:n.pathway) {
    for(i in 1:length(Nc_path1)) {
      temp <- pathway[w, ]*I - I
      indices <- NULL
      for(s in 1:nstg) {
        indices <- c(indices, temp[s] + Vc1[[i]])
      }
      add.constraint(lprec=sim_mod, xt=rep(1, length(Vc1[[i]]) * nstg), type=">=", rhs=Nc_path1[i], indices=indices)
    }
  }
  
  # constraint: second content category
  for(w in 1:n.pathway) {
    for(i in 1:length(Nc_path2)) {
      temp <- pathway[w, ]*I - I
      indices <- NULL
      for(s in 1:nstg) {
        indices <- c(indices, temp[s] + Vc2[[i]])
      }
      add.constraint(lprec=sim_mod, xt=rep(1, length(Vc2[[i]]) * nstg), type=">=", rhs=Nc_path2[i], indices=indices)
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
  for(s in 2:nstg) {
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
          for(s in 1:nstg) {
            indices <- c(indices, temp[s, 1]:temp[s, 2])
          }
          indices <- c(indices, M)
          add.constraint(lprec=sim_mod, xt=c(rep(info_list[[g]][, i], nstg), -1), type=">=", rhs=0, indices=indices)
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
  if(eval.model > 1) {
    print(paste0("A status code for this model is ", eval.model))
    return(NULL)
  }
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



# a function for ATA with bottom-up approach
ata_mstBU <- function(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
                      lp.control=list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))) {
  
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
  cont <- item.pool$CLASS
  
  C <- length(unique(cont)) # number of content categories
  I <- length(cont) # total number of items in a pool
  panel.info <- panel_info(constraints$route.map)
  config.info <- panel.info$config
  pathway <- panel.info$pathway
  nstg <- panel.info$n.stage
  nmod <- panel.info$n.module
  n.pathway <- nrow(pathway)
  m.size <- constraints$size.module
  M <- I * sum(nmod) + 1 # location of an real number of decision variable
  F <- sum(nmod) # total number of modules in a test
  # test.length <- constraints$test.length
  cont_require <- constraints$content
  
  ##------------------------------------------
  # warning messages
  if(sum(size.module) != sum(cont_require)) stop("A test length should be the same with the sum of content categoires across stages")
  
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
  
  ##------------------------------------------
  # assign item to each category of content
  Vc <- list()
  for(i in 1:C) {
    Vc[[i]] <- c(1:I)[cont == i]
  }
  
  ##--------------------------------------------------------------------
  # simultaneous assembly
  # make a model
  sim_mod <- make.lp(nrow=0, ncol=M)
  
  # objective function
  set.objfn(lprec=sim_mod, obj=1, indices=M)
  
  # set control parameters
  lp.control(lprec=sim_mod, sense="min", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)
  
  # constraint: type of decision variable
  set.type(lprec=sim_mod, columns=1:(M-1), type="binary")
  set.type(lprec=sim_mod, columns=M, type="real")
  set.bounds(lprec=sim_mod, lower=rep(0, M-1), upper=rep(1, M-1), columns=1:(M-1))
  
  # constraint: module length
  for(f in 1:F) {
    idx <- sapply(1:nstg, function(s) f %in% unique(pathway[, s]))
    add.constraint(lprec=sim_mod, xt=rep(1, I), type="=", rhs=m.size[idx], indices=(I*(f-1)+1):(I*f))
  }
  
  # constraint: content category 
  for(f in 1:F) {
    idx <- sapply(1:nstg, function(s) f %in% unique(pathway[, s]))
    Nc_temp <- cont_require[, idx]
    for(i in 1:length(Nc_temp)) {
      add.constraint(lprec=sim_mod, xt=rep(1, length(Vc[[i]])), type=">=", rhs=Nc_temp[i], indices=(f*I-I) + Vc[[i]])
    }
  }
  
  # constraint: no overlap between stages
  for(w in 1:n.pathway) {
    for(i in 1:I) {
      indices <- (pathway[w, ]*I-I) + i
      add.constraint(lprec=sim_mod, xt=rep(1, nstg), type="<=", rhs=1, indices=indices)
    }
  }
  
  # constraint: target information (absolute target)
  for(f in 1:F) {
    idx <- sapply(1:nstg, function(s) f %in% unique(pathway[, s]))
    stg_temp <- which(idx)
    mod_temp <- which(config.info[[stg_temp]] == f)   
    thetas <- targetMIF$theta[[stg_temp]][[mod_temp]]
    info_temp <- test.info(x=df_bank, theta=thetas, D=D)$itemInfo
    MIF_temp <- targetMIF$info[[stg_temp]][[mod_temp]]
    for(i in 1:length(thetas)) {
      add.constraint(lprec=sim_mod, xt=c(info_temp[, i], -1), type="<=", rhs=(MIF_temp[i]), indices=c((f*I-I+1):(f*I), M))
      add.constraint(lprec=sim_mod, xt=c(info_temp[, i], 1), type=">=", rhs=(MIF_temp[i]), indices=c((f*I-I+1):(f*I), M))
    }
  }
  
  ##------------------------------------------
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
  class_list <- vector('list', F)
  for(f in 1:F) {
    df_sim_list[[f]] <- df_bank[sim_opt_list[[f]] == 1, ]
    x <- item.pool$CLASS[sim_opt_list[[f]] == 1]
    df_sim_list[[f]]$CLASS <- x
  }
  
  ##------------------------------------------
  # test information for each module
  info_mod <- lapply(1:F, function(f) test.info(x=df_sim_list[[f]], theta=theta, D=D)$testInfo)
  names(info_mod) <- paste0("module.", 1:F)
  
  ## return results 
  rr <- list(var.opt=sim_opt, itemvar.list=sim_opt_list, itemvar.df=list(select_mat, select_df), 
             prm.df=df_sim_list, info.mod=info_mod, theta=theta, solution=solution, 
             metainfo=list(targetMIF=targetMIF, route.map=route.map)
  )
  
  rr
  
}

 
 
