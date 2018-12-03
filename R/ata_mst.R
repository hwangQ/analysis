# a function for automated test assembly of MST

ata_mst <- function(item.pool, constraints, target.dens, target.max, relative=FALSE, rel.tol=c(0.7, -0.5),
					theta, D=1.702, divide.D=FALSE, lp.control=list(timeout=2, epsint=0.1, mip.gap=c(0.1, 0.05))) {

	# set conditions
	cats <- item.pool$CATEGORY
	model <- item.pool$MODEL
	item.id <- item.pool$ID
	df_dc <- item.pool[item.pool$CATEGORY == 2, ]
	a <- df_dc$PARAM.1
	b <- df_dc$PARAM.2
	g <- df_dc$PARAM.3
	prm_dc <- list(a, b, g)
	cont <- item.pool$CLASS_RE

	C <- length(unique(cont)) # number of content categories
	I <- length(cont) # total number of items in a pool
	J <- length(theta) # total number of ability points to be evaluated
	nmod <- constraints$n.module
	nstg <- constraints$n.stage
	M <- I + I * nmod[1] + 1 # location of an real number of decision variable 
	F <- 1 + nmod[1] # total number of modules in a test
	m.size <- constraints$size.module
	cont_require <- constraints$content
	Nc_rt <- cont_require[, 2]
	Nc_stg2 <- cont_require[, 3]


	# decompose the target values of a list 
	target_rt <- target.dens[[1]]
	for(i in 1:nmod[1]) {
		assign(paste0("target_stg2_", i), target.dens[[2]][[i]])
	}

	if(relative) {
	
		# rescale the densities
		R_rt <- target_rt / target_rt[1]
		for(i in 1:nmod[1]) {
			assign(paste0("R_stg2_", i), get(paste0("target_stg2_", i)) / target_rt[1])
		}
		R_stg2 <- lapply(1:nmod[1], function(i) get(paste0("R_stg2_", i)))
	
	} else {
	
		# decompose the target maximum values of a list 
		max_rt <- target.max[[1]]
		for(i in 1:nmod[1]) {
			assign(paste0("max_stg2_", i), target.max[[2]][i])
		}
	
		# find constant to be multiplied
		const_rt <- max_rt / max(target_rt)
		for(i in 1:nmod[1]) {
			assign(paste0("const_stg2_", i), get(paste0("max_stg2_", i)) / max(get(paste0("target_stg2_", i))))
		}

		# rescale the densities
		target_rt_sc <- target_rt * const_rt
		for(i in 1:nmod[1]) {
			assign(paste0("target_stg2_", i, "_sc"), get(paste0("target_stg2_", i)) * get(paste0("const_stg2_", i)))
		}
		target_stg2_sc <- vector('list', nmod[1])
		for(i in 1:nmod[1]) {
			target_stg2_sc[[i]] <- get(paste0("target_stg2_", i, "_sc"))
		}
		names(target_stg2_sc) <- paste0("m", 1:nmod[1])

	}
	
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

	# create a matrix of item information 
	info <- test.info(x=df_bank, theta=theta, D=D)$itemInfo

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

	# set control parameters: maximization problem; integer tolerance is et to 0.1; 
	# absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
	sense <- ifelse(relative == TRUE, "max", "min")
	lp.control(lprec=sim_mod, sense=sense, timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)

	# constraint: type of decision variable
	set.type(lprec=sim_mod, columns=1:(M-1), type="binary")
	set.type(lprec=sim_mod, columns=M, type="real")
	set.bounds(lprec=sim_mod, lower=rep(0, M-1), upper=rep(1, M-1), columns=1:(M-1))

	# constraint: module length
	for(i in 1:F) {
		if(i == 1) {
			add.constraint(lprec=sim_mod, xt=rep(1, I), type="=", rhs=m.size[1], indices=(I*(i-1)+1):(I*i))
		} else {
			add.constraint(lprec=sim_mod, xt=rep(1, I), type="=", rhs=m.size[2], indices=(I*(i-1)+1):(I*i))
		}
	}

	# constraint: content category
	for(f in 1:F) {
		if(f == 1) {
			for(i in 1:length(Nc_rt)) {
				add.constraint(lprec=sim_mod, xt=rep(1, length(Vc[[i]])), type=">=", rhs=Nc_rt[i], indices=(f*I-I) + Vc[[i]])
			}
		} else {
			for(i in 1:length(Nc_stg2)) {
				add.constraint(lprec=sim_mod, xt=rep(1, length(Vc[[i]])), type=">=", rhs=Nc_stg2[i], indices=(f*I-I) + Vc[[i]])
			}
		}
	}

	# constraint: no overlap with a router module
	for(f in 1:nmod[1]) {
		for(i in 1:I) {
			add.constraint(lprec=sim_mod, xt=rep(1, 2), type="<=", rhs=1, indices=c(i, (f*I+i)))
		}
	}

	if(relative) {
		
		# constraint: target information (relative target)
		for(f in 1:F) {
			if(f == 1) {
				for(i in 1:length(theta)) {
					add.constraint(lprec=sim_mod, xt=c(info[, i], -R_rt[i]), type="<=", rhs=rel.tol[1], indices=c((f*I-I+1):(f*I), M))
					add.constraint(lprec=sim_mod, xt=c(info[, i], -R_rt[i]), type=">=", rhs=rel.tol[2], indices=c((f*I-I+1):(f*I), M))
				}
			} else {
				for(i in 1:length(theta)) {
					add.constraint(lprec=sim_mod, xt=c(info[, i], -R_stg2[[f-1]][i]), type="<=", rhs=rel.tol[1], indices=c((f*I-I+1):(f*I), M))
					add.constraint(lprec=sim_mod, xt=c(info[, i], -R_stg2[[f-1]][i]), type=">=", rhs=rel.tol[2], indices=c((f*I-I+1):(f*I), M))
				}
			}
		}
		
	} else {
	
		# constraint: target information (absolute target)
		for(f in 1:F) {
			if(f == 1) {
				for(i in 1:length(theta)) {
					add.constraint(lprec=sim_mod, xt=c(info[, i], -1), type="<=", rhs=(target_rt_sc[i]), indices=c((f*I-I+1):(f*I), M))
					add.constraint(lprec=sim_mod, xt=c(info[, i], 1), type=">=", rhs=(target_rt_sc[i]), indices=c((f*I-I+1):(f*I), M))
				}
			} else {
				for(i in 1:length(theta)) {
					add.constraint(lprec=sim_mod, xt=c(info[, i], -1), type="<=", rhs=(target_stg2_sc[[f-1]][i]), indices=c((f*I-I+1):(f*I), M))
					add.constraint(lprec=sim_mod, xt=c(info[, i], 1), type=">=", rhs=(target_stg2_sc[[f-1]][i]), indices=c((f*I-I+1):(f*I), M))
				}
			}
		}
	
	}

	# solve the model
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
	rownames(select_df) <- c(paste0("item.", 1:I), "Total")
	colnames(select_df) <- c("router", paste0("stage2_m", 1:nmod[1]), "sum")
	
	df_sim_list <- vector('list', F)
	names(df_sim_list) <- c("router", paste0("stage2_m", 1:nmod[1]))
	class_list <- vector('list', F)
	for(f in 1:F) {
		df_sim_list[[f]] <- df_bank[sim_opt_list[[f]] == 1, ]
		x <- item.pool$CLASS[sim_opt_list[[f]] == 1]
		df_sim_list[[f]]$CLASS <- x
	}

	# test information for each module
	info_rt <- test.info(x=df_sim_list[[1]], theta=theta, D=D)$testInfo
	for(i in 1:nmod[1]) {
		assign(paste0("info_stg2_", i), test.info(x=df_sim_list[[i+1]], theta=theta, D=D)$testInfo)
	}
	info_list <- vector('list', F)
	info_list[[1]] <- info_rt
	for(i in 1:nmod[1]) {
		info_list[[1+i]] <- get(paste0("info_stg2_", i))
	}
	names(info_list) <- c("router", paste0("stage2_m", 1:nmod[1]))
	info.obs=list(router=info_rt, stage2=info_list[2:F])
	names(info.obs[[2]]) <- paste0("m", 1:nmod[1])
	
	
	if(relative) {
	
		# a list of rescaled target distribution
		target_rt_sc <- sim_opt[M] * R_rt
		target_stg2_sc <- lapply(R_stg2, "*", sim_opt[M])
		target_list <- c(list(target_rt_sc), target_stg2_sc)
		names(target_list) <- c("router", paste0("stage2_m", 1:nmod[1]))
	
	} else {
	
		# a list of rescaled target distribution
		target_list <- c(list(target_rt_sc), target_stg2_sc)
		names(target_list) <- c("router", paste0("stage2_m", 1:nmod[1]))
	
	}
	
	# sum of squre differences betwwen targe and module 
	ssd <- sum(sapply(1:length(info_list), function(i) sum((target_list[[i]] - info_list[[i]])^2)))
	
	# a type of target distribution
	target.type <- ifelse(relative == TRUE, "relative", "absolute")
	
	rr <- list(var.opt=sim_opt, itemvar.list=sim_opt_list, itemvar.df=select_df, 
			prm.df=df_sim_list, info.obs=info.obs, target.sc=list(router=target_rt_sc, stage2=target_stg2_sc), 
			theta=theta, ssd=ssd, solution=solution, metainfo=list(target.type=target.type, n.module=nmod, n.stage=nstg, size.module=m.size)
			)

	class(rr) <- "ata_mst"
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

 
 
