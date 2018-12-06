######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 07/25/2018
#Title: Find the optimal values of target maximums for 1-3 MST design structure 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(cacIRT)
library(tidyverse)

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

# set condtions for ATA
theta <- seq(-4, 4, .1)
RDP <- c(-0.44, 0.44)
post <- c(-1.5, -0.44, 0.44, 1.5) 
item.pool <- bank_info
n.module <- c(3, 3)
n.stage <- 3
test.length <- 45
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))
pathway <- matrix(c(1, 2, 5,
					1, 2, 6,
					1, 3, 5,
					1, 3, 6,
					1, 3, 7,
					1, 4, 6,
					1, 4, 7), ncol=n.stage, byrow=TRUE)
# path.group <- list(group.1=c(1, 2), group.2=c(3, 4, 5), group.3=c(6, 7))
path.group <- list(group.1=1, group.2=4, group.3=7)
constraints <- list(n.module=n.module, pathway=pathway, path.group=path.group, test.length=test.length, content=cont_require)
min.module <- test.length * 0.3 # minimum module size

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
	nmod <- constraints$n.module
	nstg <- ncol(constraints$pathway)
	n.pathway <- nrow(constraints$pathway)
	path.group <- constraints$path.group
	n.group <- length(path.group)
	M <- sum(I, I * nmod, 1) # location of an real number of decision variable
	F <- sum(1, nmod) # total number of modules in a test
	test.length <- constraints$test.length
	cont_require <- constraints$content
	Nc_path <- cont_require[, "Path"]

	# warning messages
	if(test.length != sum(Nc_path)) stop("A test length should be the same with the sum of content categoires for a pathway")

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
	theta_list <- vector('list', n.group)
	info_list <- vector('list', n.group)
	for(i in 1:n.group) {
		theta_list[[i]] <- seq(post[i], post[i + 1], length.out=5)
		# theta_list[[i]] <- theta_list[[i]][c(-1, -5)]
		info_list[[i]] <- test.info(x=df_bank, theta=theta_list[[i]], D=D)$itemInfo
	}

	# assign item to each category of content
	Vc <- list()
	for(i in 1:C) {
		Vc[[i]] <- c(1:I)[cont == i]
	}

	##--------------------------------------------------------------------
	# simultaneous assembly
	# make a model
	sim_mod <- make.lp(nrow=0, ncol=M+4+18)

	# objective function
	set.objfn(lprec=sim_mod, obj=1, indices=M)

	# set control parameters: maximization problem; integer tolerance is et to 0.1; 
	# absolute MIP gaps is set to 0.1, relative MIP gap is set to 0.5
	lp.control(lprec=sim_mod, sense="max", timeout=lp.control$timeout, epsint=lp.control$epsint, mip.gap=lp.control$mip.gap)

	# constraint: type of decision variable
	set.type(lprec=sim_mod, columns=1:(M-1), type="binary")
	set.type(lprec=sim_mod, columns=M:(M+4+18), type="real")
	set.bounds(lprec=sim_mod, lower=rep(0, M-1), upper=rep(1, M-1), columns=1:(M-1))

	# constraint: test length
	for(w in 1:n.pathway) {
		temp <- cbind((pathway[w, ]*I-I + 1), (pathway[w, ]*I))
		indices <- NULL
		for(s in 1:n.stage) {
			indices <- c(indices, temp[s, 1]:temp[s, 2])
		}
		add.constraint(lprec=sim_mod, xt=rep(1, I * n.stage), type="=", rhs=test.length, indices=indices)
	}

	# constraint: content category
	for(w in 1:n.pathway) {
		for(i in 1:length(Nc_path)) {
			temp <- pathway[w, ]*I - I
			indices <- NULL
			for(s in 1:n.stage) {
				indices <- c(indices, temp[s] + Vc[[i]])
			}
			add.constraint(lprec=sim_mod, xt=rep(1, length(Vc[[i]]) * n.stage), type=">=", rhs=Nc_path[i], indices=indices)
		}
	}

	# constraint: no overlap between stages
	for(w in 1:n.pathway) {
		for(i in 1:I) {
			indices <- (pathway[w, ]*I-I) + i
			add.constraint(lprec=sim_mod, xt=rep(1, nstg), type="<=", rhs=1, indices=indices)
		}
	}

	# constraint: no overlap between modules at the same stage
	# for(s in 2:n.stage) {
	# 	for(i in 1:I) {
	# 		indices <- unique(pathway[, s]) * I - I + i 
	# 		add.constraint(lprec=sim_mod, xt=rep(1, nmod[s-1]), type="<=", rhs=1, indices=indices)
	# 	}
	# }

	# constraint: two adjacent modules have the same module information at RDP
	info.RDP <- test.info(x=df_bank, theta=RDP, D=D)$itemInfo
	i <- 1
	for(s in 2:n.stage) {
		for(m in 1:(nmod[s-1]-1)) {
			index.1 <- (unique(pathway[, s])[m] * I - I + 1):(unique(pathway[, s])[m] * I)
			index.2 <- (unique(pathway[, s])[m+1] * I - I + 1):(unique(pathway[, s])[m+1] * I)
			indices.1 <- c(index.1, M + i)
			indices.2 <- c(index.2, M + i)
			add.constraint(lprec=sim_mod, xt=c(info.RDP[, m], -1), type="=", rhs=0, indices=indices.1)
			add.constraint(lprec=sim_mod, xt=c(info.RDP[, m], -1), type="=", rhs=0, indices=indices.2)
			i <- i + 1
		}
	}

	# constraint: For each group, the module corresponding to that group should have module information greater than eqaul to other modules
	# k <- 1
	# for(s in 2:n.stage) {
	# 	for(m in 1:(nmod[s-1])) {
	# 		for(i in 1:length(theta_list[[m]])) {
	# 			pick.num <- which(1:nmod[s-1] %in% m)
	# 			nopick.num <- which(!1:nmod[s-1] %in% m)
	# 			index <- (unique(pathway[, s])[pick.num] * I - I + 1):(unique(pathway[, s])[pick.num] * I)
	# 			indices <- c(index, M + 4 + k)
	# 			add.constraint(lprec=sim_mod, xt=c(info_list[[pick.num]][, i], -1), type=">=", rhs=0, indices=indices)
	# 			for(j in nopick.num) {
	# 				index <- (unique(pathway[, s])[j] * I - I + 1):(unique(pathway[, s])[j] * I)
	# 				indices <- c(index, M + 4 + k)
	# 				add.constraint(lprec=sim_mod, xt=c(info_list[[j]][, i], -1), type="<", rhs=0, indices=indices)
	# 			}
	# 			k <- k + 1
	# 		}
	# 	}
	# }
	
	
	# constraint: minimum module length
	for(f in 1:F) {
		add.constraint(lprec=sim_mod, xt=rep(1, I), type=">=", rhs=min.module, indices=(I*(f-1)+1):(I*f))
	}

	# constraint: target information (relative target)
	for(w in 1:n.pathway) {
		for(g in 1:n.group) {
			if(w %in% path.group[[g]]) {
				for(i in 1:length(theta_list[[g]])) {
					temp <- cbind((pathway[w, ]*I-I + 1), (pathway[w, ]*I))
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


	# constraint: target information (relative target)
	# for(w in 1:n.pathway) {
	# 	for(g in 1:n.group) {
	# 		if(w %in% path.group[[g]]) {
	# 			for(i in 1:length(theta_list[[g]])) {
	# 				temp <- cbind((pathway[w, ]*I-I + 1), (pathway[w, ]*I))
	# 				indices <- NULL
	# 				for(s in 1:n.stage) {
	# 					indices <- c(indices, temp[s, 1]:temp[s, 2])
	# 				}
	# 				indices <- c(indices, M)
	# 				add.constraint(lprec=sim_mod, xt=c(rep(info_list[[g]][, i], n.stage), -1), type=">=", rhs=0, indices=indices)
	# 			}
	# 		}
	# 	}
	# }

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
	rownames(select_df) <- c(paste0("item.", 1:I), "Total")
	colnames(select_df) <- c("router", paste0("stage2_m", 1:nmod[1]), paste0("stage3_m", 1:nmod[2]), "sum")
	
	df_sim_list <- vector('list', F)
	names(df_sim_list) <- c("router", paste0("stage2_m", 1:nmod[1]), paste0("stage3_m", 1:nmod[2]))
	class_list <- vector('list', F)
	for(f in 1:F) {
		df_sim_list[[f]] <- df_bank[sim_opt_list[[f]] == 1, ]
		x <- item.pool$CLASS[sim_opt_list[[f]] == 1]
		df_sim_list[[f]]$CLASS <- x
	}

	# test information for each module
	info_rt <- test.info(x=df_sim_list[[1]], theta=theta, D=D)$testInfo
	s=1
	for(i in (s+1):(s+nmod[1])) {
		assign(paste0("info_stg2_", i-1), test.info(x=df_sim_list[[i]], theta=theta, D=D)$testInfo)
	}
	s=4
	for(i in (s+1):(s+nmod[2])) {
		assign(paste0("info_stg3_", i-nmod[2]-1), test.info(x=df_sim_list[[i]], theta=theta, D=D)$testInfo)
	}


myinfo.1 <- test.info(x=df_sim_list[[5]], theta=theta, D=D)$testInfo
myinfo.2 <- test.info(x=df_sim_list[[6]], theta=theta, D=D)$testInfo
myinfo.3 <- test.info(x=df_sim_list[[7]], theta=theta, D=D)$testInfo

plot(myinfo.1 ~ theta, ylim=c(0, 15), type='l', col=1)
lines(myinfo.2 ~ theta, col=2)
lines(myinfo.3 ~ theta, col=3)

info_path1 <- info_rt + info_stg2_1 + info_stg3_1
info_path2 <- info_rt + info_stg2_1 + info_stg3_2
info_path3 <- info_rt + info_stg2_2 + info_stg3_1
info_path4 <- info_rt + info_stg2_2 + info_stg3_2
info_path5 <- info_rt + info_stg2_2 + info_stg3_3
info_path6 <- info_rt + info_stg2_3 + info_stg3_2
info_path7 <- info_rt + info_stg2_3 + info_stg3_3


plot(info_rt ~ theta, ylim=c(0, 25), type='l', col=1)
lines(info_stg2_1 ~ theta, col=2, lty=1)
lines(info_stg2_2 ~ theta, col=2, lty=2)
lines(info_stg2_3 ~ theta, col=2, lty=3)
lines(info_stg3_1 ~ theta, col=3, lty=1)
lines(info_stg3_2 ~ theta, col=3, lty=2)
lines(info_stg3_3 ~ theta, col=3, lty=3)


plot(info_stg3_1 ~ theta, ylim=c(0, 10), type='l', col=3, lty=1)
lines(info_stg3_2 ~ theta, col=3, lty=2)
lines(info_stg3_3 ~ theta, col=3, lty=3)


plot(info_path1 ~ theta, ylim=c(0, 28), type='l', col=1, lwd=3, ylab="TIF")
lines(info_path2 ~ theta, col=2, lwd=2)
lines(info_path3 ~ theta, col=3, lwd=2)
lines(info_path4 ~ theta, col=4, lwd=3)
lines(info_path5 ~ theta, col=5, lwd=2)
lines(info_path6 ~ theta, col=6, lwd=2)
lines(info_path7 ~ theta, col=7, lwd=3)
abline(v=-.44)
abline(v=.44)
legend("topright", legend=c("1M-2E-3E", "1M-2M-3M", "1M-2H-3H"), col=c(1, 4, 7), lwd=3,
		title="Primary Pathway")



	info_list <- vector('list', F)
	info_list[[1]] <- info_rt
	for(i in 1:nmod[1]) {
		info_list[[1+i]] <- get(paste0("info_stg2_", i))
	}
	for(i in 1:nmod[2]) {
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
	
	# a type of target distribution
	target.type <- ifelse(relative == TRUE, "relative", "absolute")
	
	rr <- list(var.opt=sim_opt, itemvar.list=sim_opt_list, itemvar.df=select_df, 
			prm.df=df_sim_list, info.obs=info.obs, target.sc=list(router=target_rt_sc, stage2=target_stg2_sc), 
			theta=theta, ssd=ssd, solution=solution, metainfo=list(target.type=target.type, n.module=nmod, n.stage=nstg, size.module=m.size)
			)

	class(rr) <- "ata_mst"
	rr

}

