######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 09/27/2018
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
ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
         lp.control=list(timeout=30, epsint=0.2, mip.gap=c(0.1, 0.05)))
 
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
   cont <- item.pool$CLASS_RE

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
   test.length <- constraints$test.length
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
       add.constraint(lprec=sim_mod, xt=rep(1, I), type="=", rhs=m.size[nmod[idx]], indices=(I*(f-1)+1):(I*f))
   }
   
   # constraint: content category # from here 
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
 
 
   
   
   

theta <- seq(-4, 4, 0.1)
	# test information for each module
	info_rt <- test.info(x=df_sim_list[[1]], theta=theta, D=D)$testInfo
	s=1
	for(i in (s+1):(s+nmod[2])) {
		assign(paste0("info_stg2_", i-1), test.info(x=df_sim_list[[i]], theta=theta, D=D)$testInfo)
	}
	s=4
	for(i in (s+1):(s+nmod[3])) {
		assign(paste0("info_stg3_", i-nmod[3]-1), test.info(x=df_sim_list[[i]], theta=theta, D=D)$testInfo)
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


plot(info_path1 ~ theta, ylim=c(0, 50), type='l', col=1, lwd=3, ylab="TIF")
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


















