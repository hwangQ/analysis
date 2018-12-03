library(parallel)

# set directory
mydir <- "E:/Dissertation/Analysis/Temp"
setwd(mydir)

# Read source code
src.dir <- "E:/Dissertation/Analysis/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set condtions
theta <- seq(-4, 4, 0.1)
range.theta <- c(-4, 4)
interval <- c(-2.5, 2.5)
n.module <- c(3, 3)
n.stage <- 3
D <- 1.702

# read a routing map
route.map <- read.csv("Routing_map.csv", header=FALSE)
panel.info <- panel_info(route.map)
config.info <- panel.info$config

# read the assembled test forms
ata_forms <- df_sim_list

# create a list of item parametr data.frames across all modules at all stages
df_forms <- vector('list', n.stage)
names(df_forms) <- paste0("stage", 1:n.stage)
module.label <- vector('list', n.stage)
for(i in 1:n.stage) module.label[[i]] <- unique(pathway[, i])
for(i in 1:n.stage) {
	df_forms[[i]] <- vector('list', length(module.label[[i]]))
	names(df_forms[[i]]) <- paste0("m", 1:length(module.label[[i]]))
	for(j in 1:length(df_forms[[i]])) {
		df_forms[[i]][[j]] <- ata_forms[[module.label[[i]][[j]]]]
	}
}

# estimate the observed equated scores across all (sub) pathways
x <- ata_forms
eos_list <- est_eos(x, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)

# find points where two adjacent TIFs intersect across all stages
module.form <- ata_forms
RDP <- vector('list', 3)
RDP[[1]] <- NA
RDP[[2]] <- c(-.44, .44)
RDP[[3]] <- c(-.44, .44)

cutoff(module.form, route.map, RDP, D, range.theta, interval)

cutoff <- function(module.form, route.map, RDP, D, range.theta, interval) {

   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway
   config.info <- panel.info$config.info
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module

   # create an empty list to include cutscore information
   cut.score <- vector('list', n.stage)
   names(cut.score) <- paste0("stage", 1:n.stage)

   for(s in 1:n.stage) {
       idx <- config.info[[s]]
       forms <- lapply(idx, function(i) module.form[[i]])
       if(length(forms) == 1) {
          cut.score[[s]] <- NA
       } else {
          cut.score[[s]] <- cross_info(forms, RDP=RDP[[s]], range.theta=range.theta, D=D, interval=interval, n=1000)$cut.score
       }
   }

   cut.score

}

module.form <- ata_forms
theta <- seq(-4, 4, 0.1)
source.files <- paste0(src.dir, src.files)
nrep <- 1000

x <- runMST_rep(nrep, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=FALSE)

w <- dnorm(theta, 0, 1)
w <- w/sum(w)
myeval <- evalMST(x, w)

point <- myeval$true.theta
vars <- myeval$cond.eval$variance

plot(sqrt(vars) ~ point)

lines(cond_sigma ~ theta)

# a function to calculate (average) bias, (average) variance, (average) mse
evalMST <- function(x, weights=NULL) {

   ## Arguments
   # x: an object obtained from "runMST_rep" function
   # weights: a vector of weights for each ability point. Default is NUL:

   # extract true abilities and estimates
   true.theta <- x$true.theta
   est.theta <- x$est.theta
   
   # difference between true and estimates
   diff.ability  <- est.theta - true.theta
   mean.est <- apply(est.theta, 1, mean)
   
   # conditional evaluation
   cond.bias <- apply(diff.ability, 1, mean)
   cond.var <- apply((est.theta - mean.est)^2, 1, mean)
   cond.mse <- cond.bias^2 + cond.var
   cond.eval <- data.frame(bias=cond.bias, variance=cond.var, mes=cond.mse)
   
   # average evaluation
   if(is.null(weights)) {
      mean.bias <- mean(cond.bias)
      mean.var <- mean(cond.var)
      mean.mse <- mean(cond.mse)
   } else {
      mean.bias <- weighted.mean(x=cond.bias, w=weights)
      mean.var <- weighted.mean(x=cond.var, w=weights)
      mean.mse <- weighted.mean(x=cond.mse, w=weights)
   }
   mean.eval <- c(bias=mean.bias, variance=mean.var, mse=mean.mse)
   
   # summary
   if(is.null(weights)) weights <- NA
   
   rr <- list(true.theta=true.theta, cond.eval=cond.eval, mean.eval=mean.eval, weights=weights)
   rr
 
}


# a function to administer MST for replications
runMST_rep <- function(nrep, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=TRUE) {

   # set the number of cpu cores to n - 1 cores.
   numCores <- detectCores() - 1

   # create a parallel processesing cluster
   cl = makeCluster(numCores, type="PSOCK")

   # load some specific variable names into processing cluster
   clusterExport(cl, c("theta", "module.form", "route.map", "cut.score", "D", "range.theta"))
   
   # load source files used for the replications into processing cluster
   lapply(1:length(source.files), function(i) clusterCall(cl, fun=source, source.files[i]))

   # run MST for all replications
   f <- function(i) runMST(theta, module.form, route.map, cut.score, D, range.theta)
   admin <- parLapplyLB(cl, 1:nrep, f)

   # finish
   stopCluster(cl)
   
   # reorganize the results
   true.theta <- theta
   est.theta <- sapply(admin, '[[', c(1, 2))
   colnames(est.theta) <- paste0("rep.", 1:nrep)
   sum.score <- sapply(admin, '[[', c(1, 3))
   colnames(sum.score) <- paste0("rep.", 1:nrep)
   modules <- lapply(admin, '[[', 2)
   names(modules) <- paste0("rep.", 1:nrep)
   responses <- lapply(admin, '[[', 3)
   names(responses) <- paste0("rep.", 1:nrep)

   if(full.out) {
      rr <- list(true.theta=true.theta, est.theta=est.theta, sum.score=sum.score, modules=modules, responses=responses)
   } else {
      rr <- list(true.theta=true.theta, est.theta=est.theta, sum.score=sum.score, modules=modules)
   }
   
   rr

}


runMST(theta, module.form, route.map, cut.score, D, range.theta)

# a function to administer MAPT for all examinees
runMST <- function(theta, module.form, route.map, cut.score, D, range.theta) {

   # Arguments
   # theta: a vector of true ability values
   # module.form: a list including a data.frame of items of modules across all stages
   # route.map: a matrix of the routing information
   # cut.score: a list of cut scores across all stages
   
   # checek the number of examinees
   theta <- as.numeric(theta)
   nstd <- length(theta)
   
   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # estimate the equated observed score for all pathways across all stages
   eos.path <- est_eos(module.form, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)$eos_path
   
   # make empty cells to contain results
   nitem <- nrow(eos.path[[length(eos.path)]]) - 1 
   est.theta <- rep(NA, nstd)
   sum.score <- rep(NA, nstd)
   modules <- array(NA, c(nstd, n.stage))
   colnames(modules) <- paste0("stage", 1:n.stage)
   responses <- array(NA, c(nstd, nitem))
   colnames(responses) <- paste0("item", 1:nitem)
   
   # run MST for all examinees
   for(i in 1:nstd) {
       temp <- runMST_indv(theta[i], module.form, panel.info, cut.score, eos.path, ini.module=1, D, range.theta)
       est.theta[i] <- temp$est.theta
       sum.score[i] <- temp$sum.score
       modules[i, ] <- temp$modules
       responses[i, ] <- temp$responses
   }

   # summary
   score.table <- data.frame(true.theta=theta, est.theta=est.theta, sum.score=sum.score)
   rr <- list(score.table=score.table, modules=modules, responses=responses)
   rr
   
}



eos.path <- eos_list$eos_path
panel.info <- panel_info(route.map)
module.form <- ata_forms
theta <- 2

# without eos information 
runMST_indv(theta=theta, module.form=module.form, panel.info=panel.info, cut.score=cut.score, ini.module=1, D=D, range.theta=range.theta)

# with eos information
runMST_indv(theta, module.form, panel.info, cut.score, eos.path, ini.module=1, D, range.theta)

runMST_indv <- function(theta, module.form, panel.info, cut.score, eos.path, ini.module=1, D, range.theta, ...) {

   # Arguments
   # theta: a scala of true ability value
   # module.form: a list including a data.frame of items of modules across all stages
   # panel.info: an object obtained from "panel_info" function. This includes all informations regarding a panel configuration
   # cut.score: a list of cut scores across all stages
   # eos.path: a list of the equated observed (or number correct) scores for every paths across all stages
   # ini.module: a scala indicating a module to be chosen at the first stage.

   responses <- NULL
   modules <- ini.module
   df_params <- data.frame()
   
   # extract panel information
   pathway <- panel.info$pathway
   config.info <- panel.info$config
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # simulate item responses and estimate an examinee ability
   sel.mod <- ini.module
   th <- theta
   for(s in 1:n.stage) {
      
       # find item parameter informaion corresponding to the selected module
       params <- module.form[[sel.mod]]
       a <- params$PARAM.1
       b <- params$PARAM.2
       g <- params$PARAM.3
       cats <- params$CATEGORY

       # simulate item responses and calculate a summed score for the selected module 
       responses <- c(responses, simdat(theta=th, a.dc=a, b.dc=b, g.dc=g, cats=cats, D=D))
       sum.score <- sum(responses)
      
       # estimate ability parameters using the inverse TCC method 
       df_params <- rbind(df_params, params)
       
       if(missing(eos.path)) {
          thetas <- convert_eos(df_params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
       } else {
          unique.path <- data.matrix(unique(pathway[, 1:s]))
          idx.logic <- sapply(1:nrow(unique.path), function(k) all(modules == unique.path[k, ]))
          thetas <- eos.path[[s]][, idx.logic]
       }
       th <- thetas[sum.score + 1]

       # find a next module to route an examinee
       if(s < n.stage) {
          selected <- give_path(th, cut.score[[s+1]])$path
          sel.mod <- config.info[[s+1]][selected]
          modules <- as.numeric(c(modules, sel.mod))
       }
    
   }
   
   # get results
   rr <- list(est.theta=th, true.theta=theta, sum.score=sum.score, modules=modules, responses=responses, path.form=df_params)
   rr
   
}





