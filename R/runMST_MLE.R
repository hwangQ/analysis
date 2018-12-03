# a function to administer MST for individual examinee using ML estimation
runMST_indvMLE <- function(theta, itemBank, module.mat, ini.module=1, route.map, cut.mat=NULL, D=1.702, route.scoring="EAP", 
                           final.scoring="ML", range=c(-4, 4)){
   
   # use true theta values
   true.theta <- theta
   
   # set argumetns
   start <- list(fixModule=ini.module, D=D)
   if(is.null(cut.mat)) {
      test <- list(method=route.scoring, D=D, moduleSelect = "MFI")
   } else {
      test <- list(method=route.scoring, D=D, cutoff = cut.mat)    
   }
   final <- list(method=final.scoring, range=range, D=D)
   
   # implement MST
   admin <- mstR::randomMST(trueTheta=true.theta, itemBank=itemBank, modules=module.mat, 
                            transMatrix=route.map, start=start, test=test, final=final)
   
   # extract results 
   th <- admin$thFinal
   modules <- admin$selected.modules
   responses <- admin$pattern
   test.items <- admin$testItems
   
   # make a list of results 
   rr <- list(est.theta=th, true.theta=true.theta, modules=modules, responses=responses)
   
   rr
   
}

# a function to administer MST for all examinees using ML estimation

runMST_MLE <- function(theta, item.pool, module.mat, ini.module=1, route.map, cut.mat=NULL, D=1.702, route.scoring="EAP", 
                        final.scoring="ML", range=c(-4, 4), full.out=TRUE) {

   # make sure that true thetas are in a vectoc
   true.theta <- as.numeric(theta)

   # create an item bank
   itemBank <- cbind(item.pool[,c("PARAM.1", "PARAM.2", "PARAM.3")], d=rep(1, nrow(item.pool)))
   colnames(itemBank) <- c("a", "b", "c", "d")
   
   # the number of examinees
   nstd <- length(true.theta)

   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway
   n.stage <- panel.info$n.stage
 
   # make empty cells to contain results
   nitem <- sum(module.mat[, pathway[1, ]])
   est.theta <- rep(NA, nstd)
   modules <- array(NA, c(nstd, n.stage))
   colnames(modules) <- paste0("stage", 1:n.stage)
   responses <- array(NA, c(nstd, nitem))
   colnames(responses) <- paste0("item", 1:nitem)
   
   # run MAPT for all examinees
   for(i in 1:nstd) {
       temp <- runMST_indvMLE(true.theta[i], itemBank, module.mat, ini.module=ini.module, route.map, cut.mat, D, 
                              route.scoring, final.scoring, range)
       est.theta[i] <- temp$est.theta
       modules[i, ] <- temp$modules
       responses[i, ] <- temp$responses
   }

   # summary
   score.table <- data.frame(true.theta=theta, est.theta=est.theta)
   if(full.out) {
      rr <- list(score.table=score.table, modules=modules, responses=responses)
   } else {
      rr <- list(score.table=score.table)
   }
   
   rr
 
}

# a function to administer MST for replications using ML estimation
runMST_MLErep <- function(nrep, theta, item.pool, module.mat, ini.module=1, route.map, cut.mat=NULL, D, route.scoring="EAP", 
                          final.scoring="ML", range=c(-4, 4), full.out=FALSE, iseed=0L, n.core=3) {

   # set the number of cpu cores
   if(n.core > 3) stop("'n.core' must be less than 4.", call.=FALSE)
   left.core <- 4 - n.core
   numCores <- parallel::detectCores() - left.core

   # create a parallel processesing cluster
   cl = parallel::makeCluster(numCores, type="PSOCK")
   
   # load mstR package
   parallel::clusterEvalQ(cl, library(mstR))
   
   # load some specific variable names into processing cluster
   parallel::clusterExport(cl, c("theta", "item.pool", "module.mat", "ini.module", "route.map", "cut.mat", "D", "route.scoring", 
                                 "final.scoring", "range", "full.out"), envir = environment())
   
   # load source files used for the replications into processing cluster
   lapply(1:length(source.files), function(i) parallel::clusterCall(cl, fun=source, source.files[i]))
   
   # run MAPT for all replications
   f <- function(i) {
      runMST_MLE(theta, item.pool, module.mat, ini.module, route.map, cut.mat, D, route.scoring, final.scoring, range, full.out)
   }
   admin <- pbapply::pblapply(X=1:nrep, FUN=f, cl=cl) # to see the progress bar
   
   # finish
   parallel::stopCluster(cl)
   
   # reorganize the results
   true.theta <- theta
   est.theta <- sapply(admin, '[[', c(1, 2))
   colnames(est.theta) <- paste0("rep.", 1:nrep)
   if(full.out) {
      modules <- lapply(admin, '[[', 2)
      names(modules) <- paste0("rep.", 1:nrep)
      responses <- lapply(admin, '[[', 3)
      names(responses) <- paste0("rep.", 1:nrep)
   }
   
   # return results
   if(full.out) {
      rr <- list(true.theta=true.theta, est.theta=est.theta, modules=modules, responses=responses)
   } else {
      rr <- list(true.theta=true.theta, est.theta=est.theta)
   }
   
   rr
   
}
