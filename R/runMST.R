# a function to administer MSTS for an individual examinee
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
   true.theta <- theta
   for(s in 1:n.stage) {
      
       # find item parameter informaion corresponding to the selected module
       params <- module.form[[sel.mod]]
       a <- params$PARAM.1
       b <- params$PARAM.2
       g <- params$PARAM.3
       cats <- params$CATEGORY

       # simulate item responses and calculate a summed score for the selected module 
       responses <- c(responses, simdat(theta=true.theta, a.dc=a, b.dc=b, g.dc=g, cats=cats, D=D))
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
   rr <- list(est.theta=th, true.theta=true.theta, sum.score=sum.score, modules=modules, responses=responses, path.form=df_params)
   rr
   
}

# a function to administer MST for all examinees
runMST <- function(theta, module.form, route.map, cut.score, eos.path, D, range.theta, full.out=TRUE) {

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
   if(missing(eos.path)) {
      eos.path <- est_eos(module.form, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)$eos_path
   }
   
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
   if(full.out) {
      rr <- list(score.table=score.table, modules=modules, responses=responses)
   } else {
      rr <- list(score.table=score.table)
   }

   rr
   
}


# a function to administer MST for replications
runMST_rep <- function(nrep, theta, module.form, route.map, cut.score, D, range.theta, source.files, full.out=TRUE, iseed=0L) {

   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway

   # estimate the equated observed score for all pathways across all stages
   eos.path <- est_eos(module.form, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)$eos_path

   # set the number of cpu cores to n - 1 cores.
   numCores <- parallel::detectCores() - 1

   # create a parallel processesing cluster
   cl = parallel::makeCluster(numCores, type="PSOCK")

   # load some specific variable names into processing cluster
   parallel::clusterExport(cl, c("theta", "module.form", "route.map", "cut.score", 
                                 "eos.path", "D", "range.theta", "full.out"), envir = environment())
   
   # load source files used for the replications into processing cluster
   lapply(1:length(source.files), function(i) parallel::clusterCall(cl, fun=source, source.files[i]))

   parallel::clusterSetRNGStream(cl, iseed=iseed)
   # run MST for all replications
   f <- function(i) runMST(theta, module.form, route.map, cut.score, eos.path, D, range.theta, full.out)
   # admin <- parallel::parLapplyLB(cl, 1:nrep, f)
   admin <- pbapply::pblapply(X=1:nrep, FUN=f, cl=cl) # to see the progress bar
   
   # finish
   parallel::stopCluster(cl)
   
   # reorganize the results
   true.theta <- theta
   est.theta <- sapply(admin, '[[', c(1, 2))
   colnames(est.theta) <- paste0("rep.", 1:nrep)
   sum.score <- sapply(admin, '[[', c(1, 3))
   colnames(sum.score) <- paste0("rep.", 1:nrep)
   if(full.out) {
      modules <- lapply(admin, '[[', 2)
      names(modules) <- paste0("rep.", 1:nrep)
      responses <- lapply(admin, '[[', 3)
      names(responses) <- paste0("rep.", 1:nrep)
   }
   
   # return results
   if(full.out) {
      rr <- list(true.theta=true.theta, est.theta=est.theta, sum.score=sum.score, modules=modules, responses=responses)
   } else {
      rr <- list(true.theta=true.theta, est.theta=est.theta, sum.score=sum.score)
   }
   
   rr

}



