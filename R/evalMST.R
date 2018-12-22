# This function analytically computes measurement precision and three objective functions 
anal_evalMST <- function(x, route.method=c("AMI", "DPI"), theta=seq(-4, 4, 0.1), range.theta=c(-5, 5), 
                         interval=c(-2.5, 2.5), D=1) {
   
   # Arg
   # x: An object of class "atamst"
   # RDP: A numeric vector of routing decision points
   # route.method: A character string specifying the routing decision method
   # theta: A numeric vector of abiltiy points where conditional bias and CSEE are computed
   # range.theta: A numeric vector of two values specifying the lower and upper bounds of ability estimates
   
   
   # return NA when the ATA result of MST is NULL
   if(is.null(x)) {
      
      # return results
      res <- list(eos_list=NA, cond_moments=NA, object.fn=list(mrel=NA, ave_csee=NA, max_csee=NA))
      return(res)
      
   }
   
   # read required information from the ATA results
   ata_forms <- x$prm.df
   post <- x$metainfo$post
   route.map <- x$metainfo$route.map
   
   # find information from a routing map
   panel.info <- panel_info(route.map)
   config.info <- panel.info$config
   n.module <- panel.info$n.module
   pathway <- panel.info$pathway
   nstg <- panel.info$n.stage
   
   # make the RDP list
   RDP <- post[c(-1, -length(post))] 
   if(nstg == 2) RDPList <- list(NA, RDP)
   if(nstg == 3) RDPList <- list(NA, RDP, RDP)
   
   # estimate the observed equated scores across all (sub) pathways
   eos_list <- est_eos(ata_forms, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)
   
   # conditional raw score distribution of each module given an ability value
   cond_dist <- ability_dist_sep(ata_forms, theta, route.map, D)
   
   # conditional ability distribution of total-test ability estimates given a true ability
   eos.path <- eos_list$eos_path
   n.path <- eos_list$n_path
   
   # conditional joint distribution
   if(route.method == "AMI") {
      # find points where two adjacent TIFs intersect across all stages
      cut.score <- cutoff(ata_forms, route.map, RDP=RDPList, D, range.theta, interval)
      
      # comute the joint conditional distribution
      joint.dist <- ability_dist(x=cond_dist, eos.path, n.path, cut.score=cut.score)  
   }
   if(route.method == "DPI") {
      # comute the joint conditional distribution
      joint.dist <- ability_dist(x=cond_dist, eos.path, n.path, cut.score=RDPList)  
   }
   
   # calculate a mean and sd of conditional ability distribution at each ability value
   if(nstg == 2) {
      cond_moments <- sapply(joint.dist$joint.dist$stage2, cal_moments, node=eos_list$eos_path$stage2)
   }
   if(nstg == 3) {
      cond_moments <- sapply(joint.dist$joint.dist$stage3, cal_moments, node=eos_list$eos_path$stage3)
   }
   
   # Compute three objective functions
   # generate weights
   w <- dnorm(theta, 0, 1)
   w <- w/sum(w)
   
   # compute objective function values
   var.cond <- cond_moments[2, ]
   mrel <- objfn(obj="mrel", var.cond, var.pop=1, w=w, range=c(-2, 2))
   ave_csee <- objfn(obj="ave.se", var.cond, var.pop=1, w=w, range=c(-2, 2))
   max_csee <- objfn(obj="max.se", var.cond, var.pop=1, w=w, range=c(-2, 2))
   
   # return results
   res <- list(eos_list=eos_list, cond_moments=cond_moments, 
               object.fn=list(mrel=mrel, ave_csee=ave_csee, max_csee=max_csee))
   return(res)
   
}


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

