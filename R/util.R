# This function creates a matrix of cutoff scores used for routing
cutoff_mat <- function(cut.score, n.stage, n.module) {
  
  cut.mat <- matrix(NA, 1, 2)
  for(s in 2:n.stage) {
    cut.tmp <- c(-Inf, cut.score[[s]], Inf)
    tmp <- matrix(NA, n.module[s], 2)
    for(i in 1:n.module[s]) {
      tmp[i, ] <- cut.tmp[i:(i + 1)]
    }
    cut.mat <- rbind(cut.mat, tmp)
  }
  
  cut.mat
  
}


# a function to create target module information function (MIF)
info_mif <- function(item.pool, target.theta, route.map, test.length, D=1) {

   # extract needed information from a routing map
   panel.info <- panel_info(route.map)
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # create target MIF
   target.info <- vector('list', n.stage)
   names(target.info) <- paste0("stage", 1:n.stage)
   for(s in 1:n.stage) {
       for(m in 1:n.module[s]) {
           target.info[[s]][[m]] <- rep(NA, length(target.theta[[s]][[m]]))
           for(i in 1:length(target.theta[[s]][[m]])) {
               tmp.1 <- test.info(x=item.pool, theta=target.theta[[s]][[m]][i], D=D)$itemInfo
               tmp.2 <- sort(tmp.1, decreasing=TRUE)[1:test.length]
               target.info[[s]][[m]][i] <- size.module[s] * mean(tmp.2)
           }
       }
   }

   # return results
   rr <- list(theta=target.theta, info=target.info)

   rr
   
}

# a function to assign a path to each router score
give_path <- function(score, cut.score) { 

	# number of scores
	nscore <- length(score)

	# number of categories (or levels)
	cats <- length(cut.score) + 1

	# set total lables to be given to each examinee
	level <- 1:cats

	# create path vectors
	path <- rep(level[1], nscore)

	# give a level correspoding each score
	for(i in 1:length(cut.score)) {
		path[which(score >= cut.score[i])] <- level[i+1]
	}

	rr <- list(path=path, score=score, cut.score=cut.score)
	rr

}

# calculate moments
cal_moments <- function(node, weight) {

	mu <- sum(node * weight)
	sigma2 <- sum(node^2 * weight) - mu^2
	rr <- c(mu=mu, sigma2=sigma2)
	rr

}

# marginal reliability
mrel <- function(var.pop, var.cond, wieghted=TRUE, w) {
	
	if(wieghted) {
		rr <- (var.pop - weighted.mean(x=var.cond, w=w)) / var.pop
	} else {
		rr <- (var.pop - mean(x=var.cond)) / var.pop
	}
	
	rr

}


# calculation of four moments of normal distribution
norm_moments <- function(mu, sigma) {

  m1 <- mu
  m2 <- mu^2 + sigma^2
  m3 <- mu * (mu^2 + 3 * sigma^2)
  m4 <- mu^4 + 6 * mu^2 * sigma^2 + 3 * sigma^4
  
  rr <- c(m1=m1, m2=m2, m3=m3, m4=m4)
  rr
  
}