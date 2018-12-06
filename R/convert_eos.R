# a function to estimate the EOS corresponding to observed raw score for a list of data

est_eos <- function(x, pathway, ...) {
   
   ## Arguement
   # x: a list of item paramameter data.frame for all modules across all stages
   # pathway: a maxtrix including all pathways
   
   # check the number of stages
   n.stage <- ncol(pathway)
   
   # create a list of item parameter data.frame of for all possible (sub) paths across all stages
   df_path <- vector('list', n.stage)
   names(df_path) <- paste0("stage", 1:n.stage)
   n_path <- NA # the number of paths at each stage
   for(i in 1:n.stage) {
       path.temp <- unique(as.matrix(pathway[, 1:i]))
       n.temp <- nrow(path.temp)
       df_path[[i]] <- lapply(1:n.temp, function(j) do.call('rbind', x[path.temp[j, ]]))
       n_path[i] <- n.temp
   }
   
   # find observed equated scores using TCC method corresponding to each path at each stage
   eos_path <- vector('list', n.stage)
   names(eos_path) <- paste0("stage", 1:n.stage)
   for(i in 1:n.stage) {
      eos_path[[i]] <- sapply(1:length(df_path[[i]]), function(j) convert_eos(x=df_path[[i]][[j]], ...)$theta)
      rownames(eos_path[[i]]) <- 0:(nrow(eos_path[[i]]) - 1)
      colnames(eos_path[[i]]) <- paste0("path.", 1:ncol(eos_path[[i]]))
   }

   rr <- list(eos_path=eos_path, df_path=df_path, n_path=n_path)
   rr

}


# a function to convert raw scores to equated observed scores

convert_eos <- function(x, D=1, range.theta=c(-5, 5), intpo=TRUE, constraint=TRUE, lower=-10, upper=10) { 
  
  prm_list <- paramList(x)
  
  ##########################################
  ## Equip all information
  any.dc <- any(length(prm_list$dicho.item$id) > 0) # if there are dichotomous items
  any.py <- any(length(prm_list$poly.item$id) > 0) # if there are polytomous items
  
  # Find a maximum raw score
  max.score <- 0
  if(any.dc) {
    max.score <- sum(max.score, length(prm_list$dicho.item$a))
  }
  if(any.py) {
    sum.score <- prm_list$poly.item$score.cats - 1
    max.score <- sum(max.score, sum.score)
  }
  
  # Find a minimum raw score
  if(any.dc) {
    g <- prm_list$dicho.item$g
    min.g <- sum(g)
    min.score <- ifelse(min.g != 0, ceiling(min.g), 1)
  } else {
    min.score <- 1
  }
  
  ##########################################
  ## Condunct 3 steps of TSE 
  # 1) Find the range of raw scores
  scores <- seq(min.score, max.score, 1)
  
  # 2) Find the theta values correspoding to the raw scores 
  theta <- sapply(scores, raw2theta, df.scaled=x, D=D, lower=lower, upper=upper)
  if(constraint) {
    theta <- ifelse(theta < range.theta[1], range.theta[1], theta)
    theta <- ifelse(theta > range.theta[2], range.theta[2], theta)
  }
  
  # Create conversion table
  if(intpo) {
    
    if(min.score > 1) {
      lessg.score <- 0:(min.score-1)
      
      # calculate slope and intercept for linear interpolation
      if(range.theta[1] > theta[1]) stop(paste0("A lower limit of theta is greater than a minimum theta value of ", 
                                                round(theta[1], 3), " corresponding to the minimum raw score. Use the different a range of theta."))
      if(range.theta[2] < tail(theta, 1)) stop(paste0("A upper limit of theta is less than a maximum theta value of ", 
                                                      round(tail(theta, 1), 3), " corresponding to the maximum raw score. Use the different a range of theta."))
      slope <- scores[1] / (theta[1] - range.theta[1])
      intercept <- - slope * range.theta[1]
      lessg.theta <- (lessg.score - intercept) / slope
      
      theta <- c(lessg.theta, theta)
      theta <- ifelse(is.nan(theta), range.theta[1], theta)
      scale.table <- data.frame(score=0:max.score, theta=theta)
    } else {
      if(range.theta[1] > theta[1]) stop(paste0("A lower limit of theta is greater than a minimum theta value of ", 
                                                round(theta[1], 3), " corresponding to the minimum raw score. Use the different a range of theta."))
      if(range.theta[2] < tail(theta, 1)) stop(paste0("A upper limit of theta is less than a maximum theta value of ", 
                                                      round(tail(theta, 1), 3), " corresponding to the maximum raw score. Use the different a range of theta."))
      theta <- c(range.theta[1], theta)
      theta <- ifelse(is.nan(theta), range.theta[1], theta)
      scale.table <- data.frame(score=0:max.score, theta=theta)
    }
    
  } else {
    
    if(min.score_new > 1) {
      lessg.score <- 0:(min.score-1)
      lessg.theta <- rep(NA, length(lessg.score))
      scale.table <- data.frame(score=0:max.score, theta=c(lessg.theta, theta))
    } else {
      scale.table <- data.frame(score=0:max.score, theta=c(NA, theta))
    }
    
  }
  
  scale.table
  
}


# "raw2theta" function 
# This is a function to find theta corresponding to raw score.
# This function is used inside "convert_eos" function
raw2theta <- function(rawScore, df.scaled, D = 1, value = .5, lower = -10, upper = 10){

	params <- paramList(df.scaled)

	any.dc <- any(length(params$dicho.item$id) > 0) # if there are dichotomous items
	any.py <- any(length(params$poly.item$id) > 0) # if there are polytomous items

	# Find a maximum raw score 
	max.score <- 0
	if(any.dc) {
		max.score <- sum(max.score, length(params$dicho.item$a))
	}
	if(any.py) {
		sum.score <- params$poly.item$score.cats - 1
		max.score <- sum(max.score, sum.score)
	}

	# Adjust perfect and zero raw scores with a certain value 
	if(rawScore == max.score) rawScore <- rawScore - value
	if(rawScore == 0) rawScore <- rawScore + value

	# Extract item parameters
	if(any.dc) {
		a <- params$dicho.item$a
		b <- params$dicho.item$b
		g <- params$dicho.item$g
		min.g <- sum(g)
	} else {
		a <- b <- g <- NULL
	}
	if(any.py){
		aa <- params$poly.item$a
		d <- params$poly.item$d
		pModel <- params$poly.item$pModel
		Npoly <- length(d)
	} else {
		aa <- d <- pModel <- NULL
	}

	# Set low raw scores greater than sum of guessing parameters 
	# to something that can be estimated
	if(rawScore <= min.g) rawScore <- ceiling(min.g)

	# Set a loss function to estimate theta
	f1 <- function(theta) sum(pl.fn(theta, a, b, g, D)) - rawScore
	f2 <- function(theta) sum(mapply(poly.exp, theta, aa, d, D, pModel)) - rawScore
	f3 <- function(theta) sum(sum(pl.fn(theta, a, b, g, D)),
							sum(mapply(poly.exp, theta, aa, d, D, pModel))) - rawScore

	# Estimate the ability corresponding to a raw score using "uniroot" function 
	if(!any.py) root <- uniroot(f1, c(lower,upper), extendInt="yes", maxiter = 5000)$root #only binary items
	if(!any.dc) root <- uniroot(f2, c(lower,upper), extendInt="yes", maxiter = 5000)$root #only poly items
	if(any.dc & any.py) root <- uniroot(f3, c(lower,upper),extendInt="yes", maxiter = 5000)$root # mixed-format
	
	root
}