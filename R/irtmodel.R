# Functions for IRT models 

# IRT 1, 2, 3PL models
pl.fn <- function(theta, a, b, g, D) {
   
   # check the numbers of examinees and items
   nstd <- length(theta)
   nitem <- length(a)
   
   # when the numbers of examiness and items are greater than 1
   if(nstd > 1 & nitem > 1) {
      a <- matrix(a, nrow=nstd, ncol=nitem, byrow=TRUE)
      b <- matrix(b, nrow=nstd, ncol=nitem, byrow=TRUE) 
      g <- matrix(g, nrow=nstd, ncol=nitem, byrow=TRUE) 
   }
   
   # calculate probability of correct answer
   Da <- D * a
   z <- Da * (theta - b)
   P <- g + (1 - g) / (1 + exp(-z))
   
   P
 
}


# IRT GPC model
gpcm.fn <- function (theta, score, a, steps, D) {

   # check the numbers of examinees and items
   nstd <- length(theta)
   
   # when the number of examiness is greater than 1
   if(nstd > 1) {
      a <- matrix(a, nrow=nstd, ncol=1, byrow=TRUE) 
      steps <- matrix(steps, nrow=nstd, ncol=length(steps), byrow=TRUE) 

      # calculate probability of correct answer
      Da <- as.numeric(D * a)
      P <- exp(rowSums(Da * (theta - as.matrix(steps[, 1:(score+1)]))))/ 
               rowSums(exp(t(apply(Da * (theta - steps), 1, cumsum))))
   
   } else {
   
    # calculate probability of correct answer
     Da <- D * a
     P <- exp(sum(Da * (theta - steps[1:(score+1)])))/ 
             sum(exp(cumsum(Da * (theta - steps))))
   }
   
   P

}


# IRT GRM model
grm.fn <- function(theta, score, a, steps, D) {

   # check the number of step parameters
   m <- length(steps)
   
   # check the numbers of examinees and items
   nstd <- length(theta)
   
   # calculate all the probabilities greater than equal to each threshold
   allP <- sapply(1:m, function(i) pl.fn(theta, a, b=steps[i], g=0, D))
   
   # all possible scores
   allScores <- seq(0, m)

   # when the number of examiness is greater than 1
   if(nstd > 1) { 
      if(score == min(allScores)) {
         P <- 1 - allP[, 1]
      }
      if(score > min(allScores) & score < max(allScores)) {
         P <- allP[, score] - allP[, score+1] 
      }
      if(score == max(allScores)) {
         P <- allP[, score]	 
      }
   } else {
      if(score == min(allScores)) {
         P <- 1 - allP[1]
      }
      if(score > min(allScores) & score < max(allScores)) {
         P <- allP[score] - allP[score+1] 
      }
      if(score == max(allScores)) {
         P <- allP[score]
      }
   }
   
   P

}


# Trace Line for polytomous IRT model
poly.fn <- function(theta, a, d, D, pModel) {
	if(pModel == "GPCM") {
		P.vec <- sapply(1:length(d), function(i) gpcm.fn(theta, score=(i-1), a, steps=d, D))
	}
	if(pModel == "GRM") {
		P.vec <- sapply(1:(length(d)+1), function(i) grm.fn(theta, score=(i-1), a, steps=d, D))
	}

	return(P.vec)
}

# TCC for polytomous IRT models
poly.exp <- function(theta, a, d, D, pModel) {
	if(pModel == "GPCM") {
		sum.w <- sum(sapply(1:length(d), function(i) gpcm.fn(theta, score=(i-1), a, steps=d, D) * (i-1)))
	}
	if(pModel == "GRM") {
		sum.w <- sum(sapply(1:(length(d)+1), function(i) grm.fn(theta, score=(i-1), a, steps=d, D) * (i-1)))
	}
	return(sum.w)
}
