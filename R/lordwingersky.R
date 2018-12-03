# "lwRecursive" function
# Arg:
# prob.list: (list) each element of a list includes a probability matrix (or data.frame) for an item.
#					In each probability matrix, rows indicate theta values (or weights) and columns indicate
#					score categories. Thus, each element in the matrix represents the probability that a person
#					who has a certain ability earns a certain score.
# score.cats: (vector) containes score categories for all items 

lwRecurive <- function(prob.list, score.cats) {

	if(length(unique(sapply(prob.list, nrow))) != 1L) {
		stop("At least, one probability matrix (or data.frame) has the difference number of rows across all marices in a list")
	}
	
	if(length(prob.list) != length(score.cats)) {
		stop(paste0("There are ", nrow(prob.list), " items in the probability matrix (or data.frame) at each element of a list, ", 
			"whereas there are ", length(score.cats), " items in the category vector."))
	}

	if(min(score.cats) < 2) {
		stop("Minimum number of categories for each item is 2")
	}

	# Probabilities for each category at 1st item
	p <- prob.list[[1]]
	
	# Possible observed score range for a first item
	obs.range <- 0:(score.cats[1]-1)
	
	# Create a temporary matrix to contain all probabilities 
	tScore.range <- 0:sum(score.cats-1)
	tmp <- matrix(0, nrow=nrow(prob.list[[1]]), ncol=length(tScore.range))

	# Caculate probabilities to earn an observed score by accumulating over remaining items
	for(j in 2:length(prob.list)) {

		# Probability to earn zero score. This is a special case.
		tmp[,1] <- p[,1]*prob.list[[j]][,1]
	
		# The range of category scores for an added item
		cat.range <- 0:(score.cats[j]-1)
	
		# Possible minimum and maximum observed scores when the item is added
		# but, except zero and perfect score
		min.s <- min(obs.range) + min(cat.range) + 1
		max.s <- max(obs.range) + max(cat.range) - 1

		# Calculate probability to obtain possible observed score after the item is added
		# but, still except zero and perfect score
		for(score in min.s:max.s) {

			# Difference between the possible observed score and each category score if the added item
			poss.diff <- score - min(cat.range):max(cat.range)

			# The difference above should be greater than or equal to the minimum observed score
			# where the item is added, and less than or equal to the maximum observed score where 
			# the item is added
			cols <- poss.diff >= min(obs.range) & poss.diff <= max(obs.range)
			poss.diff <- poss.diff[cols]  # Now, these are remained difference 

			# Final probability to earn the observed score
			# "p[ ,(poss.diff + 1)]" is the probabilities for possible observed score
			prob.score <- apply(p[ ,(poss.diff + 1)] * prob.list[[j]][, cols], 1, sum)
			tmp[,1+score] <- prob.score
		
		}

		
		# Probability to earn perfect score. This is a special case.
		tmp[ ,score+2] <- p[ ,ncol(p)]*prob.list[[j]][ ,score.cats[j]]
				
		# Update the range of possible observed scores 
		obs.range <- (min.s-1):(max.s+1)
		
		# Update probabilities
		p <- tmp[, 1:length(obs.range)]

	}
	
	colnames(p) <- paste0("score.", 0:sum(score.cats-1))
	rownames(p) <- paste0("theta.", 1:nrow(prob.list[[1]]))
	
	t(p)
	
}

# "prep4lw" function
# This function is for preparing data to be used in "lw4Recursive" function
# Arg:
# dframe: (data.frame) contains infomation about eath item (e.g, parameter estimates, model, score caterogy...) 
# score.cats: (vector) containes score categories for all items 

prep4lw <- function(dframe, theta, D) {

	prmList <- paramList(dframe)

	any.dc <- any(length(prmList$dicho.item$id) > 0) # if there are dichotomous items
	any.py <- any(length(prmList$poly.item$id) > 0) # if there are polytomous items

	# Extract item parameters
	if(any.dc) { # For dichotomous items
		a <- prmList$dicho.item$a
		b <- prmList$dicho.item$b
		g <- prmList$dicho.item$g
	} else {
		a <- b <- g <- NULL
	}
	if(any.py){ # For polytomous items
		aa <- prmList$poly.item$a
		d <- prmList$poly.item$d
		pModel <- prmList$poly.item$pModel
		Npoly <- length(aa)
	} else {
		aa <- d <- pModel <- NULL
	}
	
	# Create a list including probability data.frames as elements
	# Create a score category vector
	cats <- c()
	if(any.dc) { # For dichotomous items
		score.1 <- lapply(1:length(a), function(i) pl.fn(theta, a[i], b[i], g[i], D))
		score.0 <- lapply(1:length(a), function(i) 1-score.1[[i]])
		prob.dc <- lapply(1:length(a), function(i) data.frame(s0=score.0[[i]], s1=score.1[[i]]))

		cats <- c(cats, rep(2, length(a)))
	} else {
		prob.dc <- NULL
	}
	
	if(any.py) { # For polytomous items
		cats.py <- prmList$poly.item$score.cats
		prob.py <- lapply(1:length(aa), function(i) matrix(0, nrow=length(theta), ncol=cats.py[i]))
	
		for(i in 1:length(theta)) {
			prob.tmp <- mapply(poly.fn, theta[i], aa, d, D, pModel, SIMPLIFY = FALSE)
				for(j in 1:length(aa)) {
					prob.py[[j]][i,] <- prob.tmp[[j]]
				}
		}
		prob.py <- lapply(1:length(aa), function(i) data.frame(prob.py[[i]]))
		
		cats <- c(cats, cats.py)
	} else {
		prob.py <- NULL
	}
	
	# Total probability list
	prob.list <- c(prob.dc, prob.py)
	
	# Score category vector
	score.cats <- cats

	# Ordering the probability list and score category vector according to the origincal order of data
	pos <- c(prmList$dicho.item$dPos, prmList$poly.item$pPos)
	prob.list <- prob.list[order(pos)]
	score.cats <- score.cats[order(pos)]
	
	# Return results
	list(prob.list=prob.list, score.cats=score.cats)

}


# "lord.wingersky" function
# Arg:
# prob: (matrix or data.frame) probability matrix for all categories for all items
# 		when data.frame is used all elements should be numeric.
# score.cats: (vector) containes score categories for all items 

lord.wingersky <- function(prob, score.cats) {

	if(is.data.frame(prob)) {
		logic <- sapply(1:ncol(prob), function(i) is.numeric(prob[[i]]))
		if(!all(logic)) stop("Data frame of probabilies should be all numeric") 
	}
	
	if(nrow(prob) != length(score.cats)) {
		stop(paste0("There are ", nrow(prob), " items in the probability matrix (or data.frame), ", 
			"whereas there are ", length(score.cats), " items in the category vector."))
	}
	
	if(min(score.cats) < 2) {
		stop("Minimum number of categories for each item is 2")
	}
	
	# Probabilities for each category at 1st item
	p <- as.double(prob[1, 1:score.cats[1]])
	
	# Create a temporary vector to contain probabilities 
	tmp <- c()
	
	# Possible observed score range for a first item
	obs.range <- 0:(score.cats[1]-1)

	# Caculate probabilities to earn an observed score by accumulating over remaining items
	for(j in 2:nrow(prob)) {  # should start from a 2nd item
		
		# Probability to earn zero score. This is a special case
		tmp[1] <- p[1]*prob[j, 1]
		
		# The range of category scores for an added item
		cat.range <- 0:(score.cats[j]-1)
		
		# Possible minimum and maximum observed scores when the item is added
		# but, except zero and perfect score
		min.s <- min(obs.range) + min(cat.range) + 1
		max.s <- max(obs.range) + max(cat.range) - 1
		
		# Calculate probability to obtain possible observed score after the item is added
		# but, still except zero and perfect score
		for(score in min.s:max.s) {
		
			# Difference between the possible observed score and each category score if the added item
			poss.diff <- score - min(cat.range):max(cat.range)
			
			# The difference above should be greater than or equal to the minimum observed score
			# where the item is added, and less than or equal to the maximum observed score where 
			# the item is added
			cols <- poss.diff >= min(obs.range) & poss.diff <= max(obs.range)
			poss.diff <- poss.diff[cols]  # Now, these are remained difference 
			
			# Find columns including the probability of getting the category score of added item
			# in which magnitude is the same as differences obtained from above
			cols.t <- rep(FALSE, ncol(prob))
			cols.t[1:length(cols)] <- cols
			
			# Final probability to earn the observed score
			# "p[rev(poss.diff + 1)]" is the probabilities for possible observed score
			# before the item is added
			prob.score <- sum(p[rev(poss.diff + 1)] * rev(prob[j, cols.t]))
			tmp[1+score] <- prob.score
		}
		
		# Probability to earn perfect score. This is a special case.
		tmp[score + 2] <- p[length(p)]*prob[j, score.cats[j]]
		
		# Update probabilities
		p <- tmp
		
		# Update the range of possible observed scores 
		obs.range <- (min.s-1):(max.s+1)
		
		# Reset the temporary vector
		tmp <- c()
		
	}
	
	p
}
