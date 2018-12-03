# a function used for evaluation
paramList <- function(x){

	modelGood <- all(x$MODEL %in% c('1PLM', '2PLM', '3PLM', 'GRM', 'GPCM'))
	catsGood <- all(x$CATEGORY >= 1)
	if(!modelGood) stop("At least, one of model is mis-specified. Available models are 1PLM, 2PLM, 3PLM, GRM, and GPCM")
	if(!catsGood) stop("At least, one of score category is less than 2. Score category should be greater than 1")
	
	x$ID <- as.character(x$ID) # change ID as character
	x$MODEL <- as.character(x$MODEL) # change MODEL as character
	any.dc <- any(x$CATEGORY == 2) # if there are binary items
	any.py <- any(x$CATEGORY > 2) # if there are polytomous items
	
	paramList <- list()
	
	if(any.dc) { # binary item
		
		dicho.dat <- x[x$CATEGORY == 2,]
		Ndicho <- nrow(dicho.dat)
		dicho.item <- list()
		dicho.item$a <- 0
		dicho.item$b <- 0
		dicho.item$g <- 0
		dicho.item$id <- 0
		dicho.item$dModel <- 0 
		dicho.item$dPos <- which( x$CATEGORY == 2 )
		for(i in 1:Ndicho) {
			
			start.par <- which(names(dicho.dat) == "PARAM.1")
			end.par <- ncol(dicho.dat)
			dicho.item$a[i] <- dicho.dat[i, start.par]
			dicho.item$b[i] <- dicho.dat[i, start.par + 1]
			if(dicho.dat$MODEL[i] == "3PLM") {
				dicho.item$g[i] <- dicho.dat[i, start.par + 2]
			} else {
				dicho.item$g[i] <- 0
			}
			dicho.item$id[i] <- dicho.dat$ID[i]
			dicho.item$dModel[i] <- dicho.dat$MODEL[i]
		}
	
	paramList$dicho.item <- dicho.item 
	
	}

	if(any.py) { # polytomous item
		
		poly.dat <- x[x$CATEGORY > 2,]
		Npoly <- nrow(poly.dat)
	
		poly.item <- list()
		poly.item$a <- 0
		poly.item$d <- list()
		poly.item$score.cats <- 0
		poly.item$id <- 0
		poly.item$pModel <- 0
		poly.item$pPos <- which( x$CATEGORY > 2 )
		
		for(i in 1:nrow(poly.dat)) {
			
			start.par <- which(names(poly.dat) == "PARAM.1")
			end.par <- ncol(poly.dat)
			cats <- poly.dat$CATEGORY[i]
			pModel <- poly.dat$MODEL[i]
			poly.item$a[i] <- poly.dat[i, start.par]
			if(pModel == "GRM") poly.item$d[[i]] <- unlist(poly.dat[i, (start.par+1):((start.par+1)+cats-2)])
			if(pModel == "GPCM") poly.item$d[[i]] <- c(0, unlist(poly.dat[i, (start.par+1):((start.par+1)+cats-2)]))
			poly.item$score.cats[i] <- cats
			poly.item$pModel[i] <- pModel
			poly.item$id[i] <- poly.dat$ID[i]
		}
	
	paramList$poly.item <- poly.item 
	}
	
	return(paramList)

}


