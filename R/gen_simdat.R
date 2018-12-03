# A function for generating response data for a test 
simdat <- function(theta, a.dc, b.dc, g.dc, a.py, b.py, cats, pModel, D) { 

   # check whether argument is correctly specified
   if(missing(cats)) stop("Category of each item is missing") 
   
   # Set conditions
   nstd <- length(theta)
   nitem <- length(cats)

   # Create an empty matrix
   res <- matrix(NA, nrow=nstd, ncol=nitem)

   # Set initial numbers
   idx.dc <- which(cats == 2)
   idx.py <- which(cats > 2)

   # Simulate data
   # (1) for dichotomous items
   if(any(cats == 2)) {
      res.dc <- simdat_dc(theta, a=a.dc, b=b.dc, g=g.dc, D=D)
      res[, idx.dc] <- res.dc 
   }

   # (2) for polytomous items
   if(any(cats > 2)) {
      for(i in 1:length(idx.py)) {
         res[, idx.py[i]] <- simdat_py(theta, a=a.py[i], d=b.py[[i]], D=D, pModel=pModel[i])
      }
   }

   if(nstd == 1 | nitem == 1) res <- as.numeric(res)
   
   res

}

# A function for generating binary data for one item 
simdat_dc <- function(theta, a, b, g, D) {

   # Number of examinees
   nstd <- length(theta)
   
   # Number of items
   nitem <- length(a)

   # Calculate true probability for each category
   fit <- pl.fn(theta, a, b, g, D)

   # Sample random variables from uniform dist
   tmp <- runif(nstd * nitem, 0, 1)
   rv_unif <- matrix(tmp, nrow=nstd, ncol=nitem)
   
   # Simulated Response data for one item
   res <- ifelse(fit >= rv_unif, 1, 0)
   if(nstd == 1 | nitem == 1) res <- as.numeric(res)
   
   res

}

# A function for generating categorical data for one item 
simdat_py <- function(theta, a, d, D, pModel) {

   # Number of examinees
   nstd <- length(theta)

   # add zero values for the first category for GPCM
   if(pModel == "GPCM") d <- c(0, d)

   # Calculate true probability for each category
   fit <- poly.fn(theta, a, d, D, pModel)

   # Sample random variables from uniform dist
   rv_unif <- runif(nstd, 0, 1)

   # Simulated Response data for one item
   ncat <- length(d)
   cumprob <- t(apply(fit, 1, cumsum))
   matTF <- rv_unif >= cumprob 
   res <- apply(matTF, 1, sum)
  
   res
 
}