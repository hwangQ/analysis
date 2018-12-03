# a function to find roots where two test information functions are crossed
cross_info <- function(x, range.theta=c(-5, 5), D, interval=c(-3, 3), ...) {

   ## Arguement
   # x: a list of data.frame of item paramameters. Each element in a list contains a dafa.frame of item parameter information 
   #    for each module at the same stage
   # range.theta: a vector containing the range of theta values to be used for calculating test infomatmation
   # D: a scalar of scaling factor in IRT model

   # check the number of modules
   nmod <- length(x)

   # estimate cut scores
   cut.score <- NA
   for(i in 1:(nmod-1)) {
       cut.score[i] <- rootSolve::uniroot.all(f=diff_info, interval=interval, df.x=x[[i]], df.y=x[[i+1]], D=D, ...)
   }
   
   # crate a matrix of test information function for each module
   theta <- seq(range.theta[1], range.theta[2], 0.1)
   info.mat <- array(NA, c(length(theta), nmod))
   colnames(info.mat) <- paste0("module.", 1:nmod)
   for(i in 1:nmod) {
       info.mat[, i] <- test.info(x=x[[i]], theta=theta, D=D)$testInfo
   }
   
   rr <- list(testInfo=info.mat, cut.score=cut.score)
   rr
   
}

# a function to calculate the differnce of two test information values at the same ability point
diff_info <- function(theta, df.x, df.y, D) {
  
  rr <- test.info(x=df.x, theta=theta, D=D)$testInfo - test.info(x=df.y, theta=theta, D=D)$testInfo
  rr <- as.numeric(rr)
  rr
  
}

