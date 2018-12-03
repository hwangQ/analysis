# a function to find roots where two test information functions are crossed
cross_info <- function(x, RDP=NULL, range.theta=c(-5, 5), D, interval=c(-3, 3), ...) {

   ## Arguement
   # x: a list of item paramameter data.frames for modules
   # RDP: a scala or vector including the rouding decision points
   # range.theta: a vector containing the range of theta values to be used for calculating test infomatmation
   # D: a scalar of scaling factor in IRT model

   # check the number of modules
   nmod <- length(x)

   # estimate cut scores
   cut.score <- NA
   for(i in 1:(nmod-1)) {
       temp.cut <- rootSolve::uniroot.all(f=diff_info, interval=interval, df.x=x[[i]], df.y=x[[i+1]], D=D, ...)
       if(!is.null(RDP)) {
          num <- which.min(abs(temp.cut - RDP[i]))  
          cut.score[i] <- temp.cut[num]
       } else {
          cut.score[i] <- rootSolve::uniroot.all(f=diff_info, interval=interval, df.x=x[[i]], df.y=x[[i+1]], D=D, ...)
       }
   }
   
   # crate a matrix of test information function for each module
   theta <- seq(range.theta[1], range.theta[2], 0.1)
   info.mat <- array(NA, c(length(theta), nmod))
   colnames(info.mat) <- paste0("m", 1:nmod)
   for(i in 1:nmod) {
       info.mat[, i] <- test.info(x=x[[i]], theta=theta, D=D)$testInfo
   }
   
   rr <- list(testInfo=info.mat, cut.score=cut.score, theta=theta)
   rr

}

# a function to calculate the differnce of two test information values at the same ability point
diff_info <- function(theta, df.x, df.y, D) {
  
  rr <- test.info(x=df.x, theta=theta, D=D)$testInfo - test.info(x=df.y, theta=theta, D=D)$testInfo
  rr <- as.numeric(rr)
  rr
  
}

# a function to obtain cut scores across all stages
cutoff <- function(module.form, route.map, RDP, D, range.theta, interval) {

   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway
   config.info <- panel.info$config.info
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module

   # create an empty list to include cutscore information
   cut.score <- vector('list', n.stage)
   names(cut.score) <- paste0("stage", 1:n.stage)

   for(s in 1:n.stage) {
       idx <- config.info[[s]]
       forms <- lapply(idx, function(i) module.form[[i]])
       if(length(forms) == 1) {
          cut.score[[s]] <- NA
       } else {
          cut.score[[s]] <- cross_info(forms, RDP=RDP[[s]], range.theta=range.theta, D=D, interval=interval, n=1000)$cut.score
       }
   }

   cut.score

}
