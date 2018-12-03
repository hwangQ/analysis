# This functon computes the simple moving average of test information function

mov_ave <- function(x, theta, cut.score, D, simple.move = TRUE) {
  
  # number of primary routes
  nform <- length(x)
  
  if(simple.move) {
    # find five theta points to compute the simple moving average
    theta <- theta + c(1:5) - 3
  } else {
    theta <- theta
  }

  # add -Inf and Inf to cut scores
  cuts <- c(-Inf, cut.score, Inf)
  
  # compute test information for all routes
  info.mat <- array(NA, c(length(theta), nform))
  colnames(info.mat) <- paste0("f", 1:nform)
  for(i in 1:nform) {
    info.mat[, i] <- test.info(x=x[[i]], theta=theta, D=D)$testInfo
  }
  
  # give labels to know which route is used for test information at each theta point   
  lv <- cut(theta, cuts, labels=LETTERS[1:nform])
  
  # select test information values from the corresponding routes at each theta point
  info_t <- NULL
  for(i in 1:nform) {
    idx <- lv == levels(lv)[i]
    info_t <- c(info_t, info.mat[idx, i])
  }
  
  # return a simple moving average of test information
  mean(info_t)
  
}
