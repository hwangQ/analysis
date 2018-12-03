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

