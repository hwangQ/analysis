# This function computes classfication accuracy from the simulation results
ccr <- function(x, cut_thread) {

  ## Arg
  # x: a list of simulation results
  # cut_thread: a vector of cut scores
    
  # extract required infomation from x object
  true.theta <- x$true.theta
  est.theta <- x$est.theta  
  
  # creat an empty matrix to contain classification results
  ccr_mat <- matrix(NA, length(cut_thread), 4)
  colnames(ccr_mat) <- c("CCR", "FP", "FN", "TER")
  rownames(ccr_mat) <- paste0("CUT (", cut_thread, ")")
  
  # compute the classification accuracy
  for(i in 1:length(cut_thread)) {
    
    # indexing for classification
    class_true <- ifelse(true.theta >= cut_thread[i], 1, 0)
    class_est <- ifelse(est.theta >= cut_thread[i], 1, 0)
    class_diff <- class_true - class_est
    
    # classification accuracy
    CCR <- mean(class_diff == 0)
    FP <- mean(class_diff == -1)
    FN <- mean(class_diff == 1)  
    TER <- FP + FN 
    ccr_mat[i, ] <- c(CCR, FP, FN, TER)
    
  }
  
  # return the results
  ccr_mat
  
}
