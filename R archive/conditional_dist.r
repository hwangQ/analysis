# a function to estimate conditional ability probability given the true ability for each module

ability_dist_sep <- function(x, theta, nstg, nmod, D) {

   ## Arguement
   # x: a list of item paramameter data.frame. Each element in a list contains a dafa.frame of item parameter information 
   #    obtained from shape_df function
   # theta: a vector containing theta (or ability) values
   # nstg: a scalar indicating the number of stage in MST
   # nmod: a scalar or vector indicating the number of modules in each stage
   # D: a scalar of scaling factor in IRT model
   
   if(nstg != length(x)) stop("The number of stages should be the same with the length of list of x") 
   if(nmod[1] != length(x[[2]])) stop("The number of modules in the second stage should be the same with the number of item parameter data frames")
   if(length(nmod) == 2) { 
      if(nmod[2] != length(x[[3]])) stop("The number of modules in the third stage should be the same with the number of item parameter data frames")
   }

   ## prepare the calculation of Lord-wingersky recursive algorithm
   # (a) Router
   preplw_rt <- prep4lw(x[[1]], theta, D)
   
   # (b) Second stage
   for(i in 1:nmod[1]) {
       assign(paste0("preplw_stg2_", i), prep4lw(x[[2]][[i]], theta, D))
   }

   # (c) Third stage
   if(length(nmod) == 2) { 
      for(i in 1:nmod[2]) {
          assign(paste0("preplw_stg3_", i), prep4lw(x[[3]][[i]], theta, D))
      }
   }

   ## conditional probability distribution given each theta value
   # (a) Router
   cond.prob_rt <- lwRecurive(preplw_rt$prob.list, preplw_rt$score.cats) 

   # (b) Second stage
   cond.prob_stg2 <- vector('list', nmod[1])
   for(i in 1:nmod[1]) {
       cond.prob_stg2[[i]] <- lwRecurive(get(paste0("preplw_stg2_", i))$prob.list, get(paste0("preplw_stg2_", i))$score.cats)
   }
   
   # (c) Third stage
   if(length(nmod) == 2) { 
      cond.prob_stg3 <- vector('list', nmod[2])
      for(i in 1:nmod[2]) {
          cond.prob_stg3[[i]] <- lwRecurive(get(paste0("preplw_stg3_", i))$prob.list, get(paste0("preplw_stg3_", i))$score.cats)
      }
   }

   ## create a list containing conditional probability distribution given each true ability
   # (a) When the number of stages are two
   if(length(nmod) == 1) {
      cond.dist_stg2_list <- lapply(1:length(theta), function(i) sapply(cond.prob_stg2, '[', 1:nrow(cond.prob_stg2[[1]]), i))
      cond.dist_list <- lapply(1:length(theta), function(i) rowr::cbind.fill(cond.prob_rt[, i], cond.dist_stg2_list[[i]], fill=0L))
      col.names <- c("router", paste0("stage2.", 1:nmod[1]))
      cond.dist_list <- lapply(1:length(theta), function(i) setNames(cond.dist_list[[i]], col.names))
   }

   # (b) When the number of stages are three
   if(length(nmod) == 2) {
      cond.dist_stg2_list <- lapply(1:length(theta), function(i) sapply(cond.prob_stg2, '[', 1:nrow(cond.prob_stg2[[1]]), i))
      cond.dist_stg3_list <- lapply(1:length(theta), function(i) sapply(cond.prob_stg3, '[', 1:nrow(cond.prob_stg3[[1]]), i))
      cond.dist_list <- lapply(1:length(theta), function(i) rowr::cbind.fill(cond.prob_rt[, i], cond.dist_stg2_list[[i]], cond.dist_stg3_list[[i]], fill=0L))
      col.names <- c("router", paste0("stage2.", 1:nmod[1]), paste0("stage3.", 1:nmod[2]))
      cond.dist_list <- lapply(1:length(theta), function(i) setNames(cond.dist_list[[i]], col.names))
  }
   
   names(cond.dist_list) <- theta
   if(length(nmod) == 1) {
      scores <- list(router=0:(nrow(cond.prob_rt)-1), stage2=0:(nrow(cond.prob_stg2[[1]])-1))
   }
   if(length(nmod) == 2) {
      scores <- list(router=0:(nrow(cond.prob_rt)-1), stage2=0:(nrow(cond.prob_stg2[[1]])-1), stage3=0:(nrow(cond.prob_stg3[[1]])-1))
   }
   
   rr <- list(cond_dist=cond.dist_list, scores=scores, theta=theta, n.stage=nstg, n.module=nmod)
   rr
   
}


# a function to estimate conditional ability probability given the true ability for total test score

ability_dist <- function(x, eos.router, cut.score) {

   ## Arguement
   # x: an object obtained from ability_dist_sep function
   # cut.score: a vector containing cut scores on the theta (abilit) scale
   # eos.router: a vector containing the equated observed scores of a router
   
   # set conditions
   nstg <- x$n.stage
   nmod <- x$n.module
   theta <- x$theta
   score_rt <- x$scores$router # router scores
   score_stg2 <- x$scores$stage2 # stage2 scores
   score_total <- 0:sum(max(score_rt), max(score_stg2)) # total scores
   cond_list <- x$cond_dist
   
   # give a router score a corresponding path
   path <- give_path(score=eos.router, cut.score=cut.score)$path

   # (a) When the number of stages are two: need to make this as a function!
   if(nstg == 2) { 
      theta_condist <- vector('list', length(theta))
      for(j in 1:length(theta)) {
          comb_dist <- array(0, c(length(score_total), nmod[1]))
          rownames(comb_dist) <- score_total
          for(i in 1:length(score_rt)) {
              comb_dist[i:(i+max(score_stg2)), path[i]] <- comb_dist[i:(i+max(score_stg2)), path[i]] + 
                                                          (cond_list[[j]][1:length(score_stg2), (path[i]+1)] * cond_list[[j]][i, 1])
          }
          theta_condist[[j]] <- comb_dist
      }
   }
   
   # (b) When the number of stages are three (developing now)
   if(nstg == 3) { 
      next
   }
   
   names(theta_condist) <- round(theta, 5)
   theta_condist
   
}


