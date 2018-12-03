# a function to calculate joint probability distribution for all paths given ability points
joint_dist <- function(cond.dist, path.label, n.route, max.score, score.pre, score.post) {
   
   ## Arguement
   # cond.dist: a list including conditional probability of raw scores for modules across all ability points
   # path.label: a vector including labels of routings corresponding to raw scores
   # n.route: a scalar including the number of routings
   # max.score: a scalar including the max raw scores for the cumulated test form
   # score.pre: a vector including all score points for pre-stage module
   # score.post: a vector including all score points for prost-stage module
   
   cond_joint <- vector('list', length(theta))
   for(j in 1:length(theta)) {
       jointdist <- array(0, c(max.score + 1, n.route))
       rownames(jointdist) <- 0:max.score
       for(i in 1:length(score.pre)) {
           jointdist[i:(i+max(score.post)), path.label[i]] <- jointdist[i:(i+max(score.post)), path.label[i]] + 
                                                             (cond.dist[[j]][1:length(score.post), (path.label[i]+1)] * cond.dist[[j]][i, 1])
       }
       cond_joint[[j]] <- jointdist
   }
   
   cond_joint
   
}


# a function to estimate conditional ability probability given the true ability for each module
ability_dist_sep <- function(module.form, theta, route.map, D) {

   ## Arguement
   # module.form: a list including a data.frame of items of modules across all stages
   # theta: a vector containing theta (or ability) values
   # route.map: a matrix of the routing information
   # D: a scalar of scaling factor in IRT model

   # obtain a panel information using a routing map
   panel.info <- panel_info(route.map)
   pathway <- panel.info$pathway
   config.info <- panel.info$config.info
   nstg <- panel.info$n.stage
   
   x <- vector('list', nstg)
   names(x) <- paste0("stage", 1:nstg)
   for(s in 1:nstg) {
       idx <- config.info[[s]]
       x[[s]] <- lapply(idx, function(i) module.form[[i]])
   }
   
   # warning messages
   if(nstg != length(x)) stop("The number of stages should be the same with the length of list of x") 

   # prepare the calculation of Lord-wingersky recursive algorithm
   preplw_list <- vector('list', nstg)
   names(preplw_list) <- paste0("stage", 1:nstg)
   for(s in 1:nstg) {
       preplw_list[[s]] <- lapply(x[[s]], prep4lw, theta=theta, D=D)
   }

   # conditional probability raw scores given each theta value across all modules
   cond.prob_list <- vector('list', nstg)
   names(cond.prob_list) <- paste0("stage", 1:nstg)
   for(s in 1:nstg) {
      num <- length(preplw_list[[s]])
      cond.prob_list[[s]] <- lapply(1:num, function(i) lwRecurive(preplw_list[[s]][[i]]$prob.list, preplw_list[[s]][[i]]$score.cats))
      names(cond.prob_list[[s]]) <- paste0("m", 1:num)
   }

   # a list containing a matrix of conditional probability of raw scores at each module given each ability across all stages
   cond.dist_stg_list <- vector('list', nstg)
   for(s in 1:nstg) {
       cond.dist_stg_list[[s]] <- lapply(1:length(theta), function(i) sapply(cond.prob_list[[s]], '[', 1:nrow(cond.prob_list[[s]][[1]]), i))
   }
   
   # a list containing a matrix of conditional probability of raw scores for all modules across all ability points
   cond.dist_list <- vector('list', length(theta))
   names(cond.dist_list) <- round(theta, 5)
   for(i in 1:length(theta)) {
      arg <- c(lapply(cond.dist_stg_list, '[[', i), fill=0L)
      cond.dist_list[[i]] <- do.call(rowr::'cbind.fill', arg=arg)
   }
   
   # a list of raw score points for each stage
   scores <- lapply(1:nstg, function(i) 0:nrow(x[[i]][[1]]))
   names(scores) <- paste0("stage", 1:nstg)

   rr <- list(cond_dist=cond.dist_list, scores=scores, theta=theta, route.map=route.map)
   rr
   
}

# a function to estimate conditional ability probability given the true ability for total test score
ability_dist <- function(x, eos.path, n.path, cut.score) {

   ## Arguement
   # x: an object obtained from ability_dist_sep function
   # eos.path: a list including matrixes of equated observed scores for all (sub) paths
   # n.path: a vector including the number of routing (path) ways to next stage
   # cut.score: a list including routing points along with ability scale at each stage

   # obtain a panel information using a routing map
   panel.info <- panel_info(x$route.map)
   pathway <- panel.info$pathway
   nstg <- panel.info$n.stage
   nmod <- panel.info$n.module
   
   # set conditions
   theta <- x$theta
   score <- x$scores
   score_cum <- cumsum(sapply(score, max))[-1] # cumulative sum scores
   cond_list <- x$cond_dist

   # an empty list to include routing points at all routes across all stages
   cut4route <- vector('list', nstg-1)
   names(cut4route) <- paste0("stage", 2:nstg)

   # an empty list to include path labels across all stages
   pathLabel <- vector('list', nstg-1)
   names(pathLabel) <- paste0("stage", 2:nstg)
   
   # an empty list to include probabilities of joint distribution
   jointDist <- vector('list', (nstg-1))
   names(jointDist) <- paste0("stage", 2:nstg)

   # calculate probabilities of joint distribution
   for(s in 2:nstg) {
   
       if(s == 2) {
          # for 2nd stage (default)
          path.label <- give_path(score=eos.path[[s-1]], cut.score=cut.score[[s]])$path
          max.score <- score_cum[s-1]
          n.route <- nmod[s]
          score.pre <- score[[s-1]]
          score.post <- score[[s]]
          cond.dist <- lapply(cond_list, '[', 1:nrow(cond_list[[1]]), 1:(1+n.route))
          jointDist[[s-1]] <- joint_dist(cond.dist, path.label, n.route, max.score, score.pre, score.post)
          names(jointDist[[s-1]]) <- round(theta, 5)
          cut4route[[1]] <- cut.score[[s]]
          pathLabel[[1]] <- path.label
          
       }

       if(s == 3) {
   
          n_path <- n.path[s-1]
          cut4route[[s-1]] <- vector('list', n_path)
          nmod.pre <- unique(pathway[, s-1])
          joint_temp <- vector('list', n_path)
          for(i in 1:n_path) {
  
              # find cut point(s) for routing a next module
              nmod.post <- pathway[pathway[, s-1] == nmod.pre[i], s]
              temp.num <- match(nmod.post, unique(pathway[, s]))
              temp.vec <- cut.score[[s]][min(temp.num):(min(temp.num) + length(temp.num))]
              cut4route[[s-1]][[i]] <- temp.vec[1:(length(temp.num)-1)]

              path.label <- give_path(score=eos.path[[s-1]][, i], cut.score=cut4route[[s-1]][[i]])$path
              pathLabel[[s-1]][[i]] <- path.label
              max.score <- score_cum[s-1]
              n.route <- length(nmod.post)
              score.pre <- 0:score_cum[s-2]
              score.post <- score[[s]]
              temp.dist1 <- lapply(jointDist[[s-2]], '[', 1:nrow(jointDist[[s-2]][[1]]), i)
              temp.dist2 <- lapply(cond_list, '[', 1:nrow(cond_list[[1]]), nmod.post)
              cond.dist <- lapply(1:length(theta), function(j) rowr::cbind.fill(temp.dist1[[j]], temp.dist2[[j]], fill=0L))
              joint_temp[[i]] <- joint_dist(cond.dist, path.label, n.route, max.score, score.pre, score.post)

          }
     
          jointDist[[s-1]] <- lapply(1:length(theta), function(i) do.call('cbind', lapply(joint_temp, '[[', i)))
          names(jointDist[[s-1]]) <- round(theta, 5)

       }
   
   }
   
   rr <- list(joint.dist=jointDist, cut4route=cut4route, path.label=pathLabel, 
              score=score, maxscore.cum=cumsum(sapply(score, max)), theta=theta, n.stage=nstg, n.module=nmod)
   rr
   
}
