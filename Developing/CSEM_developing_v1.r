######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 07/25/2018
#Title: Calculate Conditional Standard Error of Measurement and Find optimized Design 
#       based on Three objective functions for 1-3 MST design structure 
######################################################################################################################################
# Read source code
src.dir <- ":/Dissertation/Analysis/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/Dissertation/Analysis/Temp"
setwd(mydir)

# set condtions
theta <- seq(-4, 4, 0.1)
range.theta <- c(-4, 4)
interval <- c(-2.5, 2.5)
# nmod <- 2
nmod <- 3
# nmod <- 4
nstg <- 2
m.size <- c(20, 20)
D <- 1.702

# read the assembled test forms
ata_forms <- df_sim_list


# create a list of item parametr data.frames across all modules at all stages
df_forms <- vector('list', n.stage)
names(df_forms) <- paste0("stage", 1:n.stage)
module.label <- vector('list', n.stage)
for(i in 1:n.stage) module.label[[i]] <- unique(pathway[, i])
for(i in 1:n.stage) {
	df_forms[[i]] <- vector('list', length(module.label[[i]]))
	names(df_forms[[i]]) <- paste0("m", 1:length(module.label[[i]]))
	for(j in 1:length(df_forms[[i]])) {
		df_forms[[i]][[j]] <- ata_forms[[module.label[[i]][[j]]]]
	}
}

# estimate the observed equated scores across all (sub) pathways
x <- ata_forms
eos_list <- est_eos(x, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)

est_eos <- function(x, pathway, ...) {
   
   ## Arguement
   # x: a list of item paramameter data.frame for all modules across all stages
   # pathway: a maxtrix including all pathways
   
   # check the number of stages
   n.stage <- ncol(pathway)
   
   # create a list of item parameter data.frame of for all possible (sub) paths across all stages
   df_path <- vector('list', n.stage)
   names(df_path) <- paste0("stage", 1:n.stage)
   n_path <- NA # the number of paths at each stage
   for(i in 1:n.stage) {
       path.temp <- unique(as.matrix(pathway[, 1:i]))
       n.temp <- nrow(path.temp)
       df_path[[i]] <- lapply(1:n.temp, function(j) do.call('rbind', x[path.temp[j, ]]))
       n_path[i] <- n.temp
   }
   
   # find observed equated scores using TCC method corresponding to each path at each stage
   eos_path <- vector('list', n.stage)
   names(eos_path) <- paste0("stage", 1:n.stage)
   for(i in 1:n.stage) {
      eos_path[[i]] <- sapply(1:length(df_path[[i]]), function(j) convert_eos(x=df_path[[i]][[j]], ...)$theta)
      rownames(eos_path[[i]]) <- 0:(nrow(eos_path[[i]]) - 1)
      colnames(eos_path[[i]]) <- paste0("path.", 1:ncol(eos_path[[i]]))
   }

   rr <- list(eos_path=eos_path, df_path=df_path, n_path=n_path)
   rr

}

# find points where two adjacent TIFs intersect across all stages 
cut.score <- vector('list', n.stage-1)
names(cut.score) <- paste0("stage", 2:n.stage)
for(i in 2:n.stage) {
	x <- df_forms[[i]]
	cut.score[[i-1]] <- cross_info(x, RDP=RDP, range.theta=range.theta, D=D, interval=interval, n=1000)$cut.score
}

# a function to find roots where two test information functions intersects
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
   
   rr <- list(testInfo=info.mat, cut.score=cut.score)
   rr

}

# conditional raw score distribution of each module given an ability value
x <- df_forms
cond_dist <- ability_dist_sep(x, theta=theta, n.stage=n.stage, D=D)

ability_dist_sep <- function(x, theta, n.stage, D) {

   ## Arguement
   # x: a list of item paramameter data.frames for all modules across all stages
   # theta: a vector containing theta (or ability) values
   # n.stage: a scalar indicating the number of stage in MST
   # D: a scalar of scaling factor in IRT model
   
   # the number of stages
   nstg <- n.stage
   
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

   rr <- list(cond_dist=cond.dist_list, scores=scores, theta=theta, n.stage=nstg)
   rr
   
}


# conditional ability distribution of total-test ability estimates given a true ability
x <- cond_dist
eos.path <- eos_list$eos_path
n.path <- eos_list$n_path

ability_dist(x, eos.path, n.path, cut.score, pathway)

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

ability_dist <- function(x, eos.path, n.path, cut.score, pathway) {

   ## Arguement
   # x: an object obtained from ability_dist_sep function
   # eos.path: a list including matrixes of equated observed scores for all (sub) paths
   # n.path: a vector including the number of routing (path) ways to next stage
   # cut.score: a list including routing points along with ability scale at each stage
   # pathway: a maxtrix including all pathways

   # set conditions
   nstg <- ncol(pathway)
   nmod <- sapply(2:nstg, function(i) length(unique(pathway[, i])))
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
          path.label <- give_path(score=eos.path[[s-1]], cut.score=cut.score[[s-1]])$path
          max.score <- score_cum[s-1]
          n.route <- nmod[s-1]
          score.pre <- score[[s-1]]
          score.post <- score[[s]]
          cond.dist <- lapply(cond_list, '[', 1:nrow(cond_list[[1]]), 1:(1+n.route))
          jointDist[[s-1]] <- joint_dist(cond.dist, path.label, n.route, max.score, score.pre, score.post)
          names(jointDist[[s-1]]) <- round(theta, 5)
          cut4route[[1]] <- cut.score[[s-1]]
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
              temp.vec <- cut.score[[s-1]][min(temp.num):(min(temp.num) + length(temp.num))]
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


		# calculate a mean and sd of conditional ability distribution at each ability value
		cond_moments <- sapply(theta_condist, cal_moments, node=eos_path_mat)
		cond_sigma <- sqrt(cond_moments[2, ])
		csem[, j, k] <- cond_sigma 
	
		# marginal reliability
		w <- dnorm(theta, 0, 1)
		mgrel[k, j] <- mrel(var.pop=1, var.cond=cond_moments[2, ], wieghted=TRUE, w=w)
	
		# unweighted marginal reliability  
		selected <- colnames(cond_moments) %in% as.character(seq(-2, 2, .1))
		unw_mgrel_1[k, j] <- mrel(var.pop=1, var.cond=cond_moments[2, ], wieghted=FALSE)
		unw_mgrel_2[k, j] <- mrel(var.pop=1, var.cond=cond_moments[2, selected], wieghted=FALSE)

		# maximum value of CSEM between -2 and 2 
		csem_max[k, j] <- max(cond_sigma[selected])