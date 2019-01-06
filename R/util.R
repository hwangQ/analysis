# This function creates a table for a parition of contents
content_table <- function(x, which.mst, mod.name=NULL) {
  
  # select MSTs to be used to create a table
  obj <-x[which.mst]   
  
  # a function to create a talbe for one MST
  f <- function(dat, mod.name=NULL) {
    
    # read RDPs
    post <- dat$metainfo$post
    RDP <- post[-c(1, length(post))]
    RDP_str <- paste0("RDP (", paste(RDP, collapse = ", "), ")") 
    
    # read prm.df
    prm_df <- dat$prm.df
    
    # check the number of modules
    nmod <- length(prm_df)
    
    # creat a table
    table_df <- data.frame()
    for(i in seq_len(nmod)) {
      N <- nrow(prm_df[[i]])
      cont1 <- table(factor(prm_df[[i]]$CLASS1, levels=1:4)) %>% 
        matrix(nrow=1)
      cont2 <- table(factor(prm_df[[i]]$CLASS2, levels=1:3)) %>% 
        matrix(nrow=1)
      tmp_df <- cbind(N, cont1, cont2) %>% 
        data.frame()
      table_df <- rbind(table_df, tmp_df)
    }
    
    if(is.null(mod.name)) {
      mod.name <- paste0("module.", seq_len(nmod))
      table_df <- data.frame(mod.name, table_df)
    } else {
      table_df <- data.frame(mod.name, table_df)
    }
    
    col.names <- c("Module", "N", paste0("Cont1.", 1:4), paste0("Cont2.", 1:3))
    colnames(table_df) <- col.names
    table_df <- data.frame(RDP=RDP_str, table_df, stringsAsFactors = FALSE)
    
    table_df
  } 
  
  # creat a table for all MSTs
  df <- purrr::map_dfr(obj, .f=f, mod.name=mod.name)
  
  df
  
}


# This function summarizes the results of the objective functions in the specified order
summary_obj <- function(obj_res, order=NULL, showRDP=FALSE, RDP_mat=NULL) {
  
  if(showRDP) {
    RDP_df <- 
      as.data.frame(round(RDP_mat, 5)) %>% 
      tidyr::unite(col=RDP, sep=", ")
    # RDP_df <- data.frame("(", RDP_df, ")") %>% 
    #   tidyr::unite(col=RDP, sep="")
  }
  
  f <- function(i, showRDP) {
    if(!showRDP) {
      obj_res %>% 
        dplyr::summarize(mar.rel=sort(mrel, decreasing=TRUE, na.last=TRUE)[i], 
                         loc.mrel=order(mrel, decreasing=TRUE, na.last=TRUE)[i], 
                         ave.csee=sort(ave_csee, na.last=TRUE)[i], 
                         loc.ave=order(ave_csee, na.last=TRUE)[i],
                         max.csee=sort(max_csee, na.last=TRUE)[i], 
                         loc.max=order(max_csee, na.last=TRUE)[i])
    } else {
      
      obj_res %>% 
        dplyr::summarize(mar.rel=sort(mrel, decreasing=TRUE, na.last=TRUE)[i], 
                         rdp.mrel=RDP_df[order(mrel, decreasing=TRUE, na.last=TRUE)[i], ], 
                         ave.csee=sort(ave_csee, na.last=TRUE)[i], 
                         rdp.ave=RDP_df[order(ave_csee, na.last=TRUE)[i], ],
                         max.csee=sort(max_csee, na.last=TRUE)[i], 
                         rdp.max=RDP_df[order(max_csee, na.last=TRUE)[i], ])   
    }
  }
  
  if(is.null(order)) {
    
    res <- f(1, showRDP, order=NULL, showRDP)
    res$order <- 1
    
  } else {
    
    res <- purrr::map_df(order, f, showRDP=showRDP)
    res$order <- order
    
  }
  
  res
  
}

# This function finds theta values for each targeted subpopulations
subpop <- function(post, n.stage) {
  
  if(n.stage ==2) {
    thetas <- list()
    for(i in 1:(length(post) - 1)) {
      thetas[[i]] <- seq(post[i], post[i + 1], length.out = 3)[2]
    }
  }
  
  if(n.stage == 3) {
    med <- NA
    for(i in 1:(length(post) - 1)) {
      med[i] <- median(c(post[i], post[i + 1]))
    }
    post2 <- sort(c(med, post[-c(1, length(post))]))
    thetas <- as.list(post2)
    
    # post2 <- sort(c(post, med))
    # thetas <- list()
    # for(i in 1:(length(post2) - 2)) {
    #   thetas[[i]] <- seq(post2[i], post2[i + 2], length.out = 3)[2]
    # }
  }
  
  # return results
  names(thetas) <- paste0("g", 1:length(thetas))
  thetas
  
}

# This function creates a matrix of cutoff scores used for routing
cutoff_mat <- function(cut.score, n.stage, n.module) {
  
  cut.mat <- matrix(NA, 1, 2)
  for(s in 2:n.stage) {
    cut.tmp <- c(-Inf, cut.score[[s]], Inf)
    tmp <- matrix(NA, n.module[s], 2)
    for(i in 1:n.module[s]) {
      tmp[i, ] <- cut.tmp[i:(i + 1)]
    }
    cut.mat <- rbind(cut.mat, tmp)
  }
  
  cut.mat
  
}


# a function to create target module information function (MIF)
info_mif <- function(item.pool, target.theta, route.map, test.length, D=1) {

   # extract needed information from a routing map
   panel.info <- panel_info(route.map)
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # create target MIF
   target.info <- vector('list', n.stage)
   names(target.info) <- paste0("stage", 1:n.stage)
   for(s in 1:n.stage) {
       for(m in 1:n.module[s]) {
           target.info[[s]][[m]] <- rep(NA, length(target.theta[[s]][[m]]))
           for(i in 1:length(target.theta[[s]][[m]])) {
               tmp.1 <- test.info(x=item.pool, theta=target.theta[[s]][[m]][i], D=D)$itemInfo
               tmp.2 <- sort(tmp.1, decreasing=TRUE)[1:test.length]
               target.info[[s]][[m]][i] <- size.module[s] * mean(tmp.2)
           }
       }
   }

   # return results
   rr <- list(theta=target.theta, info=target.info)

   rr
   
}

# a function to assign a path to each router score
give_path <- function(score, cut.score) { 

	# number of scores
	nscore <- length(score)

	# number of categories (or levels)
	cats <- length(cut.score) + 1

	# set total lables to be given to each examinee
	level <- 1:cats

	# create path vectors
	path <- rep(level[1], nscore)

	# give a level correspoding each score
	for(i in 1:length(cut.score)) {
		path[which(score >= cut.score[i])] <- level[i+1]
	}

	rr <- list(path=path, score=score, cut.score=cut.score)
	rr

}

# calculate moments
cal_moments <- function(node, weight) {

	mu <- sum(node * weight)
	sigma2 <- sum(node^2 * weight) - mu^2
	rr <- c(mu=mu, sigma2=sigma2)
	rr

}

# This function computes three objective functions
objfn <- function(obj=c("mrel", "ave.se", "max.se"), var.cond, var.pop=1, w, range=c(-2, 2)) {
  
  obj <- match.arg(obj)
  
  if(obj %in% c("ave.se", "max.se")) {
    selected <- names(var.cond) %in% as.character(seq(range[1], range[2], .1))
    var.cond <- var.cond[selected]
    csee <- sqrt(var.cond)
  }
  
  res <- switch(obj,
                mrel = (var.pop - stats::weighted.mean(x=var.cond, w=w)) / var.pop,
                ave.se = mean(csee),
                max.se = max(csee)
  )
  
  res
  
}

# marginal reliability
mrel <- function(var.pop, var.cond, wieghted=TRUE, w) {
	
	if(wieghted) {
		rr <- (var.pop - weighted.mean(x=var.cond, w=w)) / var.pop
	} else {
		rr <- (var.pop - mean(x=var.cond)) / var.pop
	}
	
	rr

}


# calculation of four moments of normal distribution
norm_moments <- function(mu, sigma) {

  m1 <- mu
  m2 <- mu^2 + sigma^2
  m3 <- mu * (mu^2 + 3 * sigma^2)
  m4 <- mu^4 + 6 * mu^2 * sigma^2 + 3 * sigma^4
  
  rr <- c(m1=m1, m2=m2, m3=m3, m4=m4)
  rr
  
}