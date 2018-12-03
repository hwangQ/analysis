# a function to creat a plot for test information functions for all modules

plot.ata_mst <- function(x, xlab="Theta", ylab="TIF", main="Simultaneous Assembley", col=1:4, lwd=2, ...) {
	
	# decomposition of an object
	nstg <- x$metainfo$n.stage
	nmod <- x$metainfo$n.module
	target_rt_sc <- x$target.sc$router
	target_stg2_sc <- x$target.sc$stage2
	info_rt <- x$info.obs$router
	info_stg2 <- x$info.obs$stage2
	theta <- x$theta
	
	# target information
	ylim.max <- max(c(unlist(x$target.sc), unlist(x$info.obs))) + .5
	plot(target_rt_sc ~ theta, type='l', col=col[1], lty=2, lwd=lwd, 
		ylim=c(0, ylim.max),
		ylab=ylab, xlab=xlab, 
		main=main, ...)
	for(i in 1:nmod[1]) {
		lines(target_stg2_sc[[i]] ~ theta, lty=2, lwd=lwd, col=col[1+i])
	}

	# TIF based on selected items
	lines(info_rt ~ theta, col=col[1], lwd=lwd)
	for(i in 1:nmod[1]) {
		lines(info_stg2[[i]] ~ theta, col=col[i+1], lwd=lwd)
	}
	
	# make a legend
	legend("topleft", legend=c("Router", paste0("Stage2 (M", 1:nmod[1], ")")), title="Module", lty=1, col=col)
	legend("topright", legend=c("Router", paste0("Stage2 (M", 1:nmod[1], ")")), title="Target", lty=2, col=col)

}


# a modified version of fuction to create a figure of panel configuration (see mstR package)
plot.mst<-function (x, show.path = TRUE, border.col = "red", arrow.col = "red", 
                    module.names = NULL, save.plot = FALSE, save.options = c("path", 
                                                                             "name", "pdf"), cex.m=2, ...) 
{
  internalMST <- function() {
    nr <- 0
    if (class(x)!="mst" & class(x)!="matrix") stop("'x' must be a list of class 'mst' or a transition matrix",call.=FALSE)
    if (class(x)=="matrix") {
      x<-list(transMatrix=x)
      SHOW.path<-FALSE
    }
    else SHOW.path<-show.path
    tr <- x$transMatrix
    nr.st <- NULL
    repeat {
      nr <- nr + 1
      ind <- which(colSums(tr) > 0)
      nr.st[nr] <- ncol(tr) - length(ind)
      tr <- tr[ind, ind]
      if (sum(tr) == 0) 
        break
    }
    nr.st <- c(nr.st, ncol(tr))
    height <- 2 * length(nr.st) + 3 * (length(nr.st) - 1)
    width <- 2 * (max(nr.st) * 2 - 1)
    xl <- c(-width/2, width/2)
    yl <- c(-height/2, height/2)
    plot(0, 0, xlim = xl, ylim = yl, xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "", bty = "n", col = "white")
    xcenter <- ycenter <- NULL
    allleft <- allright <- allup <- alldown <- NULL
    for (NR in 1:length(nr.st)) {
      up <- yl[2] - 5 * (NR - 1)
      down <- up - 2
      left <- seq(from = -nr.st[NR] * 2 + 1, length = nr.st[NR], 
                  by = 4)
      right <- left + 2
      rect(left, down, right, up)
      xcenter <- c(xcenter, (left + right)/2)
      ycenter <- c(ycenter, rep((up + down)/2, nr.st[NR]))
      allleft <- c(allleft, left)
      allright <- c(allright, right)
      allup <- c(allup, rep(up, nr.st[NR]))
      alldown <- c(alldown, rep(down, nr.st[NR]))
    }
    for (i in 1:nrow(x$transMatrix)) {
      for (j in 1:ncol(x$transMatrix)) {
        if (x$transMatrix[i, j] == 1) 
          arrows(xcenter[i], ycenter[i] - 1.1, xcenter[j], 
                 ycenter[j] + 1.1, length = 0.1, angle = 20)
      }
    }
    for (i in 1:length(xcenter)) {
      if (is.null(module.names)) 
        text(xcenter[i], ycenter[i], paste("Module", 
                                           i), cex=cex.m)
      else text(xcenter[i], ycenter[i], module.names[i], cex=cex.m)
    }
    if (SHOW.path) {
      for (i in 1:length(x$selected.modules)) {
        ind <- x$selected.modules[i]
        rect(allleft[ind], alldown[ind], allright[ind], 
             allup[ind], lwd = 2, border = border.col)
      }
      for (i in 1:(length(x$selected.modules) - 1)) {
        ind <- x$selected.modules[i]
        ind2 <- x$selected.modules[i + 1]
        arrows(xcenter[ind], ycenter[ind] - 1.1, xcenter[ind2], 
               ycenter[ind2] + 1.1, length = 0.1, angle = 20, 
               lwd = 2, col = arrow.col)
      }
    }
  }
  internalMST()
  if (save.plot) {
    plotype <- NULL
    if (save.options[3] == "pdf") 
      plotype <- 1
    if (save.options[3] == "jpeg") 
      plotype <- 2
    if (is.null(plotype)) 
      cat("Invalid plot type (should be either 'pdf' or 'jpeg').", 
          "\n", "The plot was not captured!", "\n")
    else {
      if (save.options[1] == "path") 
        wd <- paste(getwd(), "/", sep = "")
      else wd <- save.options[1]
      nameFile <- paste(wd, save.options[2], switch(plotype, 
                                                    `1` = ".pdf", `2` = ".jpg"), sep = "")
      if (plotype == 1) {
        {
          pdf(file = nameFile)
          internalMST()
        }
        dev.off()
      }
      if (plotype == 2) {
        {
          jpeg(filename = nameFile)
          internalMST()
        }
        dev.off()
      }
      cat("The plot was captured and saved into", "\n", 
          " '", nameFile, "'", "\n", "\n", sep = "")
    }
  }
  else cat("The plot was not captured!", "\n", sep = "")
}

