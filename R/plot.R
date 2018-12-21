
# This function draws plots for CSEE or Bias under Study 1
plot_study1 <- function(x, x.var="theta", y.var=c("CSEE", "Bias"), point.size=1.5, line.size=0.8, 
                        lab.size=15, axis.size=15, legend.size=15, strip.size=10, ylim=NULL) {
  
  if(y.var == "CSEE") ylab <- "Standard Error"
  if(y.var == "Bias") ylab <- "Bias"
  
  if(is.null(ylim)) {
    p <- x %>% 
      ggplot(mapping=aes_string(x="theta", y=y.var)) +
      geom_point(mapping=aes_string(shape="Method"), size=point.size) +
      geom_line(mapping=aes_string(color="Method"), size=line.size) +
      # geom_line(mapping=aes_string(linetype="Method"), size=0.8) + 
      labs(x = expression(theta), y = ylab) +
      theme(axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size)) +
      theme(legend.title = element_text(size=legend.size),
            legend.text = element_text(size=legend.size)) +
      theme_bw() +
      facet_grid(Length ~ Panel) +
      theme(strip.text.x = element_text(size = strip.size, face = 'bold'),
            strip.text.y = element_text(size = strip.size, face = 'bold'))
  } else {
    p <- x %>% 
      ggplot(mapping=aes_string(x="theta", y=y.var)) +
      geom_point(mapping=aes_string(shape="Method"), size=point.size) +
      geom_line(mapping=aes_string(color="Method"), size=line.size) +
      # geom_line(mapping=aes_string(linetype="Method"), size=0.8) + 
      labs(x = expression(theta), y = ylab) +
      ylim(ylim[1], ylim[2]) +
      theme(axis.title = element_text(size=lab.size),
            axis.text = element_text(size=axis.size)) +
      theme(legend.title = element_text(size=legend.size),
            legend.text = element_text(size=legend.size)) +
      theme_bw() +
      facet_grid(Length ~ Panel) +
      theme(strip.text.x = element_text(size = strip.size, face = 'bold'),
            strip.text.y = element_text(size = strip.size, face = 'bold'))
  }
  
  p
  
}

# This function draws CSEE plot for multiple MSTs 
plot_csee <- function(cond_moments, which.mst, RDP_mat, xlab.text, ylab.text, main.text=NULL, ylim=c(0, 1), lab.size=15, 
                      main.size=15, axis.size=15, line.size=1.5, legend.size=15, legend.position="right") {
  
  # select MSTs and RDP points to draw CSEE plot
  cond_moments <- cond_moments[which.mst]
  RDP <- cbind(RDP_mat[which.mst, ])
  
  # creat column names for a data.frame
  col_names <- purrr::map_chr(1:length(which.mst), .f=function(x) paste0("MST ", x, ": RDP (", paste(round(RDP[x, ], 2), collapse = ", "), ")"))
  
  # data manipulation 
  df_csee <- 
    purrr::map_dfc(cond_moments, .f=function(x) sqrt(x[2, ])) %>% 
    stats::setNames(nm=col_names)
  df_csee$theta <- as.numeric(colnames(cond_moments[[1]])) 
  df_csee2 <- reshape2::melt(data=df_csee, variable.name="MST", id.vars="theta", value.name="CSEE")
  
  # set plot conditions
  if(missing(xlab.text)) xlab.text <- expression(theta)
  if(missing(ylab.text)) ylab.text <- 'CSEE'
  
  # draw CSEE plots
  p <- df_csee2 %>% 
    ggplot(mapping=aes_string(x="theta", y="CSEE")) +
    geom_line(mapping=aes_string(color="MST"), size=line.size) +
    labs(title = main.text, x = xlab.text, y = ylab.text) +
    ylim(ylim[1], ylim[2]) +
    theme(plot.title = element_text(size=main.size),
          axis.title = element_text(size=lab.size),
          axis.text = element_text(size=axis.size)) +
    theme(legend.title = element_text(size=legend.size),
          legend.text = element_text(size=legend.size),
          legend.position = legend.position) +
    theme_bw()
  
  p
  
}


# This function creates test information functions for all routes across all assembled MSTs
plot.list <- function(x, which.mst, range.theta=c(-5, 5), D=1, xlab.text, ylab.text, main.text=NULL, lab.size=15, main.size=15, axis.size=15,
                      line.color, line.size=1.5, legend.title, legend.text, legend.size=15, legend.position="right", 
                      layout.col=3, strip.size=12) {
  
  # select MSTs to be used for a plot
  x <- x[which.mst]
  
  info_df <- NULL
  for(i in 1:length(x)) {
    
    # read the assembled test forms
    ata_forms <- x[[i]]$prm.df
    
    # read a routing map
    route.map <- x[[i]]$metainfo$route.map
    panel.info <- panel_info(route.map)
    config.info <- panel.info$config
    n.module <- panel.info$n.module
    pathway <- panel.info$pathway
    n.stage <- length(config.info)
    
    # read RDPs
    post <- x[[i]]$metainfo$post
    RDP <- post[-c(1, length(post))]
    RDP_str <- paste0("MST ", i, ": RDP (", paste(RDP, collapse = ", "), ")") 
    
    # estimate the observed equated scores across all (sub) pathways
    eos_list <- est_eos(ata_forms, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)
    
    # extract item parameter data.frame
    df_path <- eos_list$df_path[[n.stage]]
    
    # compute test information for all routes 
    theta <- seq(-4, 4, 0.1)
    testInfo <- purrr::map_dfc(1:length(df_path), .f=function(k) test.info(df_path[[k]], theta, D)$testInfo) %>% 
      stats::setNames(nm=paste0("Route.", 1:length(df_path)))
    testInfo$theta <- theta
    testInfo <- reshape2::melt(data=testInfo, variable.name="Routes", id.vars="theta", value.name="info")
    testInfo$RDP <- RDP_str
    
    info_df <- rbind(info_df, testInfo)
    
  }
  
  # plot TIF
  # creat route names
  rout_nm <- purrr::map_chr(1:nrow(pathway), .f=function(x) paste(pathway[x, ], collapse = "-"))
  
  # Set plot conditions
  if(missing(xlab.text)) xlab.text <- expression(theta)
  if(missing(ylab.text)) ylab.text <- 'Test Information'
  if(missing(line.color)) {
    line.color <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # line.color <- brewer.pal(n = 9, name = "Greys")[3:9]
  } else {
    line.color <- line.color
  }
  if(missing(legend.title)) legend.title <- "Routes"
  if(missing(legend.text)) legend.text <- rout_nm
  
  p <- info_df %>% 
    ggplot(mapping=aes_string(x="theta", y="info")) +
    geom_line(mapping=aes_string(color="Routes"), size=line.size) +
    labs(title = main.text, x = xlab.text, y = ylab.text) +
    theme(plot.title = element_text(size=main.size),
          axis.title = element_text(size=lab.size),
          axis.text = element_text(size=axis.size)) +
    theme(legend.title = element_text(size=legend.size),
          legend.text = element_text(size=legend.size),
          legend.position = legend.position) +
    theme_bw() +
    facet_wrap(~RDP, ncol=layout.col) +
    theme(strip.text.x = element_text(size = strip.size, face = 'bold')) +
    scale_colour_manual(values=line.color, name = legend.title, labels = legend.text)
  
  p
  
}


# This function creates test information functions for all routes for an MST
plot.atamst <- function(x, range.theta=c(-5, 5), D=1, xlab.text, ylab.text, main.text=NULL, lab.size=15, main.size=15, axis.size=15,
                        line.color, line.size=1.5, legend.title, legend.text, legend.size=15, legend.position="right") {
  
  # read the assembled test forms
  ata_forms <- x$prm.df
  
  # read a routing map
  route.map <- x$metainfo$route.map
  panel.info <- panel_info(route.map)
  config.info <- panel.info$config
  n.module <- panel.info$n.module
  pathway <- panel.info$pathway
  n.stage <- length(config.info)
  
  # estimate the observed equated scores across all (sub) pathways
  eos_list <- est_eos(ata_forms, pathway=pathway, range.theta=range.theta, D=D, constraint=TRUE)
  
  # extract item parameter data.frame
  df_path <- eos_list$df_path[[n.stage]]
  
  # compute test information for all routes 
  theta <- seq(-4, 4, 0.1)
  testInfo <- purrr::map_dfc(1:length(df_path), .f=function(x) test.info(df_path[[x]], theta, D)$testInfo) %>% 
    stats::setNames(nm=paste0("Route.", 1:length(df_path)))
  testInfo$theta <- theta
  testInfo <- reshape2::melt(data=testInfo, variable.name="Routes", id.vars="theta", value.name="info")
  
  # plot TIF
  # creat route names
  rout_nm <- purrr::map_chr(1:nrow(pathway), .f=function(x) paste(pathway[x, ], collapse = "-"))
  
  # Set plot conditions
  if(missing(xlab.text)) xlab.text <- expression(theta)
  if(missing(ylab.text)) ylab.text <- 'Test Information'
  if(missing(line.color)) {
    line.color <-  c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    # line.color <- brewer.pal(n = 9, name = "Greys")[3:9]
  } else {
    line.color <- line.color
  }
  if(missing(legend.title)) legend.title <- "Routes"
  if(missing(legend.text)) legend.text <- rout_nm
  
  p <- testInfo %>% 
    ggplot(mapping=aes_string(x="theta", y="info")) +
    geom_line(mapping=aes_string(color="Routes"), size=line.size) +
    labs(title = main.text, x = xlab.text, y = ylab.text) +
    theme(plot.title = element_text(size=main.size),
          axis.title = element_text(size=lab.size),
          axis.text = element_text(size=axis.size)) +
    theme(legend.title = element_text(size=legend.size),
          legend.text = element_text(size=legend.size),
          legend.position = legend.position) +
    theme_bw() +
    scale_colour_manual(values=line.color, name = legend.title, labels = legend.text)
  
  p
  
}


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

