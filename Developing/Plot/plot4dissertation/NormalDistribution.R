library(lattice)
library(latticeExtra)
library(gridExtra)

## ------------------------------------------------------------------------------
## 1. Plot for an illustration of Luo and Kim (2018, p. 256)
# set theta values and densities
theta <- seq(-4, 4, by = 0.01)            # Data to set up out normal
dens <- dnorm(theta, 0, 1)

# set a panel function 
panel.fun <- function(x, y, ...) {

  panel.grid(h=-1, v=-1)  
  
  xx1 <- c(-4, x[x>=-4 & x<=-.44], -.44)         #Color area
  yy1 <- c(0, dens[x>=-4 & x<=-.44], 0) 
  panel.xyarea(xx1, yy1, col="ivory", border=NA, lty=1)
  ltext(x=-1.5, y=0.5, labels=paste0("1M-2E-3E", "\n", "1M-2M-3E"), pos=1, offset=2, cex=1.5)   
  
  xx2 <- c(-.44, x[x>=-.44 & x<=.44], .44)         #Color area
  yy2 <- c(0, dens[x>=-.44 & x<=.44], 0) 
  panel.xyarea(xx2, yy2, col="ivory3", border=NA, lty=1)
  ltext(x=0, y=0.5, labels=paste0("1M-2E-3M", "\n", "1M-2M-3M", "\n", "1M-2E-3M"), pos=1, offset=1, cex=1.5)   
  
  xx3 <- c(.44, x[x>=.44 & x<=4], 4)         #Color area
  yy3 <- c(0, dens[x>=.44 & x<=4], 0) 
  panel.xyarea(xx3, yy3, col="ivory4", border=NA, lty=1)
  ltext(x=1.5, y=0.5, labels=paste0("1M-2M-3H", "\n", "1M-2H-3H"), pos=1, offset=2, cex=1.5)   
  
  panel.xyplot(x, y, ...)
  panel.abline(v=c(-.44, .44), col="blue", lty=2, lwd=2, reference = FALSE)
  
}

# Lattice xyplot
p <- xyplot(dens ~ theta,                   
       type = "l",
       ylim=c(0, 0.5),
       xlab=list(label=expression(theta), cex=1.5),
       ylab=list(label="Density", cex=1.5),
       scales=list(cex=c(1.5, 1.5)),
       panel = panel.fun
)

# print out a plot
print(p)

## ------------------------------------------------------------------------------
## 2. Plot for an illustration of the modified mapping strategy
# set theta values and densities
theta <- seq(-4, 4, by = 0.01)            # Data to set up out normal
dens <- dnorm(theta, 0, 1)

# set a panel function 
pfun.primary <- function(x, y, ...) {
  
  panel.grid(h=-1, v=-1)  
  
  xx1 <- c(-4, x[x>=-4 & x<=-.44], -.44)         #Color area
  yy1 <- c(0, dens[x>=-4 & x<=-.44], 0) 
  panel.xyarea(xx1, yy1, col="ivory", border=NA, lty=1)
  ltext(x=-1.5, y=0.5, labels="1M-2E-3E", pos=1, offset=1, cex=1.5)   
  
  xx2 <- c(-.44, x[x>=-.44 & x<=.44], .44)         #Color area
  yy2 <- c(0, dens[x>=-.44 & x<=.44], 0) 
  panel.xyarea(xx2, yy2, col="ivory3", border=NA, lty=1)
  ltext(x=0, y=0.5, labels="1M-2M-3M", pos=1, offset=1, cex=1.5)   
  
  xx3 <- c(.44, x[x>=.44 & x<=4], 4)         #Color area
  yy3 <- c(0, dens[x>=.44 & x<=4], 0) 
  panel.xyarea(xx3, yy3, col="ivory4", border=NA, lty=1)
  ltext(x=1.5, y=0.5, labels="1M-2H-3H", pos=1, offset=1, cex=1.5)   
  
  panel.xyplot(x, y, ...)
  panel.abline(v=c(-.44, .44), col="blue", lty=2, lwd=2, reference = FALSE)
  
}

pfun.secondary <- function(x, y, ...) {
  
  panel.grid(h=-1, v=-1)  
  
  xx1 <- c(-1.22, x[x>=-1.22 & x<=0], 0)         #Color area
  yy1 <- c(0, dens[x>=-1.22 & x<=0], 0) 
  panel.xyarea(xx1, yy1, col="palegreen", border=NA, lty=1)
  ltext(x=-0.61, y=0.5, labels=paste0("1M-2E-3M", "\n", "1M-2M-3E"), pos=1, offset=1, cex=1.5)   
  
  xx2 <- c(0, x[x>=0 & x<=1.22], 1.22)         #Color area
  yy2 <- c(0, dens[x>=0 & x<=1.22], 0) 
  panel.xyarea(xx2, yy2, col="palegreen4", border=NA, lty=1)
  ltext(x=0.61, y=0.5, labels=paste0("1M-2M-3H", "\n", "1M-2H-3M"), pos=1, offset=1, cex=1.5)   
  
  
  panel.xyplot(x, y, ...)
  panel.abline(v=c(-1.22, 0, 1.22), col="blue", lty=2, lwd=2, reference = FALSE)
  
}


# Lattice xyplot
p1 <- xyplot(dens ~ theta,                   
       type = "l",
       main=list(label="Primary Routes", cex=1.5),
       ylim=c(0, 0.5),
       xlab=list(label=expression(theta), cex=1.5),
       ylab=list(label="Density", cex=1.5),
       scales=list(cex=c(1.5, 1.5)),
       panel = pfun.primary
)

# Lattice xyplot
p2 <- xyplot(dens ~ theta,                   
       type = "l",
       main=list(label="Secondary Routes", cex=1.5),
       ylim=c(0, 0.5),
       xlab=list(label=expression(theta), cex=1.5),
       ylab=list(label="Density", cex=1.5),
       scales=list(cex=c(1.5, 1.5)),
       panel = pfun.secondary
)


# print out a plot
grid.arrange(p1, p2, ncol=1)






