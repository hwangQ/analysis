library(lattice)
library(latticeExtra)
library(gridExtra)

## ------------------------------------------------------------------------------
## 1. Plot for an illustration of systematticaly search with 1-3 and 1-3-3 MST
# set theta values and densities
theta <- seq(-4, 4, by = 0.01)            # Data to set up out normal
dens <- dnorm(theta, 0, 1)

# set a panel function 
panel.fun <- function(x, y, ...) {
  
  panel.grid(h=-1, v=-1)  
  
  xx1 <- c(-0.8, x[x>=-0.8 & x<=-0.1], -0.1)         #Color area
  yy1 <- c(0, dens[x>=-0.8 & x<=-0.1], 0) 
  panel.xyarea(xx1, yy1, col="ivory3", border=NA, lty=1)

  xx2 <- c(.1, x[x>=.1 & x<=0.8], 0.8)         #Color area
  yy2 <- c(0, dens[x>=.1 & x<=0.8], 0) 
  panel.xyarea(xx2, yy2, col="ivory3", border=NA, lty=1)

  panel.xyplot(x, y, ...)
  panel.abline(v=c(-0.8, -0.1, .1, .8), col="blue", lty=2, lwd=2, reference = FALSE)
  
}

# Lattice xyplot
p1 <- xyplot(dens ~ theta,                   
            type = "l",
            ylim=c(0, 0.5),
            main=list(label="1-3 and 1-3-3 MSTs", cex=1.5),
            xlab=list(label=expression(theta), cex=1.5),
            ylab=list(label="Density", cex=1.5),
            scales=list(cex=c(1.5, 1.5)),
            panel = panel.fun
)

# print out a plot
print(p1)


## ------------------------------------------------------------------------------
## 2. Plot for an illustration of systematticaly search with 1-2-2 MST
# set theta values and densities
theta <- seq(-4, 4, by = 0.01)            # Data to set up out normal
dens <- dnorm(theta, 0, 1)

# set a panel function 
panel.fun <- function(x, y, ...) {
  
  panel.grid(h=-1, v=-1)  
  
  xx1 <- c(-0.7, x[x>=-0.7 & x<=0.7], 0.7)         #Color area
  yy1 <- c(0, dens[x>=-0.7 & x<=0.7], 0) 
  panel.xyarea(xx1, yy1, col="ivory3", border=NA, lty=1)
  
  panel.xyplot(x, y, ...)
  panel.abline(v=c(-0.7, .7), col="blue", lty=2, lwd=2, reference = FALSE)
  
}

# Lattice xyplot
p2 <- xyplot(dens ~ theta,                   
            type = "l",
            ylim=c(0, 0.5),
            main=list(label="1-2-2 MST", cex=1.5),
            xlab=list(label=expression(theta), cex=1.5),
            ylab=list(label="Density", cex=1.5),
            scales=list(cex=c(1.5, 1.5)),
            panel = panel.fun
)

# print out a plot
print(p2)


# print out a plot
grid.arrange(p1, p2, ncol=1)

