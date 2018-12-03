######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 10/14/2018
#Title: Drawing MST panel configurations 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(reshape2)
library(parallel)
library(pbapply)
library(mstR)

##########################################################################
# Read source code
src.dir <- "E:/Dissertation/Analysis/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/Dissertation/Analysis/Temp"
setwd(mydir)

# read item bank and content category information
bank_info <- read.csv("itempar_bank_2PL.csv")
bank_info$g <- runif(nrow(bank_info), 0, 0.4)
names(bank_info) <- c("ID", "PARAM.1", "PARAM.2", "PARAM.3", "CLASS")
cont_require <- read.csv("Content_Requirement.csv")

# transform a content variable to numeric variable
x <- bank_info[, 5]
x <- as.numeric(x)
bank_info$CLASS_RE <- x

# assign other information (e.g., model, score categories) to each item in the bank
bank_info$MODEL <- "3PLM"
bank_info$CATEGORY <- 2

# reorder columns
bank_info <- bank_info[, c(1, 7, 8, 5, 6, 2, 3, 4)]


###-----------------------------------------------------------------------
## 1. 1-3-3 MST panel
# set condtions for ATA
theta <- seq(-8, 8, .1)
RDP <- list(NA, c(-0.44, 0.44), c(-0.44, 0.44))
item.pool <- bank_info
test.length <- 45
size.module <- c(15, 15, 15)
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))
route.map <- read.csv("Routing_map.csv", header=FALSE) # read a routing map
constraints <- list(route.map=route.map, size.module=size.module, content=cont_require[, 2:4])

# create a list of target module information functions (MIFs)
theta.stg1 <- list(c(median(c(-2.0, -0.44)), median(c(-0.44, 0.44)), median(c(0.44, 2.0))))
theta.stg2 <- list(seq(-2.0, -0.44, length.out=5), 
                   seq(-0.44, 0.44, length.out=5),
                   seq(0.44, 2.0, length.out=5))
theta.stg3 <- theta.stg2
target.theta <- list(stage1=theta.stg1, stage2=theta.stg2, stage3=theta.stg3)
targetMIF <- info_mif(item.pool, target.theta, route.map, test.length)

# Assemble a MST
x1 <- ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
               lp.control=list(timeout=30, epsint=0.2, mip.gap=c(0.1, 0.05)))

# set condition to use a function of randomMST in mstR package 
module.mat <- x1$itemvar.df[[1]]
route.scoring="ML"
final.scoring="ML"
range=c(-6, 6)
true.theta <- 0

# set argumetns of randomMST function
start <- list(fixModule=1, D=D)
test <- list(method=route.scoring, D=D, moduleSelect = "MFI")
final <- list(method=final.scoring, range=range, D=D)
itemBank <- cbind(item.pool[,c("PARAM.1", "PARAM.2", "PARAM.3")], d=rep(1, nrow(item.pool)))

# run randomMST
admin.stg3 <- mstR::randomMST(trueTheta=true.theta, itemBank=itemBank, modules=module.mat, 
                         transMatrix=route.map, start=start, test=test, final=final)


###-----------------------------------------------------------------------
## 2. 1-3 MST panel
# set condtions for ATA
theta <- seq(-8, 8, .1)
RDP <- list(NA, c(-0.44, 0.44))
item.pool <- bank_info
test.length <- 30
size.module <- c(15, 15)
D <- 1.702
divide.D <- FALSE
lp.control <- list(timeout=60, epsint=0.2, mip.gap=c(0.1, 0.05))
route.map <- read.csv("Routing_map_2stage.csv", header=FALSE) # read a routing map
constraints <- list(route.map=route.map, size.module=size.module, content=cont_require[, 2:3])

# create a list of target module information functions (MIFs)
theta.stg1 <- list(c(median(c(-2.0, -0.44)), median(c(-0.44, 0.44)), median(c(0.44, 2.0))))
theta.stg2 <- list(seq(-2.0, -0.44, length.out=5), 
                   seq(-0.44, 0.44, length.out=5),
                   seq(0.44, 2.0, length.out=5))
target.theta <- list(stage1=theta.stg1, stage2=theta.stg2)
targetMIF <- info_mif(item.pool, target.theta, route.map, test.length)

# Assemble a MST
x2 <- ata_mstBU(item.pool, constraints, targetMIF, theta=seq(-4, 4, .1), D=1.702, divide.D=FALSE, 
               lp.control=list(timeout=30, epsint=0.2, mip.gap=c(0.1, 0.05)))

# set condition to use a function of randomMST in mstR package 
module.mat <- x2$itemvar.df[[1]]
route.scoring="ML"
final.scoring="ML"
range=c(-6, 6)
true.theta <- 0

# set argumetns of randomMST function
start <- list(fixModule=1, D=D)
test <- list(method=route.scoring, D=D, moduleSelect = "MFI")
final <- list(method=final.scoring, range=range, D=D)
itemBank <- cbind(item.pool[,c("PARAM.1", "PARAM.2", "PARAM.3")], d=rep(1, nrow(item.pool)))

# run randomMST
admin.stg2 <- mstR::randomMST(trueTheta=true.theta, itemBank=itemBank, modules=module.mat, 
                              transMatrix=route.map, start=start, test=test, final=final)


###-----------------------------------------------------------------------
## 3. creat MST panel figues 
name1 <- c("1M", "2E", "2M", "2H")
name2 <- c("1M", "2E", "2M", "2H", "3E", "3M", "3H")

par(mfrow=c(1, 2))
plot(admin.stg2, show.path=FALSE, module.name=name1, cex.m=4, save.plot = FALSE)
plot(admin.stg3, show.path=FALSE, module.name=name2, cex.m=4, save.plot = FALSE)



