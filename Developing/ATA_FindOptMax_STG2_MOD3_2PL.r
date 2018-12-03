######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 07/25/2018
#Title: Find the optimal values of target maximums for 1-3 MST design structure 
######################################################################################################################################
# call packages
library(lpSolveAPI)
library(SuppDists)
library(reshape2)

##########################################################################
# Read source code
src.dir <- "E:/2018 summer intership/Analysis/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/2018 summer intership/Analysis/Rcode/Temp"
setwd(mydir)

# read item bank and content category information
bank_info <- read.csv("itempar_bank_2PL.csv")
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

# set condtions for ATA
sigma <- seq(1, 2, .1)
theta <- seq(-4, 4, .4)
location <- seq(0.5, 1.5, .1)
item.pool <- bank_info
# n.module <- 2
n.module <- 3
# n.module <- 4
n.stage <- 2
m.size <- c(20, 20)
constraints <- list(n.stage=n.stage, n.module=n.module, size.module=m.size, content=cont_require)
D <- 1.702
divide.D <- FALSE
rel.tol <- c(0.7, -0.5)
lp.control <- list(timeout=40, epsint=0.2, mip.gap=c(0.1, 0.05))

# set directory to save the results
dir4results <- "E:/2018 summer intership/Analysis/Rcode/Result"
dir4save <- file.path(dir4results, paste0("STG", n.stage, "_MOD", n.module[1], "_2PLM")) 

# creat an empty matrix to contain obtimized target maximum values
opt_max <- array(NA, c(length(location), (n.module+1), length(sigma)))
# find optimized target maximum values for all scales of Johnson distributions
for(k in 1:length(sigma)) {
	# find optimized target maximum values for all locations of Johnson distributions
	for(j in 1:length(location)) {

		##################################################################################
		# set target information functions using Johnson distribution
		# (a) set first four moments
		momts_rt <- norm_moments(mu=0, sigma=sigma[k])

		# (b) fit Johnson distribution using the method of mements
		if(n.module[1] == 2) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_2$xi <- location[j]
		}
		
		if(n.module[1] == 3) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- fit_stg2_3 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_3$xi <- location[j]
		}
		
		# (c) find density of Johnson distribution corresponding theta points
		target_rt <- dJohnson(theta, fit_rt)
		target_stg2_1 <- dJohnson(theta, fit_stg2_1)
		target_stg2_2 <- dJohnson(theta, fit_stg2_2)
		if(n.module[1] == 3) target_stg2_3 <- dJohnson(theta, fit_stg2_3)
		
		target_stg2 <- vector('list', n.module[1]) 
		for(i in 1:n.module[1]) target_stg2[[i]] <- get(paste0("target_stg2_", i))
		target.dens <- list(router=target_rt, stage2=target_stg2)
		
		####################################################################################
		# find the optimized target max values sequentialy
		# a) use relative target to find target max value
		x <- tryCatch({ata_mst(item.pool, constraints, target.dens, target.max, relative=TRUE, rel.tol=rel.tol,
					theta, D=D, divide.D=divide.D, lp.control=lp.control)}, 
					error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
		if(is.null(x)) {
			temp.tol <- c(rel.tol[1] + .05, rel.tol[2] - .05)
			i <- 1 
			repeat{
				x <- tryCatch({ata_mst(item.pool, constraints, target.dens, target.max, relative=TRUE, rel.tol=temp.tol,
							theta, D=D, divide.D=divide.D, lp.control=lp.control)},
							error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})
				if(!is.null(x)) break
				temp.tol[1] <- temp.tol[1] + 0.05
				temp.tol[2] <- temp.tol[2] - 0.05
				i <- i + 1
			}
		}
		max_rel <- floor(max(unlist(x$target.sc)))

		# b) first step to find optimized max values using absolute target
		# re-scale the density of Johnson distribution so that maximum of density coressponds to the maximum of target
		if(max_rel < 5) {
			target_mat <- replicate(n=(1 + n.module[1]), expr=seq(max_rel - 2 , max_rel + 2, 0.5), simplify=TRUE)
		} else {
			target_mat <- replicate(n=(1 + n.module[1]), expr=seq(max_rel - 3 , max_rel + 1, 0.5), simplify=TRUE)
		}

		ssd_vec <- c()
		for(r in 1:nrow(target_mat)){
	
			target.max <- list(router=target_mat[r, 1], stage2=c(target_mat[r, 2:(n.module[1] + 1)]))
			ssd <- tryCatch({ata_mst(item.pool=item.pool, constraints=constraints, target.dens=target.dens, target.max=target.max, 
							theta=theta, D=D, divide.D=divide.D, lp.control=lp.control)$ssd}, 
							error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})

			if(is.null(ssd)) ssd <- 99999
			ssd_vec[r] <- ssd
			cat(r, ": ", ssd, "\n", sep="")

		}

		min.num <- which.min(ssd_vec)
		optmax_1 <- as.numeric(target_mat[min.num, ])

		# c) Second step to find optimized max values using absolute target
		# re-scale the density of Johnson distribution so that maximum of density coressponds to the maximum of target
		target_mat <- replicate(n=(1 + n.module[1]), expr=seq((optmax_1[1]-1), (optmax_1[1]+1), .1), simplify=TRUE)

		ssd_vec <- c()
		for(r in 1:nrow(target_mat)){
	
			target.max <- list(router=target_mat[r, 1], stage2=c(target_mat[r, 2:(n.module[1] + 1)]))
			ssd <- tryCatch({ata_mst(item.pool=item.pool, constraints=constraints, target.dens=target.dens, target.max=target.max, 
							theta=theta, D=D, divide.D=divide.D, lp.control=lp.control)$ssd}, 
							error = function(e) {cat("ERROR :", conditionMessage(e), "\n")})

			if(is.null(ssd)) ssd <- 99999
			ssd_vec[r] <- ssd
			cat(r, ": ", ssd, "\n", sep="")

		}

		min.num <- which.min(ssd_vec)
		optmax_2 <- as.numeric(target_mat[min.num, ])

		# store optimized target max values 
		opt_max[j, ,k] <- optmax_2

		cat(k, "th scale and ", j, "th location is done", "\n", sep="")

	}

}

# give dimnames to an array
d1.name <- paste0("DPS of ", location)
d2.name <- c("router", paste0("stage2_m", 1:n.module[1]))
d3.name <- paste0("SCALE of ", location)
dimnames(opt_max) <- list(d1.name, d2.name, d3.name)

# reshape a format of array to data.frame
opt_max_df <- melt(opt_max, varnames=c("DISPERSION", "MODULE", "SCALE"))
opt_max_df <- dcast(opt_max_df, DISPERSION + SCALE ~ MODULE)

# save the optimized target max values
saveRDS(opt_max, file.path(dir4save, paste0("OptMax_2pl_dim2_mod", n.module[1], "_array.rds")))
saveRDS(opt_max_df, file.path(dir4save, paste0("OptMax_2pl_dim2_mod", n.module[1], "_df.rds")))
write.csv(opt_max_df, file.path(dir4save, paste0("OptMax_2pl_dim2_mod", n.module[1], "_df.csv")))

# read optimized target max values 
opt_max <- readRDS(file.path(dir4save, paste0("OptMax_2pl_dim2_mod", n.module[1], "_array.rds")))
opt_max_df <- readRDS(file.path(dir4save, paste0("OptMax_2pl_dim2_mod", n.module[1], "_df.rds")))

##############################################################################################
# creat a plot of TIF
ata_forms <- vector('list', length(sigma))
ata_forms <- lapply(1:length(ata_forms), function(i) vector('list', length(location)))

for(k in 1:length(sigma)) {
	cat(k, "th scale(", sigma[k], ")", "\n", sep="")
	png(file=file.path(dir4save, paste0("SD_", sigma[k], "_1.png")), width=800, height=800)
	par(mfrow=c(3, 2), oma = c(0, 0, 2, 0))
	for(j in 1:6) {

		# (a) set first four moments
		momts_rt <- norm_moments(mu=0, sigma=sigma[k])

		# (b) fit Johnson distribution using the method of mements
		if(n.module[1] == 2) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_2$xi <- location[j]
		}
		
		if(n.module[1] == 3) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- fit_stg2_3 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_3$xi <- location[j]
		}
		
		# (c) find density of Johnson distribution corresponding theta points
		target_rt <- dJohnson(theta, fit_rt)
		target_stg2_1 <- dJohnson(theta, fit_stg2_1)
		target_stg2_2 <- dJohnson(theta, fit_stg2_2)
		if(n.module[1] == 3) target_stg2_3 <- dJohnson(theta, fit_stg2_3)
		
		target_stg2 <- vector('list', n.module[1]) 
		for(i in 1:n.module[1]) target_stg2[[i]] <- get(paste0("target_stg2_", i))
		target.dens <- list(router=target_rt, stage2=target_stg2)

		# set target maximum
		target.max <- list(router=opt_max[j, 1, k], stage2=c(opt_max[j, 2:(n.module[1] + 1), k]))

		# test assembly using the optimized target maximum values
		x <- ata_mst(item.pool=item.pool, constraints=constraints, target.dens=target.dens, target.max=target.max, 
					theta=theta, D=D, divide.D=divide.D, lp.control=lp.control)
		ata_forms[[k]][[j]] <- x$prm.df
		ssd <- round(x$ssd, 3)
	
		# creat a plot
		main <- paste0("SD = ", sigma[k], " & DSP = ", location[j], " & SSD = ", ssd)
		plot(x, xlab="Theta", ylab="TIF", main=main, col=1:4, lwd=2)
	
	}
	mtext(paste0("SD = ", sigma[k]), outer = TRUE, cex = 1.5)
	dev.off()


	png(file=file.path(dir4save, paste0("SD_", sigma[k], "_2.png")), width=800, height=800)
	par(mfrow=c(3, 2), oma = c(0, 0, 2, 0))
	for(j in 7:length(location)) {

		# (a) set first four moments
		momts_rt <- norm_moments(mu=0, sigma=sigma[k])

		# (b) fit Johnson distribution using the method of mements
		if(n.module[1] == 2) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_2$xi <- location[j]
		}
		
		if(n.module[1] == 3) {
			fit_rt <- fit_stg2_1 <- fit_stg2_2 <- fit_stg2_3 <- JohnsonFit(momts_rt, moment="use")
			fit_stg2_1$xi <- -location[j]
			fit_stg2_3$xi <- location[j]
		}
		
		# (c) find density of Johnson distribution corresponding theta points
		target_rt <- dJohnson(theta, fit_rt)
		target_stg2_1 <- dJohnson(theta, fit_stg2_1)
		target_stg2_2 <- dJohnson(theta, fit_stg2_2)
		if(n.module[1] == 3) target_stg2_3 <- dJohnson(theta, fit_stg2_3)
		
		target_stg2 <- vector('list', n.module[1]) 
		for(i in 1:n.module[1]) target_stg2[[i]] <- get(paste0("target_stg2_", i))
		target.dens <- list(router=target_rt, stage2=target_stg2)

		# set target maximum
		target.max <- list(router=opt_max[j, 1, k], stage2=c(opt_max[j, 2:(n.module[1] + 1), k]))

		# test assembly using the optimized target maximum values
		x <- ata_mst(item.pool=item.pool, constraints=constraints, target.dens=target.dens, target.max=target.max, 
					theta=theta, D=D, divide.D=divide.D, lp.control=lp.control)
		ata_forms[[k]][[j]] <- x$prm.df
		ssd <- round(x$ssd, 3)
	
		# creat a plot
		main <- paste0("SD = ", sigma[k], " & DSP = ", location[j], " & SSD = ", ssd)
		plot(x, xlab="Theta", ylab="TIF", main=main, col=1:4, lwd=2)
	
	}
	mtext(paste0("SD = ", sigma[k]), outer = TRUE, cex = 1.5)
	dev.off()

}

# save the assembled test forms to a rds file
saveRDS(ata_forms, file.path(dir4save, paste0("ata_testform_2pl_dim2_mod", n.module[1], ".rds")))




