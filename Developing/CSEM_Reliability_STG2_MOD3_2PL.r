######################################################################################################################################
#Authors: Hwanggyu Lim (UMASS AMHERST)
#Date: 07/25/2018
#Title: Calculate Conditional Standard Error of Measurement and Find optimized Design 
#       based on Three objective functions for 1-3 MST design structure 
######################################################################################################################################
# Read source code
src.dir <- "E:/2018 summer intership/Analysis/Rcode/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set directory
mydir <- "E:/2018 summer intership/Analysis/Rcode/Temp"
setwd(mydir)

# set condtions
sigma <- seq(1, 2, .1)
location <- seq(0.5, 1.5, .1)
theta <- seq(-4, 4, 0.1)
range.theta <- c(-4, 4)
interval <- c(-2.5, 2.5)
# nmod <- 2
nmod <- 3
# nmod <- 4
nstg <- 2
m.size <- c(20, 20)
D <- 1.702

# set directory to save the results
dir4results <- "E:/2018 summer intership/Analysis/Rcode/Result"
dir4save <- file.path(dir4results, paste0("STG", nstg, "_MOD", nmod[1], "_2PLM")) 

# read the assembled test forms
ata_forms <- readRDS(file.path(dir4save, paste0("ata_testform_2pl_dim2_mod", nmod[1], ".rds")))

##########################################################################
# calculate CSEM and marginal reliability
csem <- array(NA, c(length(theta), length(location), length(sigma))) # to contain CSEM
mgrel <- array(NA, c(length(sigma), length(location))) # to contain the marginal reliability
unw_mgrel_1 <- array(NA, c(length(sigma), length(location))) # to contain the unweighted marginal riliability 
unw_mgrel_2 <- array(NA, c(length(sigma), length(location))) # to contain the unweighted marginal riliability given the ability from -2 and 2
csem_max <- array(NA, c(length(sigma), length(location))) # to contain the maximum CSEM
cutScore <- array(NA, c(length(location), (nmod[1]-1), length(sigma))) # to contain cut scores

# give dimnames
ab.name <- paste0("THETA.", round(theta, 3))
loc.name <- paste0("DSP.", location)
sd.name <- paste0("SD.", sigma)
dimnames(csem) <- list(ab.name, loc.name, sd.name)
dimnames(mgrel) <-list(sd.name, loc.name)
dimnames(unw_mgrel_1) <-list(sd.name, loc.name)
dimnames(unw_mgrel_2) <-list(sd.name, loc.name)
dimnames(csem_max) <- list(sd.name, loc.name)
dimnames(cutScore) <- list(loc.name, paste0("CUT.", 1:(nmod[1]-1)), sd.name)

for(k in 1:length(sigma)) {
	for(j in 1:length(location)) {

		# create a data.frame of item parameters
		# (a) for a router and each module of stage 2
		df_rt <- ata_forms[[k]][[j]][[1]]
		for(i in 1:nmod[1]) {
			assign(paste0("df_stg2_", i), ata_forms[[k]][[j]][[i+1]])
		}

		# (b) for each path
		for(i in 1:nmod[1]) {
			assign(paste0("df_path_", i), rbind(df_rt, get(paste0("df_stg2_", i))))
		}

		# estimate the equated observed scores corresponding to observed raw score
		df_path <- lapply(1:nmod[1], function(i) get(paste0("df_path_", i)))
		x <- list(router=df_rt, path=df_path)
		eos_list <- est_eos(x, range.theta=range.theta, D=D, constraint=TRUE)
		eos_rt <- eos_list$eos_router
		eos_path_mat <- eos_list$eos_path

		# find points where two adjacent TIFs cross 
		x <- df_stg2 <- lapply(1:nmod[1], function(i) get(paste0("df_stg2_", i)))
		info_t <- cross_info(x, range.theta=range.theta, D=D, interval=interval, n=1000)
	
		cut.score <- info_t$cut.score
		cutScore[j, ,k] <- cut.score
		info_modules <- info_t$testInfo

		# conditional raw score distribution of each module given an ability value
		x <- list(router=df_rt, stage2=df_stg2)
		cond_dist <- ability_dist_sep(x, theta=theta, nstg=nstg, nmod=nmod, D=D)

		# conditional ability distribution of total-test ability estimates given a true ability
		x <- cond_dist
		theta_condist <- ability_dist(x, eos.router=eos_rt, cut.score=cut.score)

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
	
		cat(k, "th scale & ", j, "th location", "\n", sep="")

	}
}

##########################################################################
# check optimized value among three objective functions
mgrel == max(mgrel)
unw_mgrel_2 == max(unw_mgrel_2)
csem_max == min(csem_max)

# save results
write.csv(mgrel, file.path(dir4save, "Marinal_Reliability.csv"))
write.csv(unw_mgrel_2, file.path(dir4save, "Unweighted_Marinal_Reliability.csv"))
write.csv(csem_max, file.path(dir4save, "Maximum_CSEM.csv"))


##########################################################################
# plot of CSEM
plot(csem[, 1, 1] ~ theta, type="l", lwd=2, col=1, lty=1,
     ylim=c(0.0, 1.0),
     xlim=range.theta,
     ylab="CSEM", 
     xlab="Proficiency",
     main="CSEM"
	 )
for(j in 2:length(location)) { 
	lines(csem[, j, 1] ~ theta, type="l", lwd=2, col=j, lty=j)
}
abline(v=c(-2, 2), lty=2)
legend("top", legend=paste0("DSP = ", location[1:j]), lwd=3,
		lty=1:length(location), col=1:length(location), cex=2)

# plot of marginal reliability and maximum CSEM
par(mfrow=c(2, 2))
# plot of marginal reliability
matplot(x=location, y=t(mgrel), type='b', lwd=2, 
		lty=1:length(sigma), 
		col=1:length(sigma), 
		pch=1:length(sigma),
		ylim=c(0.8, 1),
		ylab="Marginal Reliability", 
		xlab="Dispersion of location",
		main=paste0("Marginal Reliability")
		)
legend("top", legend=paste0("SD = ", sigma), ncol=4, lwd=2,
	lty=1:length(sigma), col=1:length(sigma), pch=1:length(sigma), cex=1.3)
	
# plot of unweighted marginal reliability
matplot(x=location, y=t(unw_mgrel_2), type='b', lwd=2, 
		lty=1:length(sigma), 
		col=1:length(sigma), 
		pch=1:length(sigma),
		ylim=c(0.8, 1),
		ylab="Unweighted Marginal Reliability", 
		xlab="Dispersion of location",
		main=paste0("Unweighted Marginal Reliability between -2 and 2")
		)
legend("top", legend=paste0("SD = ", sigma), ncol=4, lwd=2,
	lty=1:length(sigma), col=1:length(sigma), pch=1:length(sigma), cex=1.3)
		
# plot of maximum CSEM between -2 and 2
matplot(x=location, y=t(csem_max), type='b', lwd=2, 
		lty=1:length(sigma), 
		col=1:length(sigma), 
		pch=1:length(sigma),
		ylim=c(0.3, 0.7),
		ylab="CSEM", 
		xlab="Dispersion of location",
		main=paste0("Maximum of CSEM between -2 and 2")
		)
legend("top", legend=paste0("SD = ", sigma), ncol=4, lwd=2,
	lty=1:length(sigma), col=1:length(sigma), pch=1:length(sigma), cex=1.3)

	






