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
n.module <- c(3, 3)
n.stage <- 3
D <- 1.702

# read a routing map
route_map <- read.csv("Routing_map.csv", header=FALSE)
panel_info(route_map)


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

# find points where two adjacent TIFs intersect across all stages 
module.form <- ata_forms
RDP <- vector('list', 3)
RDP[[1]] <- NA
RDP[[2]] <- c(-.44, .44)
RDP[[3]] <- c(-.44, .44)
cut.score <- cutoff(module.form, route.map, RDP, D, range.theta, interval)

# conditional raw score distribution of each module given an ability value
cond_dist <- ability_dist_sep(module.form, theta, route.map, D)

# conditional ability distribution of total-test ability estimates given a true ability
x <- cond_dist
eos.path <- eos_list$eos_path
n.path <- eos_list$n_path
joint.dist <- ability_dist(x, eos.path, n.path, cut.score)


# calculate a mean and sd of conditional ability distribution at each ability value
cond_moments <- sapply(joint.dist$joint.dist$stage3, cal_moments, node=eos_list$eos_path$stage3)
cond_sigma <- sqrt(cond_moments[2, ])

# marginal reliability
w <- dnorm(theta, 0, 1)
w <- w/sum(w)
margin_rel <- mrel(var.pop=1, var.cond=cond_moments[2, ], wieghted=TRUE, w=w)

# average of CSEM within an interval of ability
selected <- colnames(cond_moments) %in% as.character(seq(-4, 4, .1))
mean.csem <- mean(cond_sigma[selected])

# maximum value of CSEM between -2 and 2 
max.csem <- max(cond_sigma[selected])


