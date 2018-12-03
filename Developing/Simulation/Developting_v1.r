# set directory
mydir <- "E:/Dissertation/Analysis/Temp"
setwd(mydir)

# Read source code
src.dir <- "E:/Dissertation/Analysis/R/"
src.files <- list.files(src.dir, pattern="*.r")
for(files in src.files) source(paste0(src.dir, files))

# set condtions
theta <- seq(-4, 4, 0.1)
range.theta <- c(-4, 4)
interval <- c(-2.5, 2.5)
n.module <- c(3, 3)
n.stage <- 3
D <- 1.702

# read a routing map
route_map <- read.csv("Routing_map.csv", header=FALSE)
panel.info <- panel_info(route_map)
config.info <- panel.info$config

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
cut.score <- vector('list', n.stage)
names(cut.score) <- paste0("stage", 1:n.stage)
for(i in 1:n.stage) {
    x <- df_forms[[i]]
    if(length(x) == 1) {
       cut.score[[i]] <- NA
    } else {
       cut.score[[i]] <- cross_info(x, RDP=RDP, range.theta=range.theta, D=D, interval=interval, n=1000)$cut.score
    }
}


th <- 2
resp <- NA
modules <- 1
cum.params <- data.frame()

eos.path <- eos_list$eos_path

panel.info <- panel_info(route_map)
x <- ata_forms
theta <- 0


runMST_indv(x, theta, panel.info, cut.score, D, range.theta, ini.module=1)

runMST_indv <- function(x, theta, panel.info, cut.score, eos.path=NULL, ini.module=1, D, range.theta, ...) {

   # Arguments
   # x: a list including a data.frame of item parameter informations of modules across all stages
   # theta: a scala of true ability value
   # panel.info: an object obtained from "panel_info" function. This includes all informations regarding a panel configuration
   # cut.score: a list of cut scores across all stages
   # eos.path: a list of the equated observed (or number correct) scores for every paths across all stages
   # ini.module: a scala indicating a module to be chosen at the first stage.

   response <- NULL
   modules <- ini.module
   df_params <- data.frame()
   
   # extract panel information
   pathway <- panel.info$pathway
   config.info <- panel.info$config
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # simulate item responses and estimate an examinee ability
   sel.mod <- ini.module
   th <- theta
   for(s in 1:n.stage) {
      
       # find item parameter informaion corresponding to the selected module
       params <- x[[sel.mod]]
       a <- params$PARAM.1
       b <- params$PARAM.2
       g <- params$PARAM.3
       cats <- params$CATEGORY

       # simulate item responses and calculate a summed score for the selected module 
       response <- c(response, simdat(theta=th, a.dc=a, b.dc=b, g.dc=g, cats=cats, D=D))
       sum.score <- sum(response)
      
       # estimate ability parameters using the inverse TCC method 
       df_params <- rbind(df_params, params)
       
       if(!eos.path == NULL) {
          unique.path <- data.matrix(unique(pathway[, 1:s]))
          idx.logic <- sapply(1:nrow(unique.path), function(k) all(modules == unique.path[k, ]))
          thetas <- eos[[s]][, idx.logic]
       } else {
          thetas <- convert_eos(df_params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
       }
       th <- thetas[sum.score + 1]

       # find a next module to route an examinee
       if(s < n.stage) {
          selected <- give_path(th, cut.score[[s+1]])$path
          sel.mod <- config.info[[s+1]][selected]
          modules <- as.numeric(c(modules, sel.mod))
       }
    
   }
   
   # get results
   rr <- list(est.theta=th, true.theta=theta, sum.score=sum.score, modules=modules, response=response, df_params=df_params)
   rr
   
}






runMST_indv <- function(x, theta, panel_info, cut.score, D, range.theta, ini.module=1, ...) {


   response <- NULL
   modules <- ini.module
   df_params <- data.frame()
   
   # extract panel information
   pathway <- panel.info$pathway
   config.info <- panel.info$config
   n.stage <- panel.info$n.stage
   n.module <- panel.info$n.module
   
   # simulate item responses and estimate an examinee ability
   sel.mod <- ini.module
   th <- theta
   for(s in 1:n.stage) {
      
       # find item parameter informaion corresponding to the selected module
       params <- x[[sel.mod]]
       a <- params$PARAM.1
       b <- params$PARAM.2
       g <- params$PARAM.3
       cats <- params$CATEGORY

       # simulate item responses and calculate a summed score for the selected module 
       response <- c(response, simdat(theta=th, a.dc=a, b.dc=b, g.dc=g, cats=cats, D=D))
       sum.score <- sum(response)
      
       # estimate ability parameters using the inverse TCC method 
       df_params <- rbind(df_params, params)
       thetas <- convert_eos(df_params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
       th <- thetas[sum.score + 1]

       # find a next module to route an examinee
       if(s < n.stage) {
          selected <- give_path(th, cut.score[[s+1]])$path
          sel.mod <- config.info[[s+1]][selected]
          modules <- as.numeric(c(modules, sel.mod))
       }
    
   }
   
   # get results
   rr <- list(est.theta=th, true.theta=theta, sum.score=sum.score, modules=modules, response=response, df_params=df_params)
   rr
   
}















































# stage1
params <- ata_forms[[1]]
cum.params <- rbind(cum.params, params)
a <- params$PARAM.1
b <- params$PARAM.2
g <- params$PARAM.3
cats <- params$CATEGORY

resp <- simdat(theta=th, a.dc=a, b.dc=b, g.dc=g, cats=cats, D=D)
sum.score <- sum(resp)
thetas <- convert_eos(cum.params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
est.th <- thetas[sum.score + 1]

selected <- give_path(est.th, cut.score[[1]])$path
sel.mod <- config.info[[2]][selected]
modules <- c(modules, sel.mod)

# stage2
params <- ata_forms[[sel.mod]]
cum.params <- rbind(cum.params, params)
a <- params$PARAM.1
b <- params$PARAM.2
g <- params$PARAM.3
cats <- params$CATEGORY

resp <- c(resp, simdat(theta=est.th, a.d=a, b.dc=b, g.dc=g, cats=cats, D=D))
sum.score <- sum(resp)
thetas <- convert_eos(cum.params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
est.th <- thetas[sum.score + 1]

selected <- give_path(est.th, cut.score[[2]])$path
sel.mod <- config.info[[3]][selected]
modules <- c(modules, sel.mod)

# stage3
params <- ata_forms[[sel.mod]]
cum.params <- rbind(cum.params, params)
a <- params$PARAM.1
b <- params$PARAM.2
g <- params$PARAM.3
cats <- params$CATEGORY

resp <- c(resp, simdat(theta=est.th, a.d=a, b.dc=b, g.dc=g, cats=cats, D=D))
sum.score <- sum(resp)
thetas <- convert_eos(cum.params, D=D, range.theta=range.theta, intpo=TRUE, constraint=TRUE)$theta
est.th <- thetas[sum.score + 1]

est.th
sum.score






