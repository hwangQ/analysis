# a function to extract panel design information from a matrix of routing map

panel_info <- function(route_map) {

  # Arg
  # route_map: a matrix or data.frame including routing information
  
  route_map <- as.data.frame(route_map)
  
  ## 1. find module numbers across stages
  end_col <- ncol(route_map)
  
  # check the modules in the first stage 
  mod_stg1 <- which(colSums(route_map) == 0)
  
  # check the remaining modules across subsequenc stages
  config_info <- list()
  config_info[[1]] <- mod_stg1
  col_1 <- mod_stg1
  i <- 1
  repeat{
    i <- i + 1
    col_2 <- which(colSums(route_map[col_1, ] == 1) > 0)
    config_info[[i]] <- col_2
    if(max(col_2) == end_col) break
    col_1 <- col_2
  }
  names(config_info) <- paste0("stage", 1:length(config_info))
  
  # check the number of stage and modules at each stage
  n.module <- sapply(config_info, length)
  n.stage <- length(n.module)

  ## 2. create a matrix of pathways
  
  # creat a matrix of grid including all possible pathways
  all_grid <- expand.grid(config_info)
  
  # delete any pathways that are not allowed to be used
  remain <- all_grid
  for(s in 1:(n.stage-1)) {
      for(m in 1:n.module[s]) {
          tmp.1 <- config_info[[s]][m]
          tmp.2 <- which(route_map[tmp.1, ] == 1)
          delete <- (remain[, s] == tmp.1 & !(remain[, s + 1] %in% tmp.2))
          remain <- remain[!delete, ]
      }
  }
  
  # sort the remained matrix
  remain <- reshape::sort_df(remain, vars=names(remain))
  remain <- data.matrix(remain)
  rownames(remain) <- 1:nrow(remain)
  colnames(remain) <- paste0("stage", 1:ncol(remain))

  rr <- list(config.info=config_info, pathway=remain, n.module=n.module, n.stage=n.stage)
  rr
  
}