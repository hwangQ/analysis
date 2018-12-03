# number of examinees
nstd <- c(1000, 5000)

# set seed
set.seed(123)

# generate true abilties
theta_1000 <- stats::rnorm(nstd[1], 0, 1)
theta_1000 <- ifelse(theta_1000 > -4, theta_1000, -4)
theta_1000 <- ifelse(theta_1000 < 4, theta_1000, 4)

theta_5000 <- stats::rnorm(nstd[2], 0, 1)
theta_5000 <- ifelse(theta_5000 > -4, theta_5000, -4)
theta_5000 <- ifelse(theta_5000 < 4, theta_5000, 4)

# save results
saveRDS(theta_1000, file="Input/theta_1000.rds")
saveRDS(theta_5000, file="Input/theta_5000.rds")

