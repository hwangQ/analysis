## Title : Create item pools 
library(tidyverse)

## Study 1 --------------------------------------------
# set seed
set.seed(135)

# set the number of items
nitem <- 300

# generate item parameters
repeat{
  x <- rlnorm(nitem, meanlog=0, sdlog=0.3)
  if(min(x) > 0.5 & max(x) < 2.5) {
    a <- x
    break
  } else {next}
}
repeat{
  x <- rnorm(nitem, mean=0, sd=1)
  if(min(x) > -3.0 & max(x) < 3.0) {
    b <- x
    break
  } else {next}
}
g <- rbeta(n=nitem, shape1=5, shape2=42)

# generate content categories
cats <- c(rep(1, nitem * .25), rep(2, nitem * .25), rep(3, nitem * .25), rep(4, nitem * .25))
class <- sample(x=cats, size=nitem)

# create a data.frame
pool_1 <- data.frame(a=a, b=b, g=g, class=class)

# descriptive statistics of item pool
mu <- round(apply(pool_1, 2, mean), 2)
sigma <- round(apply(pool_1, 2, sd), 2)
Min <- round(apply(pool_1, 2, min), 2)
Max <- round(apply(pool_1, 2, max), 2)
stat_df1 <- data.frame(mu=mu, sigma=sigma, Min=Min, Max=Max)

# export data
saveRDS(pool_1, "Input/item_pool_1.rds")
write.csv(stat_df1, "Input/stat_pool_1.csv")


## Study 2 --------------------------------------------
# (1) item pool with 200 items 
# set seed
set.seed(377)

# set the number of items
nitem <- 200

# read MAPT Math data set
mapt_dat <- read.csv("Input/Math2018_workbook.csv")

# filter out items with sample size greater than 500
dat <- mapt_dat %>% 
  filter(!is.na(a_FY18)) %>% 
  filter(N.count.in.2018 >= 500) %>% 
  select(a_FY18, b_FY18 , c_FY18)

# randomly select 400 items for the item pool
idx <- sort(sample(1:nrow(dat), nitem, FALSE))
pool_2 <- dat[idx, ]

# generate content categories
cats1 <- c(rep(1, nitem * .25), rep(2, nitem * .25), rep(3, nitem * .25), rep(4, nitem * .25))
class1 <- sample(x=cats1, size=nitem)
cats2 <- c(rep(1, nitem * .25), rep(2, nitem * .25), rep(3, nitem * .50))
class2 <- sample(x=cats2, size=nitem)

# create a data.frame
pool_2 <- data.frame(pool_2, class1=class1, class2=class2)
names(pool_2) <- c("a", "b", "c", "class1", "class2")

# descriptive statistics of item pool
mu <- round(apply(pool_2, 2, mean), 2)
sigma <- round(apply(pool_2, 2, sd), 2)
Min <- round(apply(pool_2, 2, min), 2)
Max <- round(apply(pool_2, 2, max), 2)
stat_df2 <- data.frame(mu=mu, sigma=sigma, Min=Min, Max=Max)

# export data
saveRDS(pool_2, "Input/item_pool_2_B200.rds")
write.csv(stat_df2, "Input/stat_pool_2_B200.csv")

## Study 2 --------------------------------------------
# (1) item pool with 400 items 
# set seed
set.seed(771)

# set the number of items
nitem <- 400

# read MAPT Math data set
mapt_dat <- read.csv("Input/Math2018_workbook.csv")

# filter out items with sample size greater than 500
dat <- mapt_dat %>% 
  filter(!is.na(a_FY18)) %>% 
  filter(N.count.in.2018 >= 500) %>% 
  select(a_FY18, b_FY18 , c_FY18)

# randomly select 400 items for the item pool
idx <- sort(sample(1:nrow(dat), nitem, FALSE))
pool_2 <- dat[idx, ]

# generate content categories
cats1 <- c(rep(1, nitem * .25), rep(2, nitem * .25), rep(3, nitem * .25), rep(4, nitem * .25))
class1 <- sample(x=cats1, size=nitem)
cats2 <- c(rep(1, nitem * .25), rep(2, nitem * .25), rep(3, nitem * .50))
class2 <- sample(x=cats2, size=nitem)

# create a data.frame
pool_2 <- data.frame(pool_2, class1=class1, class2=class2)
names(pool_2) <- c("a", "b", "c", "class1", "class2")

# descriptive statistics of item pool
mu <- round(apply(pool_2, 2, mean), 2)
sigma <- round(apply(pool_2, 2, sd), 2)
Min <- round(apply(pool_2, 2, min), 2)
Max <- round(apply(pool_2, 2, max), 2)
stat_df2 <- data.frame(mu=mu, sigma=sigma, Min=Min, Max=Max)

# export data
saveRDS(pool_2, "Input/item_pool_2_B400.rds")
write.csv(stat_df2, "Input/stat_pool_2_B400.csv")






