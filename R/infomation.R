# item information function for dichotomous data
info.dich <- function(theta, a, b, g, D) {

   p <- pl.fn(theta=theta, a=a, b=b, g=g, D=D)
   q <- 1 - p 
   info <- D^2 * a^2 * (q/p * ((p - g)/(1 - g))^2)
   info
   
}

# item information function for polytomous data
info.poly <- function(theta, a, d, D, model) {

   model <- toupper(model)
   if(!model %in% c("GRM", "GPCM")) stop("Specify a correct model for polytomous items")
   info <- switch(model,
              GRM = info.grm(theta, a, d, D), 
              GPCM = info.gpcm(theta, a, d, D))

   info

}

# test information function
test.info <- function(x, theta, D) {

   any.dc <- any(x$CATEGORY == 2)
   any.py <- any(x$CATEGORY > 2)
   id <- x$ID
   prm_list <- paramList(x)

   # if there are dichotomous items
   if(any.dc) {
      a <- prm_list$dicho.item$a
      b <- prm_list$dicho.item$b
      g <- prm_list$dicho.item$g
      
      # item information matrix
      infomat_dc <- array(NA, c(length(a), length(theta)))
      for(i in 1:length(a)) {
          infomat_dc[i, ] <- info.dich(theta=theta, a=a[i], b=b[i], g=g[i], D=D)
      }
   } else {
      infomat_dc <- NULL
   }

   # if there are polytomous items
   if(any.py) {
      pModel <- prm_list$poly.item$pModel
      aa <- prm_list$poly.item$a
      d <- prm_list$poly.item$d

      # item information matrix
      infomat_py <- array(NA, c(length(aa), length(theta)))
      for(i in 1:length(aa)) {
          infomat_py[i, ] <- info.poly(theta=theta, a=aa[i], d=d[[i]], D=D, model=pModel[i])
      }
   } else {
      infomat_py <- NULL
   }
   
   # creat a item infomation matrix for all items
   infomat <- rbind(infomat_dc, infomat_py)
   
   # re-order the item information maxtirx along with the original order of items
   pos <- c(prm_list$dicho.item$dPos, prm_list$poly.item$pPos)
   if(length(pos) > 1) {
      infomat <- cbind(infomat[order(pos), ])
   } else {
      infomat <- infomat
   }
   rownames(infomat) <- id
   colnames(infomat) <- paste0("theta.", 1:length(theta))
   
   # create a vector for test infomation 
   testInfo <- colSums(infomat)
   
   rr <- list(itemInfo=infomat, testInfo=testInfo)
   rr

}


# a function to calculate the second derivatives of GPCM for theta
gpcm.deriv2 <- function (theta, a, d, D){

   ### Important Note
   ### Under the parameterization of the GPCM in this code
   ### I need to get rid of the base category which has a 
   ### pseudo step of 0
   ### Also, theta should be a scalar here. 

   j <- 1:(length(d)-1)
   m <- length(d)
   d <- d[2:m]
   Da <- D * a
   kernal <- exp(cumsum(Da * (theta - d)))
   hess <- Da^2 * (sum(j^2 * kernal)/(1 + sum(kernal)) - (sum(j * kernal)/(1 + sum(kernal)))^2)

   hess

}

# item information function for GPCM
info.gpcm <- function(theta, a, d, D){

   ### Here, theta can be a vector whereas theta should ba a scalar in gpcm.deriv2
   info <- sapply(1:length(theta), function(i) gpcm.deriv2(theta[i], a, d, D))
   info

}

# item information function for GRM
info.grm <- function(theta, a, d, D) {

   theta <- as.numeric(theta)
   n <- length(as.numeric(theta))

   info <- 0 
   for(i in 1:n) {
      j <- 1:(length(d)+1)
      m <- length(j)
      ps <- pl.fn(theta=theta[i], a=a, b=d, g=0, D=D)
      ps <- c(1, ps, 0)
      qs <- 1 - ps
      info.j <- 0
      for(j in 2:(m + 1)) {
         info.j <- info.j + (a * ps[j - 1] * qs[j - 1] - a * ps[j] * qs[j])^2 / (ps[j - 1] - ps[j])
      }
      info[i] <- info.j
   }

   info

}

