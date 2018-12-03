# a function to generate item parameters  
gen_ipar <- function(nitem, type=c("dc", "py"), ThreePL=TRUE,
					a=list(params=c(0, 0.5), range=c(0.8, 2.7)), 
					b=list(params=c(0, 1), range=c(-2.5, 2.5)), 
					g=list(params=c(3, 17), range=c(0.0, 0.3)),
					d=list(ncats=5, params1=c(-0.2, 0), params2=c(0.4, 0.7))
					) {

	# generate dichotomous item parameters
	if(type==tolower("dc")) {
		rr <- gen_ipar_dc(nitem=nitem, ThreePL=ThreePL, a=a, b=b, g=g) 
	}
	
	# generate polytomous item parameters
	if(type==tolower("py")) {
		rr <- gen_ipar_py(nitem=nitem, a=a, d=d) 
	}

	rr

}

# a function to generate item parameters: dichotomous items
gen_ipar_dc <- function(nitem, ThreePL=TRUE,
					a=list(params=c(0, 0.5), range=c(0.8, 2.7)), 
					b=list(params=c(0, 1), range=c(-2.5, 2.5)), 
					g=list(params=c(3, 17), range=c(0.0, 0.3))
					) {
					
	# Create dichotomous item parameters
	repeat{
		x <- rlnorm(nitem, meanlog=a$params[1], sdlog=a$params[2])
		if(min(x) >= a$range[1] & max(x) <= a$range[2]) {
			a.dc <- x
			break
		} else {next}
	}

	repeat{
		x <- rnorm(nitem, mean=b$params[1], sd=b$params[2])
		if(min(x) >= b$range[1] & max(x) <= b$range[2]) {
			b.dc <- x
			break
		} else {next}
	}

	if(ThreePL) {
		repeat{
			x <- rbeta(nitem, shape1=g$params[1], shape2=g$params[2])
			if(min(x) >= g$range[1] & max(x) <= g$range[2]) {
				g.dc <- x
				break
			} else {next}
		}
	} else {
		g.dc <- rep(0, nitem)
	}

	prm_mat <- data.frame(a=a.dc, b=b.dc, g=g.dc)
	
	rr <- list(a=a.dc, b=b.dc, g=g.dc, prm_mat=prm_mat)
	rr
	
}

# a function to generate item parameters: polytomous items
gen_ipar_py <- function(nitem, 
					a=list(params=c(0, 0.5), range=c(0.8, 2.7)), 
					d=list(ncats, params1=c(-0.2, 0), params2=c(0.4, 0.7))
					) {

	# Create polytomous item parameters
	repeat{
		x <- rlnorm(nitem, meanlog=a$params[1], sdlog=a$params[2])
		if(min(x) >= a$range[1] & max(x) <= a$range[2]) {
			a.py <- x
			break
		} else {next}
	}

	d.py <- vector('list', nitem)
	for(i in 1:nitem) {
		x <- runif(1, min=d$params1[1], max=d$params1[2])
		y <- runif((d$ncats-2), min=d$params2[1], max=d$params2[2])
		xy <- cumsum(c(x, y))
		d.py[[i]] <- xy
	}

	rr <- list(a=a.py, d=d.py)
	rr

}