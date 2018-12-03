# A function to create a data.frame of item parameters to be used in the analysis

shape_df <- function(par.dc=list(a=NULL, b=NULL, g=NULL), par.py=list(a=NULL, d=NULL), item.id=NULL, cats, model) {

	nitem <- sum(length(par.dc$b), length(par.py$d))
	max.cat <- max(cats)
	if(is.null(item.id)) item.id <- paste0("V", 1:nitem)
	if(length(cats) == 1) cats <- rep(cats, nitem)
	if(length(model) == 1) model <- rep(toupper(model), nitem)
	
	any.dc <- any(cats == 2)
	any.py <- any(cats > 2)
	
	if(all(cats == 2)) {
		prm_mat <- array(NA, c(nitem, 3))
	} else {
		prm_mat <- array(NA, c(nitem, max.cat))
	}
	
	# if there are dichotomous items
	if(any.dc) {
		if(is.null(par.dc$g)) par.dc$g <- 0
		prm_mat[cats == 2, 1] <- par.dc$a
		prm_mat[cats == 2, 2] <- par.dc$b
		prm_mat[cats == 2, 3] <- par.dc$g
	}
	
	# if there are polytomous items
	if(any.py) {
		row.py <- which(cats > 2)
		prm_mat[row.py, 1] <- par.py$a
		for(i in 1:length(row.py)) {
			prm_mat[row.py[i], 2:(length(par.py$d[[i]])+1)] <- par.py$d[[i]]
		}
	}
	
	if(all(cats == 2)) {
		colnames(prm_mat) <- paste0("PARAM.", 1:3)
	} else {
		colnames(prm_mat) <- paste0("PARAM.", 1:max.cat)
	}
	
	full_df <- data.frame(ID=item.id, FACTOR=rep(1, length(cats)), CATEGORY=cats, MODEL=model, prm_mat)

	full_df

}




