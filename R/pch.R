#' @importFrom stats model.matrix model.response delete.response model.frame terms
#' @importFrom stats model.weights logLik coef nobs vcov .getXlevels quantile runif approxfun sd
#' @importFrom survival Surv is.Surv



#' @export
pchreg <- function(formula, breaks, data, weights){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	if(!survival::is.Surv(zyd <- model.response(mf)))
		{stop("the model response must be created with Surv()")}
	if((n <- nrow(zyd)) == 0){stop("zero non-NA cases")}
	type <- attributes(zyd)$type
	if(type == "right"){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]}
	else if(type == "counting"){z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]}
	else{stop("only type = 'right' and type = 'counting' are supported")}
	if(!(any(d == 1))){stop("all observation are censored")}
	
	x <- model.matrix(mt,mf)
	w <- model.weights(mf)
	if(is.null(w)){w <- rep.int(1,n)}

	###

	fit <- pch.fit(z,y,d,x,w,breaks)

	###
  
	attr(mf, "u") <- fit$u
	attr(mf, "y") <- fit$y
	attr(mf, "type") <- type
	attr(mf, "n") <- length(y)
	attr(mf, "n.events") <- sum(d)

	beta <- fit$beta
	r <- dim(beta)
	colnames(beta) <- 1:r[2]
	vn <- paste(rep.int(rownames(beta), r[2]), rep(1:r[2], each = r[1]), sep = "_")
	dimnames(fit$covar) <- list(vn, vn)

  
	Hy <- predF.pch(fit, x, y)
	Hz <- (if(type == "counting") predF.pch(fit,x,z)[,"Haz"] else 0)
	logLik <- sum(d*log(Hy[,"haz"]), na.rm = TRUE) - sum(Hy[,"Haz"] - Hz)
	# note: 'haz' can be 0, I set 0*log(0) = 0.
	attr(logLik, "df") <- sum(fit$beta != 0)
	
	fit <- list(call = cl, beta = beta, breaks = fit$breaks, 
		covar = fit$covar, logLik = logLik, 
		lambda = fit$lambda, Lambda = fit$Lambda, mf = mf)
	class(fit) <- "pch"
	fit
}


#' @export
predict.pch <- function(object, type = c("distr", "quantile", "sim"), 
	newdata, p, sim.method = c("quantile", "sample"), ...){
  
	if(is.na(match(type <- type[1], c("distr", "d", "quantile", "q", "sim", "s"))))
		{stop("invalid 'type'")}
	type <- strsplit(type, split = "")[[1]][1]
	if(type == "s"){
		if(is.na(match(method <- sim.method[1], c("quantile", "q", "sample", "s"))))
		{stop("invalid 'method'")}
		method <- strsplit(method, split = "")[[1]][1]
	}
	if(type == "q"){
	  if(missing(p)){stop("please indicate 'p'")}
	  if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}	
	}
	
	mf <- object$mf
	mt <- terms(mf)
	xlev <- .getXlevels(mt, mf)
	fittype <- attr(mf, "type")
	obj <- list(beta = object$beta, breaks = object$breaks, y = attr(mf, "y"), 
		lambda = object$lambda, Lambda = object$Lambda, u = attr(mf, "u"))

	if(missing(newdata)){
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	}
	else{
		if(type == "d"){
			yn <- as.character(if(fittype == "counting") mt[[2]][[3]] else mt[[2]][[2]])
			if(is.na(ind <- match(yn, colnames(newdata))))
			{stop("for 'type = distr', 'newdata' must contain the y-variable")}
			if(fittype == "right" && length(mt[[2]]) == 3)
			  {newdata[,as.character(mt[[2]][[3]])] <- 1}
			if(fittype == "counting"){
				newdata[,as.character(mt[[2]][[4]])] <- 1
				newdata[,as.character(mt[[2]][[2]])] <- -Inf
				newdata[,yn] <- pmax(newdata[,yn], obj$breaks[1] - 1)
			}
		}
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
		{stop("'newdata' must contain all x-variables")}
		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){stop("zero non-missing values in the supplied newdata")}

		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
	}

	x <- model.matrix(mt, mf)
	n <- length(miss) + length(nomiss)
	if(type == "d"){
		y <- cbind(model.response(mf))[,1 + (fittype == "counting")]
		out <- matrix(, n, 4)
		pred <- predF.pch(obj,x,y)
		out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
  	else if(type == "q"){
		out <- matrix(, n, length(p))
		pred <- predQ.pch(obj,x,p)
        	out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
	else{
		tsim <- rep.int(NA,n)
		tsim[nomiss] <- sim.pch(obj,x,method)
		return(tsim)
	}
}

# print and summary method
#' @export
print.pch <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
summary.pch <- function(object, ...){
  out <- list(call = object$call, 
    n = attr(object$mf, "n"), n.events = attr(object$mf, "n.events"), 
		n.free.par = attr(object$logLik, "df"), logLik = object$logLik)
  class(out) <- "summary.pch"
  out
}

#' @export
print.summary.pch <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of obs: ", paste(deparse(round(x$n)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of events: ", paste(deparse(round(x$n.events)), sep = " ", collapse = " "), "\n", sep = "")
  cat("n. of free parameters: ", paste(deparse(round(x$n.free.par)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", paste(deparse(round(as.numeric(x$logLik),1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Use predict() to obtain predictions from the fitted model")
  cat("\n")
}
#' @export
logLik.pch <- function(object, ...){object$logLik}
#' @export
nobs.pch <- function(object, ...){nrow(object$mf)}
#' @export
coef.pch <- function(object, ...){object$beta}
#' @export
vcov.pch <- function(object, ...){object$covar}

CumSum <- function(x){
	if((q <- ncol(x)) > 1){
		for(j in 2:q){x[,j] <- x[,j] + x[,j-1]}
	}
	x
}


makebreaks <- function(y, breaks){
  
	n <- length(y)
	r <- range(y)
	if(missing(breaks)){breaks <- max(2, min(10, ceiling(n/50)))}
  
	if(length(breaks) > 1){
		breaks <- sort(unique(breaks))
		k <- length(breaks) - 1
		if(r[1] < breaks[1] | r[2] > breaks[k + 1])
			{stop("all y values must be within the breaks")}
	}
	else{
		k <- breaks
		breaks <- quantile(y, (0:k)/k)
		a <- duplicated(breaks)
		if(any(a)){
			for(j in which(a)){
				if(a[j - 1]){breaks[j] <- NA}
				else{
					h <- abs(y - breaks[j])
					h <- min(h[h > 0])
					h <- min(h/2, (r[2] - r[1])/n/100)
					breaks[j-1] <- breaks[j-1] - h
					breaks[j] <- breaks[j] + h
				}
			}
			breaks <- breaks[!is.na(breaks)]
			k <- length(breaks)
			h <- c(Inf, breaks[2:k] - breaks[1:(k - 1)])
			breaks[h <= 0] <- NA
			breaks <- breaks[!is.na(breaks)]
			k <- length(breaks) - 1
		}
	}
	breaks[1] <- breaks[1] - (breaks[2] - breaks[1])/n
	names(breaks) <- NULL
	list(breaks = breaks, k = k)	
}




pois.newton <- function(start, f, tol = 1e-5, maxit = 200, ...){

	f0 <- f(start, ..., deriv = 2)
	g <- attr(f0, "gradient")
	h <- attr(f0, "hessian")
	conv <- FALSE
	eps <- 1
	alg <- "nr"

	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####
		
		H1 <- try(chol(h), silent = TRUE)
		if(class(H1) != "try-error"){
			if(alg == "gs"){alg <- "nr"; eps <- 1}
			delta <- chol2inv(H1)%*%g
		}
		else{
			if(alg == "nr"){alg <- "gs"; eps <- 1}
			delta <- g
		}

		####

		f1 <- Inf
		while(f1 > f0){
			new.start <- start - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			f1 <- try(f(new.start, ..., deriv = 0), silent = TRUE)
			eps <- eps*0.5
			if(class(f1) == "try-error" || is.na(f1)){f1 <- Inf}
		}

		if(conv | f0 - f1 < tol){break}
		f1 <- f(new.start, ..., deriv = 2)
		g1 <- attr(f1, "gradient")
		h1 <- attr(f1, "hessian")

		start <- new.start; f0 <- f1; g <- g1; h <- h1
		eps <- min(eps*10,1)
	}

	list(estimate = start, n.it = i, minimum = as.numeric(f0), gradient = g, hessian = h)
}




pois.loglik <- function(beta,d,x,w,off, deriv = 0){

	log.lambda <- x%*%cbind(beta)
	lambda <- exp(log.lambda)
	a <- lambda*off
	l <- sum(w*(d*log.lambda - a))
	if(deriv == 0){return(-l)}

	s <- x*c(w*(d - a))
	if(deriv == 1){return(s)}

	s <- colSums(s)
	h <- -t(x*(c(w*a)))%*%x
	out <- -l
	attr(out, "gradient") <- -s
	attr(out, "hessian") <- -h 
	out
}





poisfit <- function(d,x,w,off){

	cn <- colnames(x)
	if(!any(d != 0)){
		sx <- apply(x,2,sd)
		if(int <- (any((const <- (sx == 0))))){
			beta <- rep(NA, ncol(x))
			beta[const] <- -Inf
			names(beta) <- cn
			return(list(beta = beta, vcov = 0))
		}
		else{stop("zero-risk can not be fitted unless an intercept is included")}
	}

	xx <- qr(x)
	sel <- xx$pivot[1:xx$rank]
	x <- x[,sel, drop = FALSE]
	q <- ncol(x)
	V <- list(x = x, off = off)

	mx <- colMeans(x)
	sx <- apply(x,2,sd)
	M <- 10/max(off)
	
	if(int <- (any((const <- (sx == 0))))){
		const <- which(const)
		mx[const] <- 0
		sx[const] <- x[1,const]
		off <- off*M
	}
	else{
		const <- integer(0)
		mx <- rep.int(0, q)
		M <- 1
	}
	vars <- which(sx > 0)
	x <- scale(x, center = mx, scale = sx)

	beta0 <- rep.int(0,q)
	fit <- pois.newton(beta0, pois.loglik, tol = 1e-5, maxit = 10*(1 + q), d = d, x = x, w = w, off = off)

	beta <- fit$estimate
	beta[vars] <- beta[vars]/sx[vars]
	beta[const] <- beta[const] - sum(beta[vars]*mx[vars]) + log(M)
	loglik <- pois.loglik(beta,d,V$x,w,V$off, deriv = 2)
	h <- attr(loglik, "hessian")

	Beta <- rep(NA, length(xx$pivot))
	Beta[sel] <- beta
	names(Beta) <- cn

	list(beta = Beta, vcov = chol2inv(chol(h)))
}






pch.fit <- function(z,y,d,x,w,breaks){ 

	n <- length(y)
	q <- ncol(x)
	if(missing(z)){z <- rep.int(-Inf,n)}
	if(missing(d)){d <- rep.int(1,n)}
	if(missing(w)){w <- rep.int(1,n)}
	Breaks <- makebreaks(y,breaks)

	k <- Breaks$k
	Breaks <- Breaks$breaks
	z[z < min(y)] <- Breaks[1]

	zcut <- cut(z,Breaks, include.lowest = TRUE)
	ycut <- cut(y,Breaks, include.lowest = TRUE)
	tab.ycut <- table(ycut)
	end.z <- match(zcut, levels(zcut))
	end.y <- match(ycut, levels(ycut))
	h <- Breaks[2:(k + 1)] - Breaks[1:k]

	A <- cbind(z,end.z,y,end.y,d,w,x)
	beta <- NULL
	V <- matrix(0,k*q, k*q)
	for(j in 1:k){

		uy <- (A[,"end.y"] == j)
		uz <- (A[,"end.z"] == j)
		vz <- (A[,"end.z"] < j)

		zj <- yj <- rep.int(h[j],n); dj <- rep.int(0,n)
		zj[uz] <- A[uz, "z"] - Breaks[j]; zj[vz]<- 0
		yj[uy] <- A[uy, "y"] - Breaks[j]
		yz_zj <- yj - zj
		dj[uy] <- A[uy, "d"]; wj <- A[,"w"]
		xj <- A[,7:ncol(A),drop = FALSE]

		r <- (yz_zj != 0)
		dj <- dj[r]; yj_zj <- yz_zj[r]
		xj <- xj[r,, drop = FALSE]; wj <- wj[r]

		modj <- poisfit(dj, xj, wj, yj_zj)

		beta <- cbind(beta, modj$beta)
		v <- matrix(0,q,q)
		sel <- which(!is.na(modj$beta))
		ind <- ((j - 1)*q + 1):(j*q)			
		v[sel,sel] <- modj$vcov
		V[ind,ind] <- v

		n <- nrow(A <- A[!uy,, drop = FALSE])
	}

	beta[is.na(beta)] <- 0	
	lambda <- exp(x%*%beta)
	Lambda <- CumSum(t(t(lambda)*h))
	colnames(lambda) <- colnames(Lambda) <- 1:k

	# y and corresponding interval
	
	y <- cbind(y, end.y, sort(y))
	colnames(y) <- c("y","interval","sort.y")
	attr(y, "tab.ycut") <- tab.ycut
	attr(Breaks, "h") <- h; attr(Breaks, "k") <- k
  
	# approxfun for quick prediction of the interval
  
	br <- c(Breaks[1] - 1, Breaks)
	u <- approxfun(br, c(1:k, k + 1, k + 1), rule = 2, method = "constant")  

	list(beta = beta, lambda = lambda, Lambda = Lambda,
		covar = V, breaks = Breaks, y = y, u = u)
}

predF.pch <- function(obj,x,y){

  Breaks <- obj$breaks
  k <- attr(Breaks, "k")
  h <- attr(Breaks, "h")
  beta <- obj$beta
  u <- obj$u
  
  if(missing(x)){
    lambda <- obj$lambda
    Lambda <- obj$Lambda
  }
  else{
    lambda <- exp(x%*%beta)
    Lambda <- CumSum(t(t(lambda)*h))
  }
  
  if(missing(y)){y <- obj$y[,1]; end.y <- obj$y[,2]}
  else{end.y <- u(y); y <- pmax(y, Breaks[1] - 1)}
  
  n <- length(y)
  t <- y - c(0,Breaks)[end.y]
  ind <- cbind(1:n,end.y)
  lambda <- cbind(0,lambda)[ind]
  Lambda <- cbind(0,0,Lambda)[ind] + lambda*t
  SF <- exp(-Lambda)
  PDF <- lambda*SF
  out <- cbind(haz = lambda, Haz = Lambda, Surv = SF, f = PDF)
  rownames(out) <- NULL
  out
}

predQ.pch <- function(obj,x,p){ # p can be a vector (to simulate)

	Breaks <- obj$breaks
	h <- attr(Breaks, "h")
	k <- attr(Breaks, "k")
	beta <- obj$beta

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- exp(x%*%beta)
		Lambda <- CumSum(t(t(lambda)*h))
	}
	
	n <- nrow(lambda)
	if(missing(p)){p <- runif(n); sim <- TRUE; r <- 1}
	else{sim <- FALSE; r <- length(p)}
	pp <- -log(1 - p)

	out <- NULL
	for(j in 1:r){
		tau <- (if(sim) pp else pp[j])
		ind <- .rowSums(Lambda < tau, n, ncol(Lambda)) + 1
		Lambdap <- cbind(0,Lambda)[cbind(1:n, ind)]
		c <- Breaks[ind]
		t <- (tau - Lambdap)/lambda[cbind(1:n, pmin(ind,k))]
		# I use the k-th hazard after the last break
		out <- cbind(out, c + t)
	}
	if(sim){return(c(out))}
	colnames(out) <- paste("p",p, sep = "")
	rownames(out) <- NULL
	out
}

# Note: method = "sample" will only work with "many" breaks.
# On the other hand, it will be more robust to model assumptions,
# and it will correcly reproduce a mass. Moreover, it will work well
# if the data are not compact, while "empty spaces" in the support
# will be filled with method = "quantile".

sim.pch <- function(obj,x,method = c("q","s")){

	if((method <- method[1]) == "q")
		{return(predQ.pch(obj,x))}
	
	Breaks <- obj$breaks
	h <- attr(Breaks, "h")
	k <- attr(Breaks, "k")
	beta <- obj$beta
	y <- obj$y
	tab.ycut <- attr(y, "tab.ycut")
	y <- y[,"sort.y"]	

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- exp(x%*%beta)
		Lambda <- CumSum(t(t(lambda)*h))
	}
	
	n <- nrow(lambda)
	p <- -log(1 - runif(n))
	ind <- .rowSums(Lambda < p, n, ncol(Lambda)) + 1 # can be k + 1
	i1 <- 0
	t <- rep.int(NA,n)
	
	for(j in 1:k){
		uj <- which(ind == j)
		i2 <- i1 + tab.ycut[j]
		t[uj] <- sample(y[(i1 + 1):i2], size = length(uj), replace = TRUE)
		i1 <- i2
	}
	if(any(u.out <- (is.na(t)))){ # where ind = k + 1
		u.out <- which(u.out)
		p.out <- p[u.out]
		lambda.out <- lambda[u.out,k]
		Lambda.out <- Lambda[u.out,k]
		t[u.out] <- Breaks[k + 1] + (p.out - Lambda.out)/lambda.out
	}

	t
}





