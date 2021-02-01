#' @importFrom stats model.matrix model.response delete.response model.frame terms prcomp
#' @importFrom stats update model.weights logLik coef nobs vcov .getXlevels runif approxfun sd
#' @importFrom survival Surv is.Surv
#' @importFrom Hmisc wtd.quantile
#' @importFrom splines ns bs

# Note to myself: even using wtd.quantile, there will still be minor differences in how the breaks
  # are computed using a certain dataset and the "equivalent" weighted dataset.

#' @export
pchreg <- function(formula, breaks, data, weights, splinex = NULL){

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	mf$formula <- formula
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")

	if(any((w <- model.weights(mf)) < 0)){stop("negative 'weights'")}
	if(is.null(w)){w <- rep.int(1, nrow(mf)); alarm <- FALSE}
	else{
	  alarm <- (w == 0)
	  sel <- which(!alarm)
	  mf <- mf[sel,]
	  w <- w[sel]
	  w <- w/mean(w)
	}
	if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
	if((n <- nrow(mf)) == 0){stop("zero non-NA cases", call. = FALSE)}
	if(!survival::is.Surv(zyd <- model.response(mf)))
		{stop("the model response must be created with Surv()")}
	type <- attributes(zyd)$type
	
	zyd <- cbind(zyd)
	if(type == "right"){y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]}
	else if(type == "counting"){z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]}
	else if(type == "interval"){y <- convert.Surv(zyd); z <- NULL; d <- zyd[,3]}
	else{stop("only 'right', 'counting', and 'interval2' data are supported")}
	if(!(any(d %in% c(1,3)))){stop("all observation are censored")}
	
	x <- model.matrix(mt,mf); ax <- NULL
	if(!is.null(splinex)){
	  x <- build.splinex(x, splinex$method, splinex$df, splinex$degree, splinex$v)
	  ax <- attributes(x)
	 }

	###

	fit <- (if(type != "interval") pch.fit.ct(z,y,d,x,w,breaks) else pch.fit.ic(z,y,d,x,w,breaks))

	###
  
	attr(mf, "rangex") <- fit$rangex
	attr(mf, "splinex") <- ax
	attr(mf, "u") <- fit$u
	attr(mf, "y") <- fit$y
	attr(mf, "type") <- type
	attr(mf, "n") <- n
	attr(mf, "n.events") <- (if(type != "interval") sum(d) 
		else c(events = sum(d == 1), 'left-censored' = sum(d == 2),
		'right-censored' = sum(d == 0), 'interval-censored' = sum(d == 3))
	)

	beta <- fit$beta
	r <- dim(beta)
	colnames(beta) <- 1:r[2]
	vn <- paste(rep.int(rownames(beta), r[2]), rep(1:r[2], each = r[1]), sep = "_")
	dimnames(fit$covar) <- list(vn, vn)

	# log-likelihood. If type != "interval", the entire log-likelihood was never
	# calculated, and I do it here. If type == "interval", instead, the log-likelihood
	# is already known, and so is the score s.i and the hessian h (to be used in the asymptotics of ctqr).

	if(type != "interval"){
		Hy <- predF.pch(fit, x, y)
		Hz <- (if(type == "counting") predF.pch(fit,x,z)[,"Haz"] else 0)
		l1 <- log(Hy[,"haz"]); l1[l1 == -Inf] <- NA
		logLik <- sum(d*l1*w, na.rm = TRUE) - sum(w*(Hy[,"Haz"] - Hz))
		# note: 'haz' can be 0, I set 0*log(0) = 0.
	}
	else{logLik <- -fit$loglik; attr(mf, "s.i") <- fit$s.i; attr(mf, "h") <- fit$h}
	attr(logLik, "df") <- sum(fit$beta != 0)
	
	attr(cl, "all.vars") <- all.vars(formula)
	fit <- list(call = cl, beta = beta, breaks = fit$breaks, 
		covar = fit$covar, logLik = logLik, 
		lambda = fit$lambda, Lambda = fit$Lambda, 
		mf = mf, x = x, 
		conv.status = fit$conv.status)
	class2 <- (if(type == "interval") "ic" else "ct")
	class(fit) <- c("pch", class2)
	fit
}



#' @export
predict.pch <- function(object, type = c("distr", "quantile", "sim"), 
	newdata, p, sim.method = c("quantile", "sample"), ...){
  
	ic <- (inherits(object, "ic"))
	if(is.na(match(type <- type[1], c("distr", "d", "quantile", "q", "sim", "s"))))
		{stop("invalid 'type'")}
	type <- strsplit(type, split = "")[[1]][1]
	if(type == "s"){
		if(is.na(match(method <- sim.method[1], c("quantile", "q", "sample", "s"))))
		{stop("invalid 'method'")}
		method <- strsplit(method, split = "")[[1]][1]
		if(method == "s" & ic)
			{stop("simulation method 'sample' is not applicable with interval censoring")}
	}
	if(type == "q"){
	  if(missing(p)){stop("please indicate 'p'")}
	  if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}	
	}
	
	mf <- object$mf
	mt <- terms(mf)
	xlev <- .getXlevels(mt, mf)
	fittype <- attr(mf, "type")
	splinex <- attr(mf, "splinex")
	
	obj <- list(beta = object$beta, breaks = object$breaks, y = attr(mf, "y"), 
		lambda = object$lambda, Lambda = object$Lambda, u = attr(mf, "u"), 
		rangex = attr(mf, "rangex"))

	if((nodata <- missing(newdata))){
		if(ic && type == "d"){stop("with type = 'distr' and interval-censored data, 'newdata' is always required \n and must include a variable named 'time' at which to compute predictions")}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	}
	else{
		if(type == "d"){
			yn <- as.character(if(fittype == "interval") "time" else if(fittype == "counting") mt[[2]][[3]] else mt[[2]][[2]])
			if(is.na(ind <- match(yn, colnames(newdata)))){
			  if(!ic){stop("for 'type = distr', 'newdata' must contain new values of the y-variable")}
			  else{stop("with type = 'distr' and interval-censored data, 'newdata' must contain a variable named 'time' at which to compute predictions")}
			}
			if(fittype == "right" && length(mt[[2]]) == 3)
			  {newdata[,as.character(mt[[2]][[3]])] <- 1}
			if(fittype == "counting"){
				newdata[,as.character(mt[[2]][[4]])] <- 1
				newdata[,as.character(mt[[2]][[2]])] <- -Inf
				newdata[,yn] <- pmax(newdata[,yn], obj$breaks[1] - 1)
			}
			if(fittype == "interval"){mt <- update(mt, time ~ .)}
		}
		else{mt <- delete.response(mt)}
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
		{stop("'newdata' must contain all x-variables")}
		mf <- model.frame(mt, data = newdata, xlev = xlev)
		if(nrow(mf) == 0){stop("zero non-missing values in the supplied newdata")}

		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
    
		x <- model.matrix(mt, mf, contrasts.arg = attr(object$x, "contrasts"))
		if(!is.null(splinex)){
		  x <- predict.splinex(x,splinex)
		}
	}

	n <- length(miss) + length(nomiss)
	if(type == "d"){
		if(!nodata){y <- cbind(model.response(mf))[,1 + (fittype == "counting")]}
		out <- matrix(, n, 4)
		pred <- (if(nodata) predF.pch(obj) else predF.pch(obj,x,y))
		out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
  	else if(type == "q"){
		out <- matrix(, n, length(p))
		pred <- (if(nodata) predQ.pch(obj, p = p) else predQ.pch(obj,x,p))
        	out[nomiss,] <- pred
		dimnames(out) <- list(1:n, colnames(pred))
		return(as.data.frame(out))
	}
	else{
		tsim <- rep.int(NA,n)
		tsim[nomiss] <- (if(nodata) sim.pch(obj, method = method) else sim.pch(obj,x,method))
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
  out <- list(call = object$call, class = class(object)[2],
    n = attr(object$mf, "n"), n.events = attr(object$mf, "n.events"), 
		n.free.par = attr(object$logLik, "df"), logLik = object$logLik, conv.status = object$conv.status)
  class(out) <- "summary.pch"
  out
}

#' @export
print.summary.pch <- function(x, ...){
  cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("n. of obs: ", paste(deparse(round(x$n)), sep = " ", collapse = " "), "\n", sep = "")
  
  if(x$class == "ct"){
    cat("n. of events: ", paste(deparse(round(x$n.events)), sep = " ", collapse = " "), "\n", sep = "")
  }
  else{
    cat("non-censored: ", paste(deparse(round(x$n.events[[1]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("left-censored: ", paste(deparse(round(x$n.events[[2]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("right-censored: ", paste(deparse(round(x$n.events[[3]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("interval-censored: ", paste(deparse(round(x$n.events[[4]])), sep = " ", collapse = " "), "\n", sep = "")
  }
  cat("\n")
  cat("n. of free parameters: ", paste(deparse(round(x$n.free.par)), sep = " ", collapse = " "), "\n", sep = "")
  cat("Log-likelihood: ", paste(deparse(round(as.numeric(x$logLik),1)), sep = " ", collapse = " "), "\n", sep = "")
  cat("convergence status: ", paste(deparse(round(x$conv.status)), sep = " ", collapse = " "), "\n", sep = "")
  cat("\n")
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



predF.pch <- function(obj,x,y){

  Breaks <- obj$breaks
  k <- attr(Breaks, "k")
  h <- attr(Breaks, "h")
  beta <- obj$beta
  u <- obj$u
  rangex <- obj$rangex

  if(missing(x)){
    lambda <- obj$lambda
    Lambda <- obj$Lambda
  }
  else{
    lambda <- cleanlambda(exp(x%*%beta), x, rangex)
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
	rangex <- obj$rangex

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- cleanlambda(exp(x%*%beta), x, rangex)
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
	rangex <- obj$rangex
	y <- obj$y
	tab.ycut <- attr(y, "tab.ycut")
	y <- y[,"sort.y"]	

	if(missing(x)){
		lambda <- obj$lambda
		Lambda <- obj$Lambda
	}
	else{
		lambda <- cleanlambda(exp(x%*%beta), x, rangex)
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

newton <- function(start, f, tol = 1e-5, maxit = 200, safeit = 0, ...){

	start <- c(start)
	f0 <- f(start, ..., deriv = 2)
	g <- attr(f0, "gradient")
	h <- attr(f0, "hessian")
	conv <- FALSE
	eps <- 1
	alg <- "gs"
		
	for(i in 1:maxit){

		if(conv | max(abs(g)) < tol){break}

		####
		
	  if(i > safeit){
		  H1 <- try(chol(h), silent = TRUE)
		  if(!inherits(H1, "try-error")){
		  	if(alg == "gs"){alg <- "nr"; eps <- 1}
		  	delta <- chol2inv(H1)%*%g
		  }
		  else{
		  	if(alg == "nr"){alg <- "gs"; eps <- 1}
		  	delta <- g
		  }
	  }
		else{delta <- g}
	  
		####

		f1 <- Inf
		while(f1 > f0){
			new.start <- start - delta*eps
			if(max(abs(delta*eps)) < tol){conv <- TRUE; break}
			f1 <- try(f(new.start, ..., deriv = 0), silent = TRUE)
			eps <- eps*0.5
			if(inherits(f1, "try-error") || is.na(f1)){f1 <- Inf}
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



