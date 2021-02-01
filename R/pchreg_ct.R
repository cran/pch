
makebreaks.ct <- function(y,d,x,w, breaks){

  n <- length(y)
  q <- ncol(x)
  r <- range(y)
  n1 <- sum(d)
  y1 <- y[d == 1]
	
  if(missing(breaks)){breaks <- max(5, min(10, ceiling(n1/q/5)))}

  if(length(breaks) > 1){
    breaks <- sort(unique(breaks))
    k <- length(breaks) - 1
    if(r[1] < breaks[1] | r[2] > breaks[k + 1])
      {stop("all y values must be within the breaks")}
    breaks <- breaks[1:which(breaks >= r[2])[1]]
    k <- length(breaks) - 1
  }
  else{
    k <- breaks
    breaks <- wtd.quantile(y1, weights = w, probs = (0:k)/k)

    check.l <- (mean(y < breaks[1]) > 0.03)
    check.r <- (mean(y > breaks[k + 1]) > 0.03)
    if(k > 2 && (check.l | check.r)){ # "many" censored data on the left or right tail
      if(check.l & !check.r){breaks <- c(-Inf, wtd.quantile(y1, weights = w, probs = (0:(k - 1))/(k - 1)))}
      if(!check.l & check.r){breaks <- c(wtd.quantile(y1, weights = w, probs = (0:(k - 1))/(k - 1)), Inf)}
      if(check.l & check.r){breaks <- c(-Inf, wtd.quantile(y1, weights = w, probs = (0:(k - 2))/(k - 2)), Inf)}
    }
    breaks[1] <- r[1]; breaks[k + 1] <- r[2]

    a <- duplicated(round(breaks,8))
    if(any(a)){
      for(j in which(a)){
        if(a[j - 1]){breaks[j] <- NA}
        else{
          h <- abs(y1 - breaks[j])
          h <- min(h[h > 0])
          h <- max(1e-6, min(h/2, (r[2] - r[1])/n/100))
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

  # ensure no obs on the breaks
	
  if(any(y %in% breaks)){
    eps <- min(breaks[2:(k + 1)] - breaks[1:k])/n1
    breaks[1:k] <- breaks[1:k] - eps
    breaks[k + 1] <- breaks[k + 1] + eps
  }
	
  names(breaks) <- NULL
  list(breaks = breaks, k = k)	
}


pch.fit.ct <- function(z,y,d,x,w,breaks){ 

  n <- length(y)
  q <- ncol(x)
  # Unnecessary, but possibly convenient for direct calls in qrcm::test.fit
  if(missing(z)){z <- rep.int(-Inf,n)}
  if(missing(d)){d <- rep.int(1,n)}
  if(missing(w)){w <- rep.int(1,n)}
  Breaks <- suppressWarnings(makebreaks.ct(y,d,x,w,breaks))

  k <- Breaks$k
  Breaks <- Breaks$breaks
  z[z < min(y)] <- Breaks[1]

  zcut <- cut(z,Breaks, include.lowest = TRUE)
  ycut <- cut(y,Breaks, include.lowest = TRUE)
  tab.ycut <- table(ycut)
  end.z <- match(zcut, levels(zcut))
  end.y <- match(ycut, levels(ycut))
  h <- Breaks[2:(k + 1)] - Breaks[1:k]

  conv <- TRUE
  A <- cbind(z,end.z,y,end.y,d,w,x)
  beta <- NULL; rangex <- list()
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
    modj <- pois.fit(dj, xj, wj, yj_zj)

    beta <- cbind(beta, modj$beta)
    rangex[[j]] <- modj$r
    v <- matrix(0,q,q)
    sel <- which(!is.na(modj$beta))
    ind <- ((j - 1)*q + 1):(j*q)
    v[sel,sel] <- modj$vcov
    V[ind,ind] <- v

    n <- nrow(A <- A[!uy,, drop = FALSE])
    conv <- (conv & modj$converged)
  }

  beta[is.na(beta)] <- 0	
  lambda <- cleanlambda(exp(x%*%beta),x,rangex)
  Lambda <- CumSum(t(t(lambda)*h))
  colnames(lambda) <- colnames(Lambda) <- 1:k

  # y and corresponding interval
	
  y <- cbind(y, end.y + 1, sort(y))
  colnames(y) <- c("y","interval","sort.y")
  attr(y, "tab.ycut") <- tab.ycut
  attr(Breaks, "h") <- h; attr(Breaks, "k") <- k
	
  # approxfun for quick prediction of the interval
	
  br <- c(Breaks[1] - 1, Breaks)
  u <- approxfun(br, c(1:k, k + 1, k + 1), rule = 2, method = "constant")  
	
  # check1: convergence
  # check2: there should not be many survival = 1
	
  conv.status <- 0
  if(!conv){conv.status <- 1; warning("the algorithm did not converge", call. = FALSE)}
  else{
    surv <- exp(-cbind(0,Lambda)[cbind(1:nrow(Lambda), end.y)]) # approx survival
    if(mean(surv <= 1e-3) > 0.01){
      conv.status <- 2
      warning(
        "numerically zero survivals fitted: the solution may not be well-behaved", 
        call. = FALSE)
    }
  }

  list(beta = beta, lambda = lambda, Lambda = Lambda,
  covar = V, breaks = Breaks, y = y, u = u, rangex = rangex, conv.status = conv.status)
}





pois.loglik <- function(beta,d,x,w,off, zeror, deriv = 0){

  log.lambda <- tcrossprod(x, t(beta)) - zeror*1e+10
  lambda <- exp(log.lambda)
  a <- lambda*off
  l <- sum(w*(d*log.lambda - a))
  if(deriv == 0){return(-l)}

  s.i <- x*c((d - a))
  s <- .colSums(w*s.i, nrow(s.i), ncol(s.i))
  h <- -crossprod(x*(c(w*a)), x)

  out <- -l
  attr(out, "gradient") <- -s
  attr(out, "hessian") <- -h
  attr(out, "s.i") <- s.i
	
  out
}

pois.fit <- function(d,x,w,off){

  cn <- colnames(x)
  if(!any(d != 0)){
    sx <- myapply(x,sd)
    const <- (if(nrow(x) > 1){sx == 0} else {colnames(x) == "(Intercept)"})
    if(int <- any(const)){
      beta <- rep(NA, ncol(x))
      beta[const] <- -Inf
      names(beta) <- cn
      return(list(beta = beta, vcov = 0, r = cbind(rep(-Inf, ncol(x)), Inf), converged = TRUE))
    }
    else{stop("zero-risk can not be fitted unless an intercept is included")}
  }

  # Handling zero-risk regions

  r <- myapply(x[d == 1,, drop = FALSE], range)
  delta <- r[,2] - r[,1]
  r[,1] <- r[,1] - 0.2*delta
  r[,2] <- r[,2] + 0.2*delta
  zeror <- rep.int(FALSE, length(d))
  for(j in 1:ncol(x)){
    out.l <- (x[,j] < r[j,1])
    out.r <- (x[,j] > r[j,2])
    ml <- mean(out.l); mr <- mean(out.r)
    if(max(ml,mr) < 0.05){r[j,1] <- -Inf; r[j,2] <- Inf}
    else{
      if(ml > mr){outx <- out.l; r[j,2] <- Inf}
      else{outx <- out.r; r[j,1] <- -Inf}
      zeror <- (zeror | outx)
    }
  }
  x <- x*(!zeror)

  # Scaling

  xx <- qr(x)
  sel <- xx$pivot[1:xx$rank]
  x <- x[,sel, drop = FALSE]
  q <- ncol(x); n <- nrow(x)
  X <- list(x = x, off = off)

  mx <- colMeans(x)
  sx <- myapply(x,sd)
  M <- 10/max(off)

  const <- (if(n > 1){sx == 0} else {colnames(x) == "(Intercept)"})
  if(int <- any(const)){
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

  # Fit

  conv <- TRUE
  beta0 <- rep.int(0,q)
  if(int){beta0[const] <- max(-10, log(mean(d[!zeror])))}
  safeit <- 0; fit.ok <- FALSE; count <- 0
  while(!fit.ok){

    fit <- newton(beta0, pois.loglik, tol = 1e-5, maxit = 10*(1 + q), safeit = safeit,
      d = d, x = x, w = w, off = off, zeror = zeror)

    fit.ok <- all(abs(fit$gradient) < sqrt(n)/15 + count)
    count <- count + 1; safeit <- safeit + 2
    if(count == 20){conv <- FALSE; break}
  }

  beta <- fit$estimate
  beta[vars] <- beta[vars]/sx[vars]
  beta[const] <- beta[const] - sum(beta[vars]*mx[vars]) + log(M)
  loglik <- pois.loglik(beta,d,X$x,w,X$off, zeror, deriv = 2)

  Beta <- rep(NA, length(xx$pivot))
  Beta[sel] <- beta
  names(Beta) <- cn

  # covariance matrix

  h <- attr(loglik, "hessian")
  s.i <- attr(loglik, "s.i")
  eps <- 1/(n^2)
  ok <- FALSE
  while(!ok){
    omega <- try(chol2inv(chol(t(s.i*w)%*%s.i + diag(eps, q))), silent = TRUE)
    v <- try(chol2inv(chol(h%*%omega%*%h + diag(eps, q))), silent = TRUE)
    ok <- (!inherits(omega, "try-error") && !inherits(v, "try-error"))
    eps <- eps*2
  }


  list(beta = Beta, vcov = v, r = r, converged = conv)
}






