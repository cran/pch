makebreaks.ic <- function(y,d,x,w, breaks){

  n <- nrow(y)
  q <- ncol(x)
  yy <- c(y[,1], y[,2])
  r <- range(yy[abs(yy) != Inf])
  n1 <- sum(o <- (d %in% c(1,3)))
  yc <- (y[o,1] + y[o,2])/2 # centers of the intervals (excl left/right censored data)
	
  if(missing(breaks)){breaks <- max(5, min(10, ceiling(n1/q/5)))}

  if(length(breaks) > 1){
    breaks <- sort(unique(breaks))
    k <- length(breaks) - 1
    if(r[1] < breaks[1] | r[2] > breaks[k + 1])
     {stop("all finite y values must be within the breaks")}
    breaks <- breaks[1:which(breaks >= r[2])[1]]
    k <- length(breaks) - 1
  }
  else{
    k <- breaks
    breaks <- wtd.quantile(yc, weights = w, probs = (0:k)/k)
    check.l <- (mean(yy < breaks[1]) > 0.03)
    check.r <- (mean(yy > breaks[k + 1]) > 0.03)
    if(k > 2 && (check.l | check.r)){ # "many" censored data on the left or right tail
      if(check.l & !check.r){breaks <- c(-Inf, wtd.quantile(yc, weights = w, probs = (0:(k - 1))/(k - 1)))}
      if(!check.l & check.r){breaks <- c(wtd.quantile(yc, weights = w, probs = (0:(k - 1))/(k - 1)), Inf)}
      if(check.l & check.r){breaks <- c(-Inf, wtd.quantile(yc, weights = w, probs = (0:(k - 2))/(k - 2)), Inf)}
    }
    breaks[1] <- r[1]; breaks[k + 1] <- r[2]

    a <- duplicated(round(breaks,8))
    if(any(a)){
      for(j in which(a)){
        if(a[j - 1]){breaks[j] <- NA}
        else{
          h <- abs(yc - breaks[j])
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
	
  if(any(yy %in% breaks)){
    eps <- min(breaks[2:(k + 1)] - breaks[1:k])/n1
    breaks[1:k] <- breaks[1:k] - eps
    breaks[k + 1] <- breaks[k + 1] + eps
  }
	
  names(breaks) <- NULL
  list(breaks = breaks, k = k)	
}


# Take a Surv(type = "interval2") object (which has a quite useless format) and returns (y1, y2)
# Where y1 = y2, y1 will be replaced by y1 - eps, and y2 will be replaced by y2 + eps.

convert.Surv <- function(y){

  code <- y[,3]
    
  n <- nrow(y)
  y1 <- rep.int(-Inf,n)
  y2 <- rep.int(Inf,n)
  
  w1 <- (code != 2)
  y1[w1] <- y[w1,1]

  w2 <- (code != 0)
  y2[w2] <- y[w2,2]

  w2bis <- (code == 2)
  y2[w2bis] <- y[w2bis,1]

  w2ter <- (code == 1)
  y2[w2ter] <- y[w2ter,1]

  yy <- c(y1,y2); yy <- yy[is.finite(yy)]
  eps <- (max(yy) - min(yy))/1e+5
  w <- (code == 1)
  y1[w] <- y1[w] - eps
  y2[w] <- y2[w] + eps
  
  y <- cbind(y1, y2)
  attr(y, "code") <- code
  attr(y, "eps") <- eps
  y
}


# pch.fit.ic and icpch.fit could be merged into a single function.
# I keep them separated for consistency with pch.fit.ct and pois.fit.


pch.fit.ic <- function(z,y,d,x,w,breaks){ 

  # Unnecessary, but possibly convenient for direct calls in qrcm::test.fit
  if(missing(z)){z <- rep.int(-Inf,n)}
  if(missing(d)){d <- rep.int(1,n)}
  if(missing(w)){w <- rep.int(1,n)}


  Breaks <- suppressWarnings(if(all(attr(y, "code") %in% 0:1)) makebreaks.ct(y[,1],is.finite(y[,2]),x,w,breaks)
    else makebreaks.ic(y,d,x,w,breaks))
  n <- nrow(y); k <- Breaks$k
  Breaks <- Breaks$breaks
  h <- Breaks[2:(k + 1)] - Breaks[1:k]
  attr(Breaks, "h") <- h; attr(Breaks, "k") <- k

  mod <- icpch.fit(y,d,x,w,Breaks)

  beta <- mod$beta
  beta[is.na(beta)] <- 0	
  lambda <- cleanlambda(exp(x%*%beta),x, mod$rangex)
  Lambda <- CumSum(t(t(lambda)*h))
  colnames(lambda) <- colnames(Lambda) <- 1:k

	
  # approxfun for quick prediction of the interval
	
  br <- c(Breaks[1] - 1, Breaks)
  u <- approxfun(br, c(1:k, k + 1, k + 1), rule = 2, method = "constant")
	
  conv.status <- 0
  if(!mod$conv){conv.status <- 1; warning("the algorithm did not converge", call. = FALSE)}
  list(beta = beta, lambda = lambda, Lambda = Lambda, loglik = mod$loglik, s.i = mod$s.i, h = mod$h,
     covar = mod$vcov, breaks = Breaks, y = y, u = u, rangex = mod$rangex, conv.status = conv.status)
}



icpch.fit <- function(y,d,x,w,breaks){


  # Handling zero-risk regions

  k <- length(breaks) - 1
  n <- nrow(y)
  zeror <- NULL
  rangex <- list()
  for(j in 1:k){
    open <- (y[,1] >= breaks[j] | y[,1] <= breaks[j+1])
    close <- (y[,2] >= breaks[j] | y[,2] <= breaks[j+1])
    r <- myapply(x[(open | close),, drop = FALSE], range)
    delta <- r[,2] - r[,1]
    r[,1] <- r[,1] - 0.2*delta
    r[,2] <- r[,2] + 0.2*delta
    zeror.j <- rep.int(FALSE, n)
    for(h in 1:ncol(x)){
      out.l <- (x[,h] < r[h,1])
      out.r <- (x[,h] > r[h,2])
      ml <- mean(out.l); mr <- mean(out.r)
      if(max(ml,mr) < 0.05){r[h,1] <- -Inf; r[h,2] <- Inf}
      else{
        if(ml > mr){outx <- out.l; r[h,2] <- Inf}
        else{outx <- out.r; r[h,1] <- -Inf}
        zeror.j <- (zeror.j | outx)
      }
    }
    zeror <- cbind(zeror, zeror.j)
    rangex[[j]] <- r
  }

  # Scaling x

  cn <- colnames(x)
  qq <- ncol(x)
  xx <- qr(x)
  sel <- xx$pivot[1:xx$rank]
  x <- x0 <- x[,sel, drop = FALSE]
  q <- ncol(x)

  mx <- colMeans(x)
  sx <- myapply(x,sd)

  const <- (sx == 0)
  if(int <- any(const)){
    const <- which(const)
    mx[const] <- 0
    sx[const] <- x[1,const]
  }
  else{
    const <- integer(0)
    mx <- rep.int(0, q)
  }
  vars <- which(sx > 0)
  x <- scale(x, center = mx, scale = sx)


  # Fit
	
  conv <- TRUE
  U1 <- break.y(y[,1], breaks)
  U2 <- break.y(y[,2], breaks)
  beta0 <- matrix(0,q,k)
  if(int){beta0[const,] <- max(-10, -log(mean(breaks - breaks[1])))}
  safeit <- 0; fit.ok <- FALSE; count <- 0
  while(!fit.ok){

    fit <- newton(beta0, icpch.loglik, tol = 1e-5, maxit = 10*(1 + q*k), safeit = safeit,
      y = y, x = x, w = w, U1 = U1, U2 = U2, zeror = zeror)

    fit.ok <- all(abs(fit$gradient) < sqrt(n)/15 + count)
    count <- count + 1; safeit <- safeit + 2
    if(count == 20){conv <- FALSE; break}
  }

  beta <- matrix(fit$estimate, q,k)
  beta[vars,] <- beta[vars,]/sx[vars]
  beta[const,] <- beta[const,] - colSums(beta[vars,, drop = FALSE]*mx[vars])
  loglik <- icpch.loglik(beta, y, x0,w,U1,U2, zeror, deriv = 2, final = TRUE)

  Beta <- matrix(NA, qq, k)
  Beta[sel,] <- beta
  rownames(Beta) <- cn

  # covariance matrix

  h <- attr(loglik, "hessian")
  s.i <- attr(loglik, "s.i")
  eps <- 1/(n^2)
  ok <- FALSE
  while(!ok){
    omega <- try(chol2inv(chol(t(s.i*w)%*%s.i + diag(eps, q*k))), silent = TRUE)
    v <- try(chol2inv(chol(h%*%omega%*%h + diag(eps, q*k))), silent = TRUE)
    ok <- (!inherits(omega, "try-error") && !inherits(v, "try-error"))
    eps <- eps*5
  }

  V <- matrix(0, qq*k, qq*k)
  sel2 <- rep.int(sel, k)
  sel2 <- sel2 + rep(0:(k - 1), each = length(sel))*q
  V[sel2,sel2] <- v

  # Note: s.i and h are used by ctqr to compute asymptotics.
  list(beta = Beta, vcov = V, rangex = rangex, converged = conv, loglik = as.numeric(loglik), s.i = s.i, h = h)
}



# For each y, compute the contribution, in terms of time, to each interval
break.y <- function(y,breaks){
  k <- length(breaks) - 1
  U <- NULL
  for(j in 1:k){
    u <- pmin(y - breaks[j], breaks[j + 1] - breaks[j])
    U <- cbind(U, pmax(u,0))
  }
  U
}

# If final == TRUE, compute (numerical) PDF when time1 = time2.
# This is only used to compute the log-likelihood properly
# (if I use S1 - S2, this is not the PDF, because it is not a numerical derivative)

icpch.loglik <- function(beta, y, x,w,U1,U2, zeror, deriv = 0, final = FALSE){

  k <- ncol(U1)
  q <- ncol(x)
  n <- nrow(x)
  beta <- matrix(beta, q,k)

  log.lambda <- tcrossprod(x, t(beta))*(!zeror) - zeror*1e+10
  lambda <- exp(log.lambda)
  lambda <- pmin(lambda, 1e+10) # Sometimes I get lambda = Inf.

  H1 <- .rowSums(lambda*U1, n, k)
  H2 <- .rowSums(lambda*U2, n, k)
  S1 <- exp(-H1); S1[y[,1] == -Inf] <- 1
  S2 <- exp(-H2); S2[y[,2] == Inf] <- 0
  deltaS <- pmax(S1 - S2, 1e-12)
  
  dS <- deltaS
  if(final){ 
    events <- which(attr(y, "code") == 1)
    dS[events] <- dS[events]/2/attr(y, "eps")
  }

  l <- sum(w*log(dS))
  if(deriv == 0){return(-l)}

  # Building blocks

  xw <- x*w
  deltaSU <- S1*U1 - S2*U2
  deltaU <- U2 - U1
  lambda.deltaU <- lambda*deltaU
  lambda.deltaSU <- lambda*deltaSU

  A1 <- lambda.deltaSU/deltaS
  A2 <- S1*S2/deltaS^2
  A3 <- A1 + A2*lambda.deltaU^2
  xA2 <- x*A2

  # score

  s.i <- NULL
  for(j in 1:k){s.i <- cbind(s.i, -x*A1[,j])}
  s <- .colSums(w*s.i, n, q*k)
  
  # hessian

  h <- matrix(NA, q*k, q*k)
  for(j1 in 1:k){
    ind1 <- (j1*q - q + 1):(j1*q)
    for(j2 in j1:k){
      ind2 <- (j2*q - q + 1):(j2*q)
      if(j1 == j2){h[ind1,ind2] <- -crossprod(xw, x*A3[,j1])}
      else{h[ind1,ind2] <- h[ind2,ind1] <- -crossprod(xw, xA2*lambda.deltaU[,j1]*lambda.deltaU[,j2])}
    }
  }

  out <- -l
  attr(out, "gradient") <- -s
  attr(out, "hessian") <- -h
  attr(out, "s.i") <- s.i
	
  out

}















