
myapply <- function(x,fun){
  out <- NULL
  for(j in 1:ncol(x)){out <- rbind(out, fun(x[,j]))}
  if(ncol(out) == 1){out <- c(out)}
  out
}


cleanlambda <- function(lambda,x,rangex){
  for(h in 1:ncol(lambda)){
    r <- rangex[[h]]
    for(j in 1:ncol(x)){
      e <- which(x[,j] < r[j,1] | x[,j] > r[j,2])
      lambda[e,h] <- 0
    }
  }
  pmin(lambda, 1e+6)
}


# splinex functions

#' @export
splinex <- function(method = c("ns", "bs"), df = 2, degree = 2, v = 0.98, ...){
  if(is.na(match((method <- method[1]), c("ns", "bs")))){stop("invalid 'method'")}
  if((df <- round(df)) < 1){stop("invalid 'df'")}
  if((degree <- round(degree)) < 1){stop("invalid 'degree'")}
  if(method == "bs" && df < degree){df <- degree; warning(paste("'df' was too small; have used", df))}
  if(v <= 0 | v > 1){stop("0 < v <= 1 is required")}
  list(method = method, df = df, degree = degree, v = v)
}

build.splinex <- function(x, method, df, degree, v){

  m <- (if(method == "ns") function(x, df, deg) splines::ns(x, df = df) 
        else function(x, df, deg) splines::bs(x, df = df, degree = deg))
  x <- cbind("(Intercept)" = 1, as.matrix(x))
  qx <- qr(x)
  sel1 <- qx$pivot[1:qx$rank]
  x <- x[,sel1, drop = FALSE]
  cnx <- colnames(x)
  
  if(ncol(x) == 1){
    attr(x, "onlyint") <- TRUE
    return(x)
  }
  if(df == 1){
    X <- x
    df <- rep.int(1, ncol(X))
    sel2 <- 1:ncol(X)
    knots <- bknots <- NULL
  }
  else{
    u <- myapply(x, function(a){length(unique(a))})
    X <- cbind("(Intercept)" = rep.int(1, nrow(x)))
    dfX <- 1; cnX <- "(Intercept)"
    knots <- bknots <- list()
    for(j in 2:ncol(x)){
      dfj <- dfX[j] <- min(df, u[j] - 1)
      xx <- (if(dfj > 1) m(x[,j], df = dfj, deg = degree) else x[,j])
      X <- cbind(X, xx)
      cnX <- c(cnX, paste(cnx[j],1:dfj, sep = "."))
      knots[[j]] <- attr(xx, "knots")
      bknots[[j]] <- attr(xx, "Boundary.knots")
    }
    colnames(X) <- cnX
    df <- dfX
    qX <- qr(X)
    sel2 <- qX$pivot[1:qX$rank]
    X <- X[,sel2, drop = FALSE]
  }
  if(v == 1){rot <- diag(ncol(X) - 1); center <- scale <- FALSE; sel3 <- 1:(ncol(X) - 1)}
  else{
    XX <- prcomp(X[,-1, drop = FALSE], center = TRUE, scale. = TRUE)
    rot <- XX$rotation; center <- XX$center; scale <- XX$scale
    vx <- cumsum(XX$sdev^2); vx <- vx/vx[length(vx)]
    sel3 <- 1:(which(vx >= v)[1])
    X <- cbind(X[,1,drop = FALSE],XX$x[,sel3])
    colnames(X) <- c("(Intercept)", paste("pc", 1:(ncol(X) - 1), sep = ""))
  }
  
  attr(X, "onlyint") <- FALSE
  attr(X, "method") <- method[1]
  attr(X, "sel1") <- sel1
  attr(X, "df") <- df
  attr(X, "degree") <- degree
  attr(X, "knots") <- knots
  attr(X, "bknots") <- bknots
  attr(X, "sel2") <- sel2
  attr(X, "center") <- center
  attr(X, "scale") <- scale
  attr(X, "rot") <- rot
  attr(X, "sel3") <- sel3
  X
}

predict.splinex <- function(x, a){
  
  if(a$onlyint){return(x)}

  m <- (if(a$method == "ns") function(x, df, deg, kn, bkn) splines::ns(x, df = df, knots = kn, Boundary.knots = bkn) 
        else function(x, df, deg, kn, bkn) splines::bs(x, df = df, degree = deg, knots = kn, Boundary.knots = bkn))
  
  x <- cbind(1,x)
  x <- x[,a$sel1, drop = FALSE]
  X <- 1
  for(j in 2:ncol(x)){
    xx <- (if(a$df[j] == 1) x[,j] else m(x[,j], df = a$df[j], deg = a$degree, 
      kn = a$knots[[j]], bkn = a$bknots[[j]]))
    X <- cbind(X, xx)
  }
  X <- X[,a$sel2, drop = FALSE]
  X <- X[,-1,drop = FALSE]
  X <- scale(X, center = a$center, scale = a$scale)
  X <- X%*%a$rot
  X <- X[,a$sel3, drop = FALSE]
  cbind(1,X)
}
