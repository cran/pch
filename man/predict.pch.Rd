\name{predict.pch}
\alias{predict.pch}
\title{
Prediction from Fitted Piecewise Constant Hazard Models
}
\description{
This function returns predictions for an object of class \dQuote{\code{pch}}, usually the result of a call
to \code{\link{pchreg}}.
}
\usage{
\method{predict}{pch}(object, type = c("distr", "quantile", "sim"), 
   newdata, p, sim.method = c("quantile", "sample"), ...)
}
\arguments{
  \item{object}{
  a \dQuote{\code{pch}} object.
}
  \item{type}{
a character string (just the first letter can be used) indicating the type of prediction. See \sQuote{Details}.
}
  \item{newdata}{
optional data frame in which to look for variables with which to predict. It must include all the covariates that enter the model and, if \code{type = 'distr'}, also the time variable (see \sQuote{Details} for additional information of interval-censored data). 
If \code{newdata} is omitted, the original data will be used.
}
  \item{p}{
vector of quantiles, to be specified if \code{type = "quantile"}.
}
  \item{sim.method}{
a character string (just the first letter can be used) indicating the simulation method if \code{type = "sim"}. Only \code{sim.method = 'quantile'} is valid with interval-censored data.
}
  \item{\dots}{
for future methods.
}
}
\details{
If \code{type = "distr"} (the default), this function returns a data frame with columns \code{(haz, Haz, Surv, f)}
containing the fitted values of the hazard function, the cumulative hazard, the survival function, and 
the probability density function, respectively. 

If \code{type = "quantile"}, a data frame with the fitted quantiles (corresponding to the supplied 
values of \code{p}) is returned.

If \code{type = "sim"}, new data are simulated from the fitted model. Two methods are available:
with \code{sim.method = "quantile"}, data are simulated by applying the estimated quantile function
to a vector of random uniform numbers; if \code{sim.method = "sample"}, the quantile function is only used to identify the time interval, and the data are resampled from the observed values in the interval.
The second method only works properly if there is a large number of breaks. However, it is less sensitive to
model misspecification and facilitates sampling from distributions with a probability mass or non compact support. This method is not applicable to interval-censored data.

Predictions are computed at \code{newdata}, if supplied. Note that \code{newdata} 
must include all the variables that are needed for the prediction, and that if \code{type = "distr"},
new values of the response variable are also required. If the data are interval-censored between \code{time1}
and \code{time2}, these will \emph{not} be used as time-to-events and \code{newdata} must include
a variable \code{'time'} at which to compute predictions.

}
\value{
If \code{type = "distr"}, a 4-columns data frame with columns \code{(haz, Haz, Surv, f)}.
If \code{type = "quantile"}, a named data frame with a column for each value of \code{p}.
If \code{type = "sim"}, a vector of simulated data.

The presence of missing values in the response or the covariates will always cause the prediction to be \code{NA}.
}
\author{
Paolo Frumento <paolo.frumento@unipi.it>
}
\note{
If the data are right-censored, some high quantiles may not be estimated: beyond the last observable quantile,
all types of predictions (including \code{type = "sim"} with \code{sim.method = "sample"}) are
computed assuming that the hazard remains constant after the last interval.
}
\seealso{
\code{\link{pchreg}}
}
\examples{

  # using simulated data
  
  ##### EXAMPLE 1 - Continuous distribution ############################
  
  n <- 1000
  x <- runif(n)
  time <- rnorm(n, 1 + x, 1 + x) # time-to-event
  cens <- rnorm(n,2,2) # censoring variable
  y <- pmin(time,cens) # observed variable
  d <- (time <= cens) # indicator of the event
  model <- pchreg(Surv(y,d) ~ x, breaks = 20)

  # predicting hazard, cumulative hazard, survival, density

  pred <- predict(model, type = "distr")
  plot(pred$Surv, 1 - pnorm(y, 1 + x, 1 + x)); abline(0,1) 
  # true vs fitted survival
  
  
  # predicting quartiles

  predQ <- predict(model, type = "quantile", p = c(0.25,0.5,0.75))
  plot(x,time)
  points(x, qnorm(0.5, 1 + x, 1 + x), col = "red") # true median
  points(x, predQ$p0.5, col = "green")             # fitted median
  
  
  # simulating new data
  
  tsim1 <- predict(model, type = "sim", sim.method = "quantile")
  tsim2 <- predict(model, type = "sim", sim.method = "sample")

  qt <- quantile(time, (1:9)/10)  # deciles of t
  q1 <- quantile(tsim1, (1:9)/10) # deciles of tsim1
  q2 <- quantile(tsim2, (1:9)/10) # deciles of tsim2

  par(mfrow = c(1,2))
  plot(qt,q1, main = "sim.method = 'quantile'"); abline(0,1)
  plot(qt,q2, main = "sim.method = 'sample'"); abline(0,1)

  # prediction with newdata
  
  predict(model, type = "distr", newdata = data.frame(y = 0, x = 0.5)) # need y!
  predict(model, type = "quantile", p = 0.5, newdata = data.frame(x = 0.5))
  predict(model, type = "sim", sim.method = "sample", newdata = data.frame(x = c(0,1)))

  ##### EXAMPLE 2 - non-compact support ############################
  # to simulate, sim.method = "sample" is recommended ##############
  
  n <- 1000
  t <- c(rnorm(n,-5), rnorm(n,5)) 
  model <- pchreg(Surv(t) ~ 1, breaks = 30)
  
  tsim1 <- predict(model, type = "sim", sim.method = "quantile")
  tsim2 <- predict(model, type = "sim", sim.method = "sample")
  
  par(mfrow = c(1,3))
  hist(t, main = "true distribution")
  hist(tsim1, main = "sim.method = 'quantile'") # the empty spaces are 'filled'
  hist(tsim2, main = "sim.method = 'sample'")   # perfect!
}
\keyword{regression}
