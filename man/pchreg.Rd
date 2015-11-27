\name{pchreg}
\alias{pchreg}
\title{
Piecewise Constant Hazards Models
}
\description{
This function estimates piecewise exponential models on right-censored, left-truncated data. 
The effect of covariates, and not just the baseline hazard, 
varies across intervals. This differentiates this function
from the \code{phreg} function available in the \pkg{eha} package.
}
\usage{
pchreg(formula, breaks, data, weights)
}
\arguments{
  \item{formula}{
 an object of class \dQuote{\code{\link{formula}}}: a symbolic description of the regression model. 
The response must be a \kbd{Surv} object as returned by \code{\link{Surv}} (see \sQuote{Details}).
}
  \item{breaks}{
either a numeric vector of two or more unique cut points or a single number (greater than or equal to 1) 
giving the number of intervals into which the time variable is to be cut. If missing, the number of intervals 
is set to \code{max(2, min(10, ceiling(n/50)))}. See \dQuote{Notes} for further details.
}
  \item{data}{
an optional data frame containing the variables in the model.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process.
}
}
\details{
The left side of \kbd{formula} must be of the form \kbd{Surv(time, event)} if the data are right-censored, and \kbd{Surv(time0, time, event)} if the data are right-censored and left-truncated (\kbd{time0 < time}). Using \kbd{Surv(time)} is also allowed and indicates that the data are neither censored nor truncated.
Note that the response variable (and thus the \code{breaks}) can be negative.

To fit the model, the time interval is divided in sub-intervals as defined by \code{breaks}.
A different costant hazard (exponential) model is then fitted in each sub-interval, using Poisson regression
to model the log-hazard as a linear function of covariates.
The number of parameters is equal to the number of intervals, multiplied by the number of
covariates. Fitting such a large number of parameters permits achieving the flexibility of nonparametric estimators, while keeping the model parametric in practice.

The fitted model is mainly intended for prediction and simulation: see \code{\link{predict.pch}}.
}
\value{
An object of class \dQuote{\code{pch}}, which is a list with the following items:
\item{call}{the matched call.}
\item{beta}{a matrix of regression coefficients. Rows correspond to covariates, while columns correspond to different time intervals.}
\item{breaks}{the used cut points, with attributes \code{'h'} indicating the length of each interval, and 
\code{'k'} denoting the number of intervals.}
\item{covar}{the estimated asymptotic covariance matrix.}
\item{logLik}{the value of the maximized log-likelihood, with attribute \dQuote{\code{df}} indicating the number of free model parameters.}
\item{lambda}{the fitted hazard values in each interval.}
\item{Lambda}{the fitted cumulative hazard values at the end of each interval.}
\item{mf}{the model frame used}

The accessor functions \code{summary}, \code{coef}, \code{predict}, \code{nobs}, \code{logLik}, \code{AIC}, \code{BIC} can be used to extract information from the fitted model.
}
\references{
Friedman, M. (1982). Piecewise Exponential Models for Survival Data with Covariates. The Annals of Statistics,
10(1), pp. 101-113.
}
\author{
Paolo Frumento <paolo.frumento@ki.se>
}
\seealso{
\code{\link{predict.pch}}
}
\note{
When the location of \code{breaks} is not specified, the empirical quantiles of the time variable
are used as cut points. If there is a probability mass, this may result in two or more \code{breaks}
being equal: in this case, an interval that only includes the mass point is created.

Unless \code{breaks} are user-defined, the first break is set to 
\code{y1 - (y2 - y1)/n}, where \code{y1} and \code{y2}
denote the first and second smallest observed times, respectively. This reflects the assumption
that no events occur at time zero, and ensures that the fitted survival is never equal to \code{1}.
}
\examples{

  # using simulated data
  
  n <- 1000
  x <- runif(n)
  time <- rnorm(n, 1 + x, 1 + x)
  cens <- rnorm(n,2,2)
  y <- pmin(time,cens) # censored variable
  d <- (time <= cens) # indicator of the event

  model <- pchreg(Surv(y,d) ~ x, breaks = 20)

  # see the documentation of predict.pch
}
\keyword{survival}
\keyword{models}
\keyword{regression}
