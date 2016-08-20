\name{pchreg}
\alias{pchreg}
\title{
Piecewise Constant Hazards Models
}
\description{
This function estimates piecewise exponential models on right-censored, left-truncated data. 
The effect of covariates, and not just the baseline hazard, 
varies across intervals. Moreover, a special handling of zero-risk regions is implemented.
Differently
from the \code{phreg} function available in the \pkg{eha} package,
this function is mainly intended to be used as a nonparametric
maximum likelihood estimator.
}
\usage{
pchreg(formula, breaks, data, weights, splinex = NULL)
}
\arguments{
  \item{formula}{
 an object of class \dQuote{\code{\link{formula}}}: a symbolic description of the regression model. 
The response must be a \kbd{Surv} object as returned by \code{\link{Surv}} (see \sQuote{Details}).
}
  \item{breaks}{
either a numeric vector of two or more unique cut points or a single number (greater than or equal to 1) 
giving the number of intervals into which the time variable is to be cut. If missing, the number of intervals 
is set to \code{max(5, min(50, ceiling(n1/q/5)))}, where \code{n1} is the number of events, and \code{q} is the
number of predictors.
}
  \item{data}{
an optional data frame containing the variables in the model.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process.
}
  \item{splinex}{
either \code{NULL}, or an object created with \code{\link{splinex}} (see \sQuote{Details}). 
}
}
\details{
The left side of \kbd{formula} must be of the form \kbd{Surv(time, event)} if the data are right-censored, and \kbd{Surv(time0, time, event)} if the data are right-censored and left-truncated (\kbd{time0 < time}). Using \kbd{Surv(time)} is also allowed and indicates that the data are neither censored nor truncated.
Note that the response variable (and thus the \code{breaks}) can be negative.

To fit the model, the time interval is first divided in sub-intervals as defined by \code{breaks}.
When the location of \code{breaks} is not specified, the empirical quantiles of \code{time[event == 1]}
are used as cut points. If there is a probability mass, this may result in two or more \code{breaks}
being equal: in this case, an interval that only includes the mass point is created automatically.

A different costant hazard (exponential) model is then fitted in each sub-interval, using Poisson 
regression to model the log-hazard as a linear function of covariates.
Within each interval, the risk of the event may be zero at some covariate values.
For each covariate \code{x}, the algorithm will try to identify a threshold \code{c} 
such that all events (in any given interval) occur when \code{x < c} (\code{x > c}).
A zero risk will be automatically fitted above (below) the threshold, using an offset of \kbd{-100}
on the log-hazard. 

This type of model can be utilized to obtain a nonparametric maximum likelihood estimator 
of a conditional distribution, achieving the flexibility of nonparametric estimators
while keeping the model parametric in practice. Users unfamiliar with this approach
are recommended reading Geman and Hwang (1982) for an overview, and the paper by Ackerberg, Chen and Hahn (2012) describing how this approach can be applied to simplify inference in two-step semiparametric models.

The number of parameters is equal to the number of intervals, multiplied by the number of
covariates. Both quantities are usually supposed to increase with the sample size.
The special function \code{splinex} is a handy tool that facilitates implementing
the linear predictor. When \code{splinex} is not \code{NULL}, each column of 
the original design matrix (as defined by \code{formula}) is automatically replaced with 
the corresponding spline basis. See the documentation of \code{\link{splinex}} for details.
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
\item{mf}{the model frame used.}
\item{x}{the model matrix.}
\item{conv.status}{a code indicating the convergence status. It takes value \kbd{0} if the
algorithm has converged successfully; \kbd{1} if convergence has not been achieved; 
and \kbd{2} if, although convergence has been achieved, more than 1\% of observations
have an associated survival numerically equal to zero, indicating that the solution is not
well-behaved or the model is misspecified.}

The accessor functions \code{summary}, \code{coef}, \code{predict}, \code{nobs}, \code{logLik}, \code{AIC}, \code{BIC} can be used to extract information from the fitted model.
This function is mainly intended for prediction and simulation: see \code{\link{predict.pch}}.
}
\references{
Ackerberg, D., Chen, X., and Hahn, J. (2012). A Practical Asymptotic Variance Estimator for Two-Step Semiparametric Estimators. The Review of Economics and Statistics, 94(2), 481-498.

Friedman, M. (1982). Piecewise Exponential Models for Survival Data with Covariates. The Annals of Statistics, 10(1), pp. 101-113.

Geman, S., and Hwang, C.R. (1982). Nonparametric Maximum Likelihood Estimation by the Method of Sieves.
The Annals of Statistics,10(2), 401-414. 
}
\author{
Paolo Frumento <paolo.frumento@ki.se>
}
\seealso{
\code{\link{predict.pch}}, \code{\link{splinex}}
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
