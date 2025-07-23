\name{pchreg}
\alias{pchreg}
\title{
Piecewise Constant Hazard Models
}
\description{
This function estimates piecewise exponential models on right-censored, left-truncated, or interval-censored data.
The function is mainly intended for prediction and, unlike the \code{phreg} function available in the \pkg{eha} package, 
it allows the effect of covariates, and not just the baseline hazard, to depend on time.
}
\usage{
pchreg(formula, breaks, data, weights, splinex = NULL)
}
\arguments{
  \item{formula}{
 an object of class \dQuote{\code{\link{formula}}}: a symbolic description of the regression model. 
The response must be a \kbd{Surv} object as returned by \code{\link[survival]{Surv}} (see \sQuote{Details}).
}
  \item{breaks}{
either a numeric vector of two or more unique cut points, or a single number 
giving the number of intervals. If missing, the number and position of the \kbd{breaks}
are determined automatically.
}
  \item{data}{
an optional data frame containing the variables in the model.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process. 
The weights will be normalized to sum to the sample size. 
This implies that, for example, using double weights will not halve the standard errors.
}
  \item{splinex}{
either \code{NULL}, or an object created with \code{\link{splinex}} (see \sQuote{Details}). 
}
}
\details{
The left side of the \kbd{formula} must be specified as \kbd{Surv(time, event)}, for right-censored data; 
\kbd{Surv(time0, time, event)}, for right-censored and left-truncated data 
(\kbd{time0 < time}, \code{time0} can be \kbd{-Inf});
and \kbd{Surv(time1, time2, type = "interval2")} for interval-censored data (use \kbd{time1 = time2} 
for exact observations, \kbd{time1 = -Inf} or \kbd{NA} for left-censored, and \kbd{time2 = Inf} or \kbd{NA} 
for right-censored). Using \kbd{Surv(time)} is also allowed and indicates that the data are neither censored 
nor truncated. Note that the response variable (and thus the \code{breaks}) can be negative.

To fit the model, the time interval is first divided in sub-intervals as defined by \code{breaks}.
When the location of \code{breaks} is not specified, the empirical quantiles are used as cut points.
A different costant hazard (exponential) model is then fitted in each sub-interval, modeling 
the log-hazard as a linear function of covariates.
The special function \code{\link{splinex}} can be used to build flexible models.

This type of model can be utilized to obtain a nonparametric maximum likelihood estimator 
of a conditional distribution, achieving the flexibility of nonparametric estimators
while keeping the model parametric in practice. Users unfamiliar with this approach
are recommended to read Geman and Hwang (1982) for an overview, and the paper by Ackerberg, Chen and Hahn (2012) 
describing how this approach can be applied to simplify inference in two-step semiparametric models.


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
have an associated survival numerically equal to zero, indicating that the solution may not be
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
\note{
NOTE1. Right-censoring is a special case of interval censoring, in which exact events are identified by
\code{time2 = time1}, while censored observations have \code{time2 = Inf}.
Note, however, that \code{pchreg} will \emph{not} use the same routines for right-censored
and interval-censored data, implying that
\code{pchreg(Surv(time1, time2, type = "interval2") ~ x)}
may not be identical to \code{pchreg(Surv(time = time1, event = (time2 < Inf)) ~ x)}.
The latter is usually faster and slightly more accurate.

NOTE2. Within each interval, the risk of the event may be zero at some covariate values.
For each covariate \code{x}, the algorithm will try to identify a threshold \code{c} 
such that all events (in any given interval) occur when \code{x < c} (\code{x > c}).
A zero risk will be automatically fitted above (below) the threshold, using an offset of \kbd{-100}
on the log-hazard. 
}
\author{
Paolo Frumento <paolo.frumento@unipi.it>
}
\seealso{
\code{\link{predict.pch}}, \code{\link{splinex}}
}
\examples{

  # Simulate right-censored data
  
  n <- 1000
  x <- runif(n) # a covariate
  time <- rexp(n, exp(1 + x)) # time-to-event
  cens <- runif(n,0,2) # censoring event
  y <- pmin(time,cens) # observed variable
  d <- (time <= cens) # indicator of the event
  model <- pchreg(Surv(y,d) ~ x, breaks = 10)




  # Simulate right-censored, left-truncated data
  
  n <- 1000
  x <- runif(n) # a covariate
  time0 <- rexp(n, 10) # time at enrollment
  time <- rexp(n, exp(1 + x)) # time-to-event
  cens <- runif(n,0,2) # censoring event
  
  # y,d,x are only observed if (y > time0)
  y <- pmin(time,cens)
  d <- (time <= cens)
  u <- (y > time0)
  y <- y[u]
  d <- d[u]
  x <- x[u]
  z <- time0[u]
  model <- pchreg(Surv(z,y,d) ~ x, breaks = 10)




  # Simulate interval-censored data
  
  n <- 1000
  x <- runif(n) # a covariate
  time <- 10*rexp(n, exp(1 + x)) # time-to-event
  time1 <- floor(time)
  time2 <- ceiling(time)
  # Individuals are observed at discrete times
  # I observe (time1,time2) such that time1 <= time <= time2
  model <- pchreg(Surv(time1,time2, type = "interval2") ~ x, breaks = 10)
  
  
  
  
  # Try summary(model), predict(model)
  # See the documentation of predict.pch for more examples
}
\keyword{survival}
\keyword{models}
\keyword{regression}
