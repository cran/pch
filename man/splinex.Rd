\name{splinex}
\alias{splinex}
\title{
Including Interactions and Splines in Piecewise Constant Hazard Regression
}
\description{
This function can be used within a call to \code{\link{pchreg}} to automatically include 
spline functions in the linear predictor of the model.
}
\usage{
  splinex(method = c("ns", "bs"), df = 2, degree = 2, v = 0.98, \ldots)
}
\arguments{
  \item{method}{
  a character string indicating whether natural splines (\code{\link{ns}}) or B-splines (\code{\link{bs}})
  should be used. Default is \dQuote{\kbd{ns}}.
}
  \item{df}{
the degrees of freedom of the spline basis.
}
  \item{degree}{
the degree of the polynomial (only for \code{bs}).
}
  \item{v}{
a value between \kbd{0} and \kbd{1} determining how many principal components of the
design matrix must be used in model fitting (see \dQuote{Details}).
}
  \item{\ldots}{
for future arguments.
}
}
\details{
The piecewise constant hazard model implemented by \code{pchreg} can be used as
a nonparametric maximum likelihood estimator, in which the number of parameters
is allowed to increase with the sample size in order to achieve any desired flexibility.
Modeling the effect of covariates is as important as setting a sufficiently large
number of \code{breaks}.

By letting \code{splinex = splinex(\ldots)}, each column of the original design matrix 
is automatically replaced by the corresponding spline basis, defined by \code{method}, 
\code{df}, and \code{degree}.

This modeling approach has the drawback of generating a potentially large design matrix:
to reduce its dimension, you can choose \code{v < 1}. In this case, only the principal 
components explaining at least a proportion \code{v} of the variance
are used to fit the model (see \dQuote{Examples}).
}
\value{
The function returns its arguments, to be passed to an internal function \code{build.splinex}
that actually computes the design matrix.
}
\author{
Paolo Frumento <paolo.frumento@unipi.it>
}
\note{
  A multidimensional spline can be created by including a term like, for example
  \code{ns(x1,df)*ns(x2,df)}. This is not supported by \command{splinex},
  as it may generate a very large design matrix.
}
\seealso{
\code{\link{pchreg}}
}
\examples{

  require(splines)
  n <- 1000
  x1 <- runif(n,-2,2)
  x2 <- runif(n,-2,2)
  t <- rexp(n, exp(-abs(x1 - x2)))

  # a simple model
  model1 <- pchreg(Surv(t) ~ x1 + x2)

  # using splinex: the same as ~ ns(x1, df = 2) + ns(x2, df = 2)
  model2 <- pchreg(Surv(t) ~ x1 + x2, splinex = splinex("ns", v = 1))

  # include interaction: ~ ns(x1, df = 2) + ns(x2, df = 2) + ns(x1*x2, df = 2)
  model3 <- pchreg(Surv(t) ~ x1 * x2, splinex = splinex("ns", v = 1))

  # the same as model 3, only keep the PCs explaining at least 95 percent of the variance
  model4 <- pchreg(Surv(t) ~ x1 * x2, splinex = splinex("ns", v = 0.95))

  # true CDF vs fitted

  trueF <- pexp(t, exp(-abs(x1 - x2)))
  par(mfrow = c(2,2))
  plot(trueF, 1 - predict(model1)$Surv); abline(0,1, col = "red") # does not fit
  plot(trueF, 1 - predict(model2)$Surv); abline(0,1, col = "red") # neither
  plot(trueF, 1 - predict(model3)$Surv); abline(0,1, col = "red") # great!
  plot(trueF, 1 - predict(model4)$Surv); abline(0,1, col = "red") # almost as good
}
\keyword{regression}
