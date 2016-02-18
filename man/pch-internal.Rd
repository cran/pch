\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
CumSum(x)
myapply(x,fun)

makebreaks(y, d, x, breaks)
pois.newton(start, f, tol = 1e-5, maxit = 200, safeit = 0, \ldots)
pois.loglik(beta,d,x,w,off, zeror, deriv = 0)
poisfit(d,x,w,off)
pch.fit(z, y, d, x, w, breaks)

cleanlambda(lambda,x,rangex)
predF.pch(obj, x, y)
predQ.pch(obj, x, p)
sim.pch(obj, x, method = c("q", "s"))

build.splinex(x, method, df, degree, v)
predict.splinex(x, a)

\method{logLik}{pch}(object, \ldots)
\method{nobs}{pch}(object, \ldots)
\method{vcov}{pch}(object, \ldots)
\method{coef}{pch}(object, \ldots)
\method{print}{pch}(x, \ldots)
\method{summary}{pch}(object, \ldots)
\method{print.summary}{pch}(x, \ldots)
}
\keyword{internal}
