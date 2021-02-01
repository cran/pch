\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
CumSum(x)
myapply(x,fun)

makebreaks.ct(y, d, x, w, breaks)
makebreaks.ic(y, d, x, w, breaks)

pois.loglik(beta,d,x,w,off, zeror, deriv = 0)
icpch.loglik(beta, y, x,w,U1,U2, zeror, deriv = 0, final = FALSE)
convert.Surv(y)
break.y(y, breaks)

pois.fit(d,x,w,off)
pch.fit.ct(z, y, d, x, w, breaks)
pch.fit.ic(z, y, d, x, w, breaks)

newton(start, f, tol = 1e-5, maxit = 200, safeit = 0, \ldots)

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
