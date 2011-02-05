\name{nsem}
\alias{nsem}
\alias{coef.lava.nlin}
\alias{logLik.lava.nlin}
\alias{score.lava.nlin}
\alias{vcov.lava.nlin}
\alias{print.lava.nlin}
\title{Non-linear structural equation models}
\description{
  Fits non-linear SEM via 1. order Laplace Approximation
}
\usage{
nsem(model,data,laplace.control=list(lambda=0.3,niter=100,Dtol=1e-5),control=list(trace=1),vcov=TRUE,...)
}
\arguments{
  \item{model}{model list}
  \item{data}{data.frame}
  \item{laplace.control}{Options to the Laplace Approximation algorithm}
  \item{control}{Options parsed on to optimizer (nlminb)}
  \item{vcov}{Boolean indicating whether calculation of asymptotic
    standard errors should be disabled}
  \item{\dots}{Additional parameters parsed on to lower level functions}
}
\value{
  \code{lava.nlin} object (available methods: \code{coef},
  \code{logLik}, \code{score}, \code{vcov}, \code{print}, ...)
}
\details{
  The \code{model} arguments must be a named list:
  list( measure0, measure1, measure2, pred0, pred1, pred2, model)

  where 'model' currently can be either 'nsem2' (2-factor model) or
  'nsem3' (3-factor model). Each of the remaining elements can be either
  character vectors or formulas.

  Alternatively, a list of lists can be used as the \code{model}
  argument in which case a list of data.frames also much be given as
  \code{data} argument. A multigroup model is then fitted.
}
\examples{
\dontrun{
model1 <- list(measure1=~parc+supfc+antcin+postcin,measure2=~thS+putS+midbrainS,pred0=~age+bmi,model="nsem3")
b <- nsem(model1,data=dtrp)
}
}
\author{Klaus K. Holst}
\keyword{models}
\keyword{regression}
