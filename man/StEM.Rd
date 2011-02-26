\name{StEM}
\alias{StEM}
\title{Stochastic EM algorithm}
\description{
  Stochastic EM algorithm
}
\usage{
StEM(model, modelpar, theta0, data, cluster = 1:nrow(data), 
    eta, control = list(), Mfun, CondVarEta, update = FALSE, 
    m = 1, nsim = 50, iter = 50, stepsize = 1, burnin = nsim - 
        1, plot = FALSE, idx = 1:length(theta0), printidx = idx, 
    printvar = FALSE, ...) 
}
\arguments{
  \item{model}{String defining model}
  \item{modelpar}{Model parameters}
  \item{theta0}{Initial parameter vector}
  \item{data}{data-frame}
  \item{cluster}{See \code{Estep}}
  \item{eta}{See \code{Estep}}
  \item{control}{See \code{Estep}}
  \item{Mfun}{Function defining M-step}
  \item{CondVarEta}{Variance of proposal distribution (defaults to the
    identity matrix).}
  \item{update}{Adaptive algorithm, where the variance of the proposal
    distribution is updated in each E-step from the empirical
    distribution of the latent variables given the data (obtained from
    previous E-step).}
  \item{m}{The number of imputations to base the E-step on (m=1 => StEM)}
  \item{nsim}{Number of samples to draw in the E-step (a burnin period
    is needed!)}
  \item{iter}{Number of iterations of the EM-algorithm}
  \item{stepsize}{Scaling of the variance of the proposal distribution}
  \item{burnin}{For development testing only}
  \item{plot}{Plot coordinates 'idx' of chain}
  \item{idx}{Trace these coordinates}
  \item{printidx}{Print these coordinates in each M-step}
  \item{printvar}{Boolean indicating whether variance parameters should
    be printed}
  \item{\dots}{Additional parameters parsed on to lower-level functions.}
}
\value{
  list
}
\examples{
\dontrun{
modelpar <- list(nlatent=3, ny1=4, ny2=3, npred=2)
modelpar$model <- "nsem1"
nparreg <- with(modelpar, ny1+ny2 + (ny1-1) + (ny2-1) + 2 + npred)
nparvar <- with(modelpar, ny1+ny2+nlatent)
modelpar$theta <- rep(5,nparreg+nparvar)
npar <- with(modelpar, ny1+ny2 + (ny1-1) + (ny2-1) + npred + 2 + ny1 +
ny2 + nlatent)
aa <- StEM(modelpar$model,modelpar=modelpar,theta0=modelpar$theta,data=yy,cluster=1:NROW(yy),nsim=200,iter=500,plot=FALSE,update=FALSE,stepsize=0.2,idx=15:16,m=5)
plot(aa,idx=15:16)
bb <- lava.nlin::restart(aa,stepsize=0.1,nsim=1000,update=TRUE,m=1)
bb
}
}
\author{Klaus K. Holst}
\keyword{models}
