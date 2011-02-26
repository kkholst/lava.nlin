\name{Estep}
\alias{Estep}
\title{Metropolis-Hastings sampler}
\description{
  General Metropolis-Hastings sampler. 
}
\usage{
Estep(modelpar, theta = modelpar$theta, data = matrix(), 
    cluster = 1:NROW(data), eta = with(modelpar, matrix(0, max(cluster), 
        nlatent)), CondVarEta = with(modelpar, diag(rep(1, nlatent), 
        ncol = nlatent)), fulldata = FALSE, keep = 1, control = list(), 
    ...) 
}
\arguments{
  \item{modelpar}{Model parameters}
  \item{theta}{Initial parameter}
  \item{data}{data}
  \item{cluster}{cluster indicator}
  \item{eta}{Initial state of latent variable (starting point of Markov Chain)}
  \item{CondVarEta}{Conditional variance of eta given data (i.e. the
    variance of the multivariate normal proposal distribution). Defaults
  to the unit matrix.}
  \item{fulldata}{For internal use only}
  \item{keep}{For internal use only}
  \item{control}{Parameters for the MH algorithm}
  \item{\dots}{Additional parameters parsed to lower-level functions}
}
\value{
  List with four members. 'accept' gives the number of accepted
  jumps. 'chain': a matrix containing the simulated markov chain.
  'eta': The last 'k' observations from the chain (unless thin was
  defined in the control parameters). 'data': the data and 'eta' joined
  together. 
}
\details{
  control parameters. 'nsim': Number of samples. 'thin': amount of
  thinning. 'stepsize': scale parameter of the multivariate normal
  proposal distribution...
}
\examples{
## Primary use is as a sampler for the StEM-function (Stochastic EM
## algorithm)
## For faster sampling the evaluation of the target distribution should
## be implemented in C/C++.
## Here a small example, where we sample from a bivariate normal-distribution
## without specifying the normalizing constant

normMH <- function(data,eta,modelpar,...) {
  m <- modelpar$theta[1:2]
  W <- matrix(modelpar$theta[3:6],ncol=2)
  -0.5*(eta-m)\%*\%W\%*\%t(eta-m)
}

mu <- c(0,0) ## mean
S <- matrix(c(1,-0.5,-0.5,2),ncol=2) ## Variance
theta <- c(mu,as.vector(solve(S)))

e <- Estep(modelpar=list(model="normMH",nlatent=2,theta=theta,internal=FALSE),eta=cbind(10,10),control=list(nsim=10000,stepsize=2))

print(e$accept/10000) ## Acceptance ratio

\donttest{
X <- e$chain
plot(as.mcmc(X)) ## Convergence of chain?

library(MASS)
eqscplot(e$chain,type="l")
points(X,col=heat.colors(nrow(X)),cex=0.5)
print(colMeans(e$chain[-(1:(nsim/2)),]))
print(var(e$chain[-(1:(nsim/2)),]))
}

## An a multivariate t-distribution (df=3)
tMH <- function(data,eta,modelpar,...) {
  m <- modelpar$theta[1:2]
  W <- matrix(modelpar$theta[3:6],ncol=2)
  p <- 2; df <- theta[7]
  -((df+p)/2)*log(1+1/df*(eta-m)\%*\%W\%*\%t(eta-m))
}
theta <- c(mu,as.vector(solve(S)),df=3)
e <- Estep(modelpar=list(model="tMH",nlatent=2,theta=theta,internal=FALSE),eta=cbind(10,10),control=list(nsim=10000,stepsize=2))
print(e$accept/10000) ## Acceptance ratio

}
\author{Klaus K. Holst}
\keyword{models}
