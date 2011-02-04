lap <- function(data,modelpar,control=list(niter=100,lambda=0.5,Dtol=1e-5),
                model
                ) {
  arglist <- list(name=model,
                  data=data,
                  theta=modelpar$theta,
                  Sigma=modelpar$Sigma,
                  modelpar=modelpar,
                  control=control,
                  DUP=FALSE,
                  PACKAGE="lava.nlin")
  res <- do.call(".Call",arglist)
  return(res)  
}

Lapl <- function(data,p,
                 modelpar,
                 model="nsem3",
                 regidx,...) {
  if (missing(regidx)) {
    nvari <- unlist(sapply(c("nlatent","nvar"),function(x) grep(x,names(modelpar))))
    nvar <- sum(unlist(modelpar[nvari]))
    regidx <- 1:(length(p)-nvar)
  }
  modelpar$theta <- p[regidx]; modelpar$Sigma <- diag(exp(p[-regidx]))
  lap(data,modelpar,model=model,...)$logLik
}
