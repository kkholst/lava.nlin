lap <- function(data,modelpar,
                control=list(niter=100,lambda=0.5,Dtol=1e-5),
                model,eb0=NULL) {
    if (is.null(eb0)) eb0 <- rbind(rep(0,modelpar$nlatent))
    arglist <- list(model,
                    data=data,
                    theta=modelpar$theta,
                    Sigma=modelpar$Sigma,
                    modelpar=modelpar,
                    control=control,
                    eb0=eb0,
                    DUP=FALSE,
                    PACKAGE="lava.nlin")
    res <- do.call(".Call",arglist)
    return(res)  
}

Lapl <- function(data,p,
                 modelpar,
                 model="nsem3",
                 regidx,indiv=FALSE,
                 ...) {
  if (missing(regidx)) {
    nvari <- unlist(sapply(c("nlatent","nvar"),function(x) grep(x,names(modelpar))))
    nvar <- sum(unlist(modelpar[nvari]))
    regidx <- 1:(length(p)-nvar)
  }
  modelpar$theta <- p[regidx]; modelpar$Sigma <- diag(exp(p[-regidx]))

  res <- lap(data,modelpar,model=model,...)
  if (!indiv) return(res$logLik)
  return(res[[1]])
}
