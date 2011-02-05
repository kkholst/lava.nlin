nsem <- function(model,
                 data,
                 laplace.control=list(lambda=0.3,niter=100,Dtol=1e-5),
                 control=list(trace=1),
                 vcov=TRUE,
                 ...) {

  require(numDeriv)  
  procmod <- function(M,data,...) {    
    if (class(M$measure0)[1]=="formula") M$measure0 <- all.vars(M$measure0)
    if (class(M$measure1)[1]=="formula") M$measure1 <- all.vars(M$measure1)
    if (class(M$measure2)[1]=="formula") M$measure2 <- all.vars(M$measure2)
    if (class(M$pred0)[1]=="formula") M$pred0 <- all.vars(M$pred0)
    if (class(M$pred1)[1]=="formula") M$pred1 <- all.vars(M$pred1)
    if (class(M$pred2)[1]=="formula") M$pred2 <- all.vars(M$pred2)
    mydata <- data[,c(M$measure0,M$measure1,M$measure2)]
    if (!is.null(M$pred0)) mydata <- cbind(mydata,data[,M$pred0])
    if (!is.null(M$pred1)) mydata <- cbind(mydata,data[,M$pred1])
    if (!is.null(M$pred2)) mydata <- cbind(mydata,data[,M$pred2])
    mydata <- na.omit(as.matrix(mydata))
    nn <- c(M$measure0,M$measure1,M$measure2)
    if (length(M$measure0)>0) nn <- c(nn,paste(M$measure0,"eta0",sep="<-"))
    if (length(M$measure1)>1) nn <- c(nn,paste(M$measure1[-1],"eta1",sep="<-"))
    if (length(M$measure2)>1) nn <- c(nn,paste(M$measure2[-1],"eta2",sep="<-"))
    if (length(M$pred0)>0) nn <- c(nn,paste("eta0",M$pred0,sep="<-"))
    if (length(M$pred1)>0) nn <- c(nn,paste("eta1",M$pred1,sep="<-"))
    if (length(M$pred2)>0) nn <- c(nn,paste("eta2",M$pred2,sep="<-"))
    nn <- c(nn,
            switch(M$model,
                   nsem3=c("eta2<-eta0","eta2<-eta0^2",
                          "eta0<->eta0","eta1<->eta1","eta2<->eta2"),
                   nsem2=c("eta2<-eta1","eta2<-eta1^2",
                     "eta1<->eta1","eta2<->eta2")))
    nlatent <- switch(M$model,
                      nsem3=3,
                      nsem2=2)
    if (length(M$measure0)>0) nn <- c(nn,paste(M$measure0,M$measure0,sep="<->"))
    if (length(M$measure1)>0) nn <- c(nn,paste(M$measure1,M$measure1,sep="<->"))
    if (length(M$measure2)>0) nn <- c(nn,paste(M$measure2,M$measure2,sep="<->"))

    mm <- list(nlatent=nlatent, nvar0=length(M$measure0), nvar1=length(M$measure1), nvar2=length(M$measure2), npred0=length(M$pred0), npred1=length(M$pred1), npred2=length(M$pred2))
    npar <- c()
    theta0 <- rep(0,with(mm,2*nvar0 + 2*(nvar1+nvar2-1) + npred0+npred1+npred2+
                         2+nvar0+nvar1+nvar2+nlatent))
    res <- list(data=mydata,names=nn,modelpar=mm,model=M$model)
    return(res)
  }

  models <- c()  
  if (is.list(model[[1]])) { ## multigroup
    for (i in 1:length(model)) {      
      models <- c(models, list(procmod(model[[i]],data[[i]])))
    }
  } else {
    models <- list(procmod(model,data))
  }
  allnames <- unique(unlist(lapply(models,function(x) x$names)))
  for (i in 1:length(models)) {
    models[[i]]$idx <- match(models[[i]]$names,allnames)
  }

  f <- function(p,...) { ## -log-likelihood
    val <- 0
    for (i in 1:length(models))
      val <- val + with(models[[i]],-Lapl(data,p[idx],modelpar,model=model,control=laplace.control))
    return(val)
  }
  
  theta0 <- rep(0,length(allnames)) ## Starting values
  if (!is.null(control$start)) {
    if (length(control$start)==length(theta0)) {      
      theta0 <- control$start
    } else {
      theta0[match(names(control$start),allnames)] <-
        control$start[which(control$start%in%allnames)]
    }    
  }

  res.Laplace <- tryCatch(nlminb(theta0,f,control=control),error=function(e) NULL)
  if (is.null(res.Laplace)) stop("Optimization error")

  if (vcov) {
    S0 <- numDeriv::grad(f,res.Laplace$par)
    H0 <- numDeriv::hessian(f,res.Laplace$par)
    vcov0 <- tryCatch(solve(H0),error=function(e) matrix(NA,ncol=ncol(H0),nrow=nrow(H0)))
  } else {
    S0 <- NULL
    vcov0 <- matrix(NA,nrow=length(theta0),ncol=length(theta0))
  }
  mycoef <- cbind(res.Laplace$par,diag(vcov0)^0.5)
  mycoef <- cbind(mycoef,mycoef[,1]/mycoef[,2],2*(1-pnorm(abs(mycoef[,1]/mycoef[,2]))))
  rownames(mycoef) <- allnames; colnames(mycoef) <- c("Estimate","Std.Err","Z-value","P(>|z|)")
  res <- list(coef=mycoef, vcov=vcov0, opt=res.Laplace, score=S0, data=data)
  class(res) <- "lava.nlin"
  return(res)
}


print.lava.nlin <- function(x,...) printCoefmat(x$coef)
vcov.lava.nlin <- function(object,...) object$vcov
score.lava.nlin <- function(x,...) x$score
coef.lava.nlin <- function(object,...) object$coef[,1]
logLik.lava.nlin <- function(object,...) {
  loglik <- -object$opt$objective
    if (is.null(attr(loglik, "nall"))) 
      attr(loglik, "nall") <- nrow(object$data)
  if (is.null(attr(loglik, "nobs"))) 
      attr(loglik, "nobs") <- nrow(object$data) - nrow(object$coef)
  if (is.null(attr(loglik, "df"))) 
      attr(loglik, "df") <- nrow(object$coef)
  class(loglik) <- "logLik"
  return(loglik)
}
