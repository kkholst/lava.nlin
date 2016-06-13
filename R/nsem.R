#' Non-linear structural equation models
#' 
#' Fits non-linear SEM via 1. order Laplace Approximation
#' 
#' The \code{model} arguments must be a named list: list( measure0, measure1,
#' measure2, pred0, pred1, pred2, model)
#' 
#' where 'model' currently can be either 'nsem2' (2-factor model) or 'nsem3'
#' (3-factor model). Each of the remaining elements can be either character
#' vectors or formulas.
#' 
#' Alternatively, a list of lists can be used as the \code{model} argument in
#' which case a list of data.frames also much be given as \code{data} argument.
#' A multigroup model is then fitted.
#' 
#' @aliases nsem coef.lava.nlin logLik.lava.nlin score.lava.nlin vcov.lava.nlin
#' print.lava.nlin
#' @param model model list
#' @param data data.frame
#' @param laplace.control Options to the Laplace Approximation algorithm
#' @param control Options parsed on to optimizer (nlminb)
#' @param vcov Boolean indicating whether calculation of asymptotic standard
#' errors should be disabled
#' @param \dots Additional parameters parsed on to lower level functions
#' @return \code{lava.nlin} object (available methods: \code{coef},
#' \code{logLik}, \code{score}, \code{vcov}, \code{print}, ...)
#' @author Klaus K. Holst
#' @keywords models regression
#' @examples
#' 
#' \dontrun{
#' model1 <- list(measure1=~parc+supfc+antcin+postcin,measure2=~thS+putS+midbrainS,pred0=~age+bmi,model="nsem3")
#' b <- nsem(model1,data=dtrp)
#' }
#' 
#' @export nsem
nsem <- function(model,
                 data,
                 laplace.control=list(),
                 control=list(trace=1),
                 method="ucminf",
                 vcov=TRUE,
                 p,fast=FALSE,
                 ...) {


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
    if (length(M$measure0)>0) nn <- c(nn,paste(M$measure0,"eta0",sep="~"))
    if (length(M$measure1)>1) nn <- c(nn,paste(M$measure1[-1],"eta1",sep="~"))
    if (length(M$measure2)>1) nn <- c(nn,paste(M$measure2[-1],"eta2",sep="~"))
    if (length(M$pred0)>0) nn <- c(nn,paste("eta0",M$pred0,sep="~"))
    if (length(M$pred1)>0) nn <- c(nn,paste("eta1",M$pred1,sep="~"))
    if (length(M$pred2)>0) nn <- c(nn,paste("eta2",M$pred2,sep="~"))
    nn <- c(nn,
            switch(M$model,
                   nsem3b=c("eta2~eta0","eta2~eta0^2",
                     "eta1~eta0","eta1~eta0^2",
                      "eta0,eta0",
                     "eta1,eta1","eta2,eta2"),
                   nsem3=c("eta2~eta0","eta2~eta0^2",
                          "eta0,eta0","eta1,eta1","eta2,eta2"),
                   nsem2=c("eta2~eta1","eta2~eta1^2",
                     "eta1,eta1","eta2,eta2")))
    omit <- switch(M$model,
                   nsem3b="eta0,eta0",
                   NULL)
    nlatent <- switch(M$model,
                      nsem3b=3,
                      nsem3=3,
                      nsem2=2)
    if (length(M$measure0)>0) nn <- c(nn,paste(M$measure0,M$measure0,sep=","))
    if (length(M$measure1)>0) nn <- c(nn,paste(M$measure1,M$measure1,sep=","))
    if (length(M$measure2)>0) nn <- c(nn,paste(M$measure2,M$measure2,sep=","))
    if (M$model=="nsem2") {
        nn[1] <- "eta1"
        nn[1+length(M$measure1)] <- "eta2"        
    }    
    mm <- list(nlatent=nlatent, nvar0=length(M$measure0), nvar1=length(M$measure1), nvar2=length(M$measure2), npred0=length(M$pred0), npred1=length(M$pred1), npred2=length(M$pred2))
##    npar <- c()
##    theta0 <- rep(0,with(mm,2*nvar0 + 2*(nvar1+nvar2-1) + npred0+npred1+npred2+
##                         2+nvar0+nvar1+nvar2+nlatent))
    res <- list(data=mydata,names=nn,modelpar=mm,model=M$model,omit=omit,
                measure0=M$measure0,measure1=M$measure1,measure2=M$measure2,
                pred0=M$pred0,pred1=M$pred1,pred2=M$pred2)
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
  allomit <- unique(unlist(lapply(models,function(x) x$omit)))
  pidx <- which(!(allnames%in%allomit))

  laplace.control0 <- list(lambda=0.5,niter=100,Dtol=1e-5,nq=0)
  laplace.control0[names(laplace.control)] <- laplace.control
  laplace.control <- laplace.control0

  f <- function(p,...) { ## -log-likelihood
    p0 <- rep(0,length(allnames))
    p0[pidx] <- p
    val <- 0
    for (i in 1:length(models))
      val <- val + with(models[[i]],-Lapl(data,p0[idx],modelpar,model=model,control=laplace.control))
    return(val)
  }
  if (!missing(p)) {
      logL <- f(p)
      return(logL)
  }
    
  
  theta0 <- rep(-0.2,length(setdiff(allnames,allomit))) ## Starting values
  if (!is.null(control$start)) {
    if (length(control$start)==length(theta0)) {      
      theta0 <- control$start
    } else {
      theta0[match(names(control$start),allnames)] <-
        control$start[which(control$start%in%allnames)]
    }    
  }
  control$start <- NULL
  control0 <- list(trace=1,kkt=FALSE)
  control0[names(control)] <- control
  control <- control0

  ##  res.Laplace <- tryCatch(nlminb(theta0,f,control=control),error=function(e) NULL)
  op <- tryCatch(optimx::optimx(theta0,f,method=method,control=control),error=function(e) NULL)
  res.Laplace <- list(par=as.vector(coef(op)), opt=op)
    
  if (is.null(res.Laplace)) stop("Optimization error")
    
  if (vcov) {
    S0 <- numDeriv::grad(f,res.Laplace$par)
    res.Laplace$grad <- S0     
    H0 <- numDeriv::hessian(f,res.Laplace$par)
    vcov0 <- tryCatch(solve(H0),error=function(e) matrix(NA,ncol=ncol(H0),nrow=nrow(H0)))
  } else {
    S0 <- NULL
    vcov0 <- matrix(NA,nrow=length(theta0),ncol=length(theta0))
  }
  
  colnames(vcov0) <- rownames(vcov0) <- setdiff(allnames,allomit)
  names(res.Laplace$par) <- setdiff(allnames,allomit)
  mycoef <- cbind(res.Laplace$par,diag(vcov0)^0.5)
  mycoef <- cbind(mycoef,mycoef[,1]/mycoef[,2],2*(1-pnorm(abs(mycoef[,1]/mycoef[,2]))))
  rownames(mycoef) <- setdiff(allnames,allomit); colnames(mycoef) <- c("Estimate","Std.Err","Z-value","P(>|z|)")
  res <- list(coef=mycoef, vcov=vcov0, opt=res.Laplace, score=S0, data=data, models=models)
  class(res) <- "lava.nlin"

    if (model$model=="nsem2" && length(models)==1 && !fast) {
        M <- models[[1]]
        m <- lvm()        
        m <- regression(m,from="eta1",to=M$measure1)
        m <- regression(m,from="eta2",to=M$measure2)
        m <- regression(m,from=c("eta1","eta1^2"),to="eta2")
        m <- regression(m,from=M$pred1,to="eta1")
        m <- regression(m,from=M$pred2,to="eta2")
        intercept(m,c("eta1^2",M$measure1[1],M$measure2[1])) <- 0
        covariance(m,"eta1^2") <- 0
        regression(m,from="eta1",to=M$measure1[1]) <- 1
        regression(m,from="eta2",to=M$measure2[1]) <- 1
        latent(m) <- c("eta1","eta2","eta1^2")
        transform(m, y="eta1^2",x="eta1") <- function(x) x^2
        labels(m) <- c("eta1"=expression(eta[1]),"eta2"=expression(eta[2]),"eta1^2"=expression(eta[1]^2))
        pn <- names(pars(res))
        ## pn <- gsub("\\^","",pn)
        pm <- coef(m)
        m$order <- match(pn,pm)
        ee <- structure(c(res, list(model=m,data=data)),class=c("lava.nlin","lvmfit"))
        ee$coef <- ee$coef[match(pm,pn),]
        p <- length(M$measure1)+length(M$measure2)
        rownames(ee$coef)[seq(p+1,length(pm))] <- paste0("p",seq(length(pm)-p))
        rownames(ee$coef)[seq(p)] <- paste0("m",seq(p))
        return(ee)
    }
    
  return(res)
}


plot.lava.nlin <- function(x,diag=FALSE,f,...) {
    if (inherits(x,"lvmfit") & missing(f)) {
        NextMethod("plot",diag=diag,...)
    } else {
        return(plot(estimate(NULL,coef=pars(x),vcov=vcov(x)),f=f,...))
    }
}

print.lava.nlin <- function(x,...) {
    if (inherits(x,"lvmfit")) {
        return(NextMethod("print",...))
    }
    printCoefmat(x$coef)
    return(invisible(x))
}


vcov.lava.nlin <- function(object,...) object$vcov

score.lava.nlin <- function(x,...) x$score

pars.lava.nlin <- function(x,...) x$opt$par
    
logLik.lava.nlin <- function(object,...) {
  loglik <- -object$opt$opt["value"]
    if (is.null(attr(loglik, "nall"))) 
      attr(loglik, "nall") <- nrow(object$data)
  if (is.null(attr(loglik, "nobs"))) 
      attr(loglik, "nobs") <- nrow(object$data) - length(pars(object))
  if (is.null(attr(loglik, "df"))) 
      attr(loglik, "df") <- length(pars(object))
  class(loglik) <- "logLik"
  return(loglik)
}
