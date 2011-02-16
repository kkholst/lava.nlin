##source("Mstep.R")
## mylib <- paste("mh2", .Platform$dynlib.ext, sep = "")
## if (is.loaded("MH"))
##   dyn.unload(mylib)
## try(dyn.load(mylib))  
if (!require(coda)) print("ups no coda")


###{{{ Eval

Eval <- function(modelpar,eta,data,cluster=1:NROW(data)) {
  arglist <- list(name="Eval",
                  modelpar=modelpar,
                  eta=eta,
                  data=data,
                  cluster=cluster,
                  DUP=FALSE)
  res <- do.call(".Call",arglist)
  return(res$x)
}

###}}} Eval

###{{{ as.mcmc

as.mcmc.StEM <- function(x,...) {
  require(coda)
  ##  theta <- lapply(x, function(y) y[["theta"]])
  ##  draws <- as.mcmc.list(lapply(theta, function(x) mcmc(x)))
  ##  return(draws)
  return(as.mcmc(x[["theta"]]))
}

###}}} as.mcmc

###{{{ window

window.StEM <- function(x,start,end=nrow(x$theta),n,...) {
  if (!missing(n)) {
    if (missing(start)) {
      idx <- seq(end-n+1,end)
    } else {
      idx <- seq(1,n)+start
    }
  } else {
    if (missing(start)) {
      start <- 1
    } 
    idx <- seq(start,end)
  }
  x$theta <- x$theta[idx,]
  return(x)
}

###}}} window

###{{{ print

print.StEM <- function(x,burnin=0,...) {
  with(x, {cat("Model: '",model,"' with ", modelpar$nlatent, " latent variable",sep="");
              if(modelpar$nlatent>1) cat("s",sep=""); cat(".\n", sep="") })
  with(x, cat(nrow(theta), " StEM iterations with an E-step based on ", control$nsim, " MCMC samples. Average acceptance rate: ", formatC(mean(mc$accept)/control$nsim),"\n", sep=""))
  cat("Estimated mean parameter:\n");
  res <- coef(x,burnin,...)
  print(res)
  cat("\n")
##  cat("Variance parameters:\n")
##  varpar <- coef(x,burnin,var=TRUE,...)
##  varpar <- with(x$modelpar, diag(matrix(varpar,nvar,nvar))) 
##  print(varpar, quote=FALSE)
  invisible(x)
}

###}}} print

###{{{ plot
plot.StEM <- function(x,idx,start=1,end=nrow(x$theta),coda=FALSE,lwd=2,xlab="Iteration",ylab="Parameter value",...) {
  mywin <- seq(start,end)
  if (coda) {
    require(coda)
    m <- window(as.mcmc(x),start=start,end=end)
    plot(m[,idx],...)
    return(invisible(x))
  }
  if (missing(idx)) {
    matplot(x$theta[mywin,], type="l",lwd=lwd,ylab=ylab,xlab=xlab,...)
    return(invisible(x))
  }
  matplot(x$theta[mywin,idx], type="l",lwd=lwd,ylab=ylab,xlab=xlab,...)
  invisible(x)
}
###}}} plot

###{{{ coef
coef.StEM <- function(object,burnin=0,var=FALSE,vardiag=FALSE,both=FALSE,...) {
  if (both) {
    res <- c(coef(object,burnin=burnin,...),coef(object,burnin=burnin,vardiag=TRUE,...))
    return(res)
  }
  if (burnin>0) {
    res <- colMeans(object$theta[-c(1:burnin),,drop=FALSE])
  } else {
    res <- colMeans(object$theta)
  }
  if (var | vardiag) {
    if (burnin>0) {
      res <- colMeans(object$Sigma[-c(1:burnin),,drop=FALSE])
    } else {
      res <- colMeans(object$Sigma)
    }
    if (vardiag) {
      res <- with(object$modelpar, diag(matrix(res,nvar,nvar)))
    }    
  }
  return(res)
}
###}}} coef

###{{{ sim

sim.StEM <- function(x,n=5000,burnin=min(NROW(x$theta)-1,100),theta=coef(x,burnin), control=list(), eta=x$mc$eta, CondVarEta=var(x$mc$eta), ...) {
  
  if (is.null(CondVarEta)) {
    CondVarEta <- with(x$modelpar, diag(rep(1,nlatent),ncol=nlatent))
  } 
  mymodelpar <- x$modelpar
  mymodelpar$theta <- theta[mymodelpar$theta.idx]
  opts <- x$control
  opts[names(control)] <- control
  opts$nsim <- n  
  
##  arglist <- with(x, list(modelpar=mymodelpar,data=data,eta=mc$eta,CondVarEta=CondVarEta,control=opts,cluster=cluster))  
##  draw <- do.call("Estep", arglist)  
  draw <- Estep(modelpar=mymodelpar, mymodelpar$theta,
                data=x$data, cluster=x$cluster,
                eta=eta,
                CondVarEta=CondVarEta,
                control=opts,...)
  return(draw)
}

###}}} sim

###{{{ merge

merge.StEM <- function(x,...) {
  objects <- list(x,...)
  if (length(objects)<2)
    return(x)
  K <- length(objects)
  control <- list(x$control)
  theta <- objects[[1]]$theta
  Sigma <- objects[[1]]$Sigma
  for (l in objects[-1]) {
    theta <- rbind(theta,l$theta)
    Sigma <- rbind(Sigma,l$Sigma)
    control <- c(control, list(l$control))
  }
  res <- objects[[K]]
  res$theta <- theta
  res$Sigma <- Sigma
  res$mergeinfo <- list(K=K+1, control=control)
  ##class(res) <- "StEM"
  return(res)
}

###}}} merge

###{{{ restart

"restart" <- function(x,...) UseMethod("restart")
restart.StEM <- function(x,iter=5000,theta=tail(x$theta,1),Sigma=tail(x$Sigma,1),nsim,burnin=nsim-1,stepsize,CondVarEta, Mfun=x$Mfun, ...) {
  res <- NULL
  on.exit(      
          return(res)
          )
  mycontrol <- x$control
  mycontrol$iter <- iter
  if (!missing(nsim)) {
    mycontrol$nsim <- nsim
    mycontrol$burnin <- burnin
  }
  if (!missing(stepsize)) {
    mycontrol$stepsize <- stepsize
  }
  if (!missing(stepsize)) {
    mycontrol$stepsize <- stepsize
  }
  MyCondVarEta <- x$CondVarEta
  if (!missing(CondVarEta))
    MyCondVarEta <- CondVarEta
  Mytheta <- theta
  MySigma <- Sigma
  MyMfun <- Mfun
##  browser()
  res <- with(x, StEM(model=model,
                      modelpar=modelpar,
                      theta0=Mytheta, Sigma=MySigma,
                      eta=mc$eta,
                      cluster=x$cluster,
                      CondVarEta=MyCondVarEta,
                      data=data,
                      Mfun=MyMfun,
                      control=mycontrol, ...)
              )
  return(res)
}

###}}} restart

###{{{ StEM

StEM <- function(model,
                 modelpar,
                 theta0,data,cluster=1:nrow(data),eta,
                 ##Sigma,
                 control=list(),Mfun,
                 CondVarEta, update=FALSE, m=1,nsim=50,iter=50,stepsize=1,burnin=nsim-1,
                 plot=FALSE,idx=1:length(theta0),printidx=idx,printvar=FALSE,...) {
  res <- mc <- Var <- cdata <- c()
  mycall <- match.call()
  opts <- list(iter=iter,
               nsim=nsim,
               burnin=burnin,
               thin=0,
               m=m,
               stepsize=stepsize,
               delta.var=20,
               verbose=1)  
  opts[names(control)] <- control
  n <- nrow(data)
  mymodelpar <- list(model=model,
                     npred=1,
                     nlatent=1,
                     internal=TRUE,
                     theta.idx=NULL
##                     nvar=4
                     )
  mymodelpar[names(modelpar)] <- modelpar
  if (!missing(theta0)) mymodelpar$theta <- theta0
  nlatent <- mymodelpar$nlatent

  on.exit(
          {
            val <- list(theta=res, ##Sigma=Var,
                        mc=mc, cdata=cdata, modelpar=mymodelpar, model=model, Mfun=Mfun, control=opts, call=mycall, data=data, CondVarEta=EtaVar, cluster=cluster);
            class(val) <- "StEM";
            return(val)
          }
          )
    
  if (missing(eta)) {
    eta <- matrix(0,nrow=max(cluster),ncol=nlatent)
  }

  cur <- mymodelpar$theta
  cur.idx <- 1:length(cur)
  if (!is.null(mymodelpar$theta.idx)) {
    cur <- cur[mymodelpar$theta.idx]
    u.idx <- unique(modelpar$theta.idx)
    cur.idx <- sapply(u.idx, function(x) which(modelpar$theta.idx%in%x)[1])
  }

  
  if(opts$verbose>0)
    cat("iter=0; theta=",paste(formatC(cur[cur.idx]),collapse=" "),"\n",sep="")
##  if (missing(Sigma))
##    Sigma <- diag(mymodelpar$nvar)
##  Sigma <- as.double(Sigma)

  mypars <- list(...)
  if (plot) {
    if (is.null(mypars$xlab))
      mypars$xlab <- "Iterations"
    if (is.null(mypars$ylab))
      mypars$ylab <- "Parameter values"
    if (is.null(mypars$main))
      mypars$main <- "EM algorithm"
    if (is.null(mypars$xlim))
       mypars$xlim <- c(1,opts$iter);
    plotpars <- mypars
##    plotpars[["type"]] <- NULL
    plotpars$type <- "n"
    plotpars$x <- plotpars$y <- 0
    do.call("plot",plotpars)
    if (is.null(mypars$col))
      mypars$col <- 1:length(idx);
  }

  if (is.character(model) & missing(Mfun))
    eval(parse(text=paste("Mfun <- Mstep_",model,sep="")))

  EtaVar <- diag(rep(1,nlatent),ncol=nlatent)
  if (!missing(CondVarEta))
    EtaVar <- CondVarEta

  ee <- list(theta=c(1,1,1,1))
  for (i in 1:opts$iter) {    
    if (plot) {
      count <- 0
      for (ii in idx) {
        count <- count+1
        plotpars <- mypars
        plotpars$col <- plotpars$col[count]
        plotpars$x <- i
##        if (ii>length(cur))
##          plotpars$y <- Sigma[ii-length(cur)]
##        else
        plotpars$y <- cur[ii]
        do.call("points",plotpars)
      }
    }

##    Var <- rbind(Var, as.vector(Sigma))
    res <- rbind(res, cur[cur.idx])
##    print("EtaVar=")
##    print(EtaVar)
    ##    mymodelpar$theta.var <- as.double(Sigma)
    mc <- Estep(modelpar=mymodelpar, cur,
                data=data, cluster=cluster,
                eta=eta,
                CondVarEta=EtaVar,
                control=opts)
    
    ##    print(mc$mean)
    ##    print(mc$var)

##     if (m==1) {
##       cdata <- mc$data
## ##      names(cdata)[1:nlatent] <- paste("eta",0:(nlatent-1),sep="")
##     } else {
##       d <- dim(mc[[1]])
##       midx <- seq(from=d[1],by=-1, length.out=m)
      
##       draws <- data.frame(mc$chain[midx,,drop=FALSE])
##       if (opts$thin>1)
##         draws <- Thin(draws,opts$thin)
##       require(Hmisc)          
##       etadata<- reShape(draws, base=paste("eta",1:nlatent,".",sep=""), reps=n)
##       names(etadata) <- c("seqno",paste("eta",0:(nlatent-1),sep=""))
##       ##    names(etadata) <- c(paste("eta",0:(mymodelpar$nlatent-1),sep=""))
##       cdata <- cbind(etadata, t(sapply(etadata$seqno,function(i) unlist(data[i,]))))
##     }
    cdata <- mc$data
##    browser()
    newpars <- Mfun(cdata,mymodelpar,mc$eta)
    ##    etamed <- apply(mc$chain,2,median)
    ##    opts$eta <- matrix(etamed, ncol=eta)
    ##opts$eta <- matrix(tail(mc$chain,1),ncol=mymodelpar$nlatent)
    eta <- mc$eta
    cur <- newpars$theta; ##Sigma <- newpars$theta.var
    if (missing(CondVarEta) | update) {  
      ##    EtaVar <- mc$var
      ##     EtaVar <- var(etadata[,-1])
      EtaVar <- var(mc$eta)
    }
  
    if (opts$verbose!=0 & i%%opts$verbose==0) {
      cat("iter=",i, "; Average a.ratio=", formatC(mean(mc$accept)/opts$nsim), "\n\ttheta=",paste(formatC(cur[cur.idx][printidx]),collapse=" "),"\n",sep="")
##      if (printvar)
##        cat("\tvar=", paste(formatC(diag(matrix(Sigma,mymodelpar$nvar,mymodelpar$nvar))),collapse=" "),"\n",sep="")
    }
  }
}

###}}} StEM

###{{{ Estep

Estep <- function(modelpar,
                  theta=modelpar$theta,
                  data=matrix(),
                  cluster=1:NROW(data),
                  eta=with(modelpar,matrix(0,max(cluster),nlatent)),
                  CondVarEta=with(modelpar, diag(rep(1,nlatent),ncol=nlatent)),
                  fulldata=FALSE,
                  keep=1,
                  control=list(),...) {
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
  mymodelpar <- list(internal=TRUE)
  mymodelpar[names(modelpar)] <- modelpar; modelpar <- mymodelpar
  modelpar$theta <- theta
  n <- nrow(data); k <- ncol(data)
  if (length(eta)==1) {
    modelpar$nlatent <- eta
    eta <- matrix(0,nrow=max(cluster),ncol=eta)
  }
  mycontrol <- list(stepsize=1,nsim=100,burnin=0,m=1,thin=0)
  mycontrol[names(control)] <- control
  ndata <- ncol(data)
  arglist <- list(name="MH",
                  data=data,
                  cluster=cluster,
                  etainit=eta,
                  Sigma=CondVarEta,
                  modelpar=modelpar,
                  control=mycontrol,
                  DUP=FALSE)
  save(arglist, file="arglist.rda")
  res <- do.call(".Call",arglist)

##  return(res)
##  return(res)
  ##                accept=res[[2]][idx,,drop=FALSE])
  ##res$p <- with(res, colSums(accept)/nrow(accept))
  ##  res$p <- res[[2]]
  
  ##  print(res$chain[1,])
##  onedraw <- as.data.frame(res$chain[nrow(res$chain),,drop=FALSE])
##   print(class(onedraw))
##   print(names(onedraw))
##   print(modelpar$nlatent)
##  mydraw <- reShape(onedraw, base=paste("eta",1:modelpar$nlatent,".",sep=""), reps=nrow(data))[,-1]
##  colnames(mydraw) <- paste("eta",1:ncol(mydraw),sep="")
  
  ##  res <- lapply(res,function(x) as.data.frame(x[-c(1:burnin),]))
  ##  return(res)
  if (fulldata) {
    require(Hmisc)
    cdata<- reShape(res$chain, base=paste("eta",1:ncol(eta),".",sep=""), reps=nrow(data))[,-1]
    names(cdata) <- c(paste("eta",0:(ncol(eta)-1),sep=""))
    res$mean <- colMeans(cdata)
    res$var <- var(cdata)
    mydata <- as.data.frame(data)
    cdata <- cbind(cdata, t(sapply(cdata$seqno,function(i) unlist(mydata[i,]))))
    res$cdata <- cdata
  }
  res
}

###}}} Estep


