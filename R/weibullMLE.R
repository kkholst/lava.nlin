weibullmm <- function(formula, random, id, data,
                      init, iter=100, nsim=100, stepsize=1, update=FALSE,
                      ...) {
  aa <- c()
  on.exit(return(aa))
  require(survival)
  M <- model.frame(formula,data)
  Y <- model.extract(M, "response")
  if (!inherits(Y, "Surv")) {
    stop("Needs a Surv object")
  }    
  xf <- attributes(terms(random))$term.labels
  decomp.specials <- function(x,pattern="[()]") {
    st <- gsub(" ","",x)
    vars <- rev(unlist(strsplit(st,pattern)))[1]
    unlist(strsplit(vars,","))
  }
  decompran <- lapply(xf,decomp.specials)
  newf <- as.formula(paste("~",paste(unlist(lapply(decompran,function(x) x[1])),collapse="+")))
  Z <- model.frame(newf,data)
  restr <- unlist(lapply(decompran, function(x) {
    if (length(x)==1) return(NA)
    else x[2]
  }))
  X <- model.matrix(formula,data)[,-1,drop=FALSE]
  mydata <- cbind(Y[,1],Y[,2],X,Z)
  modelpar <- list(nlatent=as.double(ncol(Z)),nx=as.double(NCOL(X)),internal=TRUE)

  idx2 <- numeric(NCOL(Z))
  used <- list()
  curval <- 0
  for (i in 1:length(restr)) {
    if (is.na(restr[i])) {
      curval <- curval+1
      idx2[i] <- curval
    }
    else {
      if (restr[i]%in%names(used)) {
        idx2[i] <- used[[restr[i]]]
      } else {
        curval <- curval+1
        idx2[i] <- curval
        used[restr[i]] <- curval
      }
    }
  }
  theta.idx <- seq(length.out=2+modelpar$nx)
  modelpar$theta.idx <- c(theta.idx,unlist(idx2)+length(theta.idx))
  if (missing(init))
    init <- c(0,0,rep(1,length(unique(modelpar$theta.idx))-2))

  if (is.character(id)) {
    id <- data[,id]
  }
  newcl <- numeric(length(id))
  uparid <- unique(id)
  for (i in 1:length(uparid)) {
    newcl[which(id%in%uparid[i])] <- i
  }  
  aa <- StEM("weibullmm",modelpar=modelpar,init,data=as.matrix(mydata),cluster=newcl,nsim=nsim,iter=iter,update=update,stepsize=stepsize,...)
  return(aa)
}


Mstep_weibullmm <- function(cdata,modelpar,eta,...) {
  require(lava)
  cdata <- as.matrix(cdata)
  z_idx <- with(modelpar, (1:nlatent)+2+nx)
  eta_idx <- z_idx+modelpar$nlatent
  eta.all <- cdata[,eta_idx,drop=FALSE];
  colnames(eta) <- colnames(eta.all) <- paste("e",1:NCOL(eta),sep="")
  Z <- cdata[,z_idx]*eta.all
  newtheta.idx <- modelpar$theta.idx
  if (is.null(newtheta.idx)) {
    newtheta.idx <- 1:length(modelpar$theta)
  }
  newtheta.idx[z_idx] <- NA
  modelpar$theta <- modelpar$theta[na.omit(unique(newtheta.idx))]
  T <- cdata[,1]; status <- cdata[,2]
  if (modelpar$nx>0) Z <- cbind(cdata[,2+(1:modelpar$nx)],Z)
  est <- weibull.maximize(T,status,Z,theta.idx=newtheta.idx,...)
  pareta_idx <- with(modelpar,2+nx+1:nlatent)
  uvar.idx <- eta_idx
  if (!is.null(modelpar$theta.idx)) {
    var.idx <- modelpar$theta.idx[pareta_idx]
    uvar.idx <- na.omit(unique(var.idx))
  } 
  if (length(uvar.idx)<length(eta_idx)) {
    m <- lvm(,silent=TRUE)
    count <- 0
    for (i in var.idx) {
      count <- count+1
      var <- paste("e",count,sep="")
      covfix(m,var,var,silent=TRUE) <- paste("v",var.idx[count],sep="")
    }
    e <- estimate(m,eta,silent=TRUE)
    theta.sd <- diag(modelVar(e)$P)^0.5
  } else 
  theta.sd <- sd(eta)
  theta <- c(est[,1],theta.sd^2)
  res <- list(theta=theta)
  res
}


logl.weibull <- function(theta,T,status,X,theta.idx=NULL,indiv=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  if (missing(X)) {
    eta <- 0
  } else {
    beta <- theta[-c(1:2)]
    eta <- X%*%beta
  }  
  val <- status*log(lambda*p) + status*(p-1)*log(lambda*T) + status*eta - (lambda*T)^p*exp(eta)
  if (indiv)
    return(val)
  sum(val)
}
obj.weibull <- function(...) -logl.weibull(...)
score.weibull <- function(theta,T,status,X,theta.idx=NULL,indiv=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  lambdaT <- lambda*T
  loglambdaT <- log(lambdaT)
  lambdaTp <- exp(loglambdaT*p)
  
  if (missing(X)) {
    eta <- 0; expeta <- 1
    dbeta <- NULL
  } else {
    beta <- theta[-c(1:2)]
    eta <- X%*%beta
    expeta <- exp(eta)
    dbeta <- ((status-expeta*lambdaTp)%x%rbind(rep(1,NCOL(X))))*X
  }
  dp <- status*(1/p + loglambdaT) - loglambdaT*lambdaTp*expeta
  dlogp <- p*dp
  dlambda <- status*(p/lambda) - p*lambdaTp/lambda*expeta
  dloglambda <- lambda*dlambda
  S <- cbind(dloglambda,dlogp,dbeta)
  if (!is.null(theta.idx)) {
    u.idx <- na.omit(unique(theta.idx))
    newS <- matrix(0,ncol=length(u.idx),nrow=nrow(S))
    for (i in u.idx) {
      newS[,i] <- cbind(rowSums(S[,which(theta.idx==i),drop=FALSE]))
    }
    S <- newS
  }
  if (indiv)
    return(S)
  colSums(S)
}

hessian.weibull <- function(theta,T,status,X,theta.idx=NULL,all=FALSE) {
  if (!is.null(theta.idx)) {
    offsets <- which(is.na(theta.idx))
    theta <- theta[theta.idx]
    theta[offsets] <- 1
  }
  lambda <- exp(theta[1])
  p <- exp(theta[2])
  lambdaT <- lambda*T
  loglambdaT <- log(lambdaT)
  Tp <- T^p
  lambdaTp <- lambda^p*Tp
    if (missing(X)) {
    eta <- 0; expeta <- 1
    d2dlogpdbeta <- d2dlogpdbeta <- d2beta <- NULL
  } else {
    beta <- theta[-c(1:2)]
    eta <- X%*%beta
    expeta <- exp(eta)
    ## D(beta,beta)
    U <- ((expeta*Tp)%x%rbind(rep(1,NCOL(X))))*X
    d2beta <- -t(lambda^p*U)%*%X
    ## D(p,beta)
    d2dpdbeta <-  colSums(-((loglambdaT*lambdaTp*expeta)%x%rbind(rep(1,NCOL(X))))*X)
    d2dlogpdbeta <- d2dpdbeta*p  
    ## D(lambda,beta)
    d2dlambdadbeta <- -U*p*lambda^(p-1)
    d2dloglambdadbeta <- colSums(d2dlambdadbeta)*lambda  
  }
  ## D(p,p)
  dp <- status*(1/p + loglambdaT) - loglambdaT*lambdaTp*expeta
  d2p <- -sum(status/(p^2) + loglambdaT^2*expeta*lambdaTp)
  dlogp <- p*dp
  d2logp <- sum(dlogp)+p^2*d2p
  ## D(lambda,lambda)
  dlambda <- status*(p/lambda) - p*lambdaTp/lambda*expeta
  d2lambda <- -sum(status*(p/lambda^2) + p*(p-1)*lambdaTp/(lambda^2)*expeta)
  dloglambda <- lambda*dlambda
  d2loglambda <- sum(dloglambda)+lambda^2*d2lambda
  ## D(p,lambda)
  d2dpdlambda <- -status/(p^2) - loglambdaT^2*expeta*lambdaTp
  d2dpdlambda <- status/lambda - lambdaTp/lambda*expeta - p*loglambdaT*lambdaTp/lambda*expeta 
  d2dlogpdloglambda <- sum(d2dpdlambda)*p*lambda  
  ## Hessian:
  H <- matrix(0,length(theta),length(theta))
  H[1,1] <- d2loglambda
  H[2,2] <- d2logp
  H[1,2] <- H[2,1] <- d2dlogpdloglambda
  if (!missing(X)) {
    H[3:length(theta),3:length(theta)] <- d2beta
    H[2,3:length(theta)] <- H[3:length(theta),2] <- d2dlogpdbeta
    H[1,3:length(theta)] <- H[3:length(theta),1] <- d2dloglambdadbeta
  }
  if (!is.null(theta.idx)) {
    u.idx <- na.omit(unique(theta.idx))
    newH <- matrix(0,length(u.idx),length(u.idx))
    for (i in u.idx) {
      for (j in u.idx) {
        newH[i,j] <- sum(H[which(theta.idx==i),which(theta.idx==j)])
      }
    }
    H <- newH
  }
  if (all) {
    ## Score:
    if (missing(X)) dbeta <- NULL else dbeta <- status*X-lambda^p*U
    S <- cbind(dloglambda,dlogp,dbeta)
    if (!is.null(theta.idx)) {
      u.idx <- na.omit(unique(theta.idx))
      newS <- matrix(0,ncol=length(u.idx),nrow=nrow(S))
      for (i in u.idx) {
        newS[,i] <- cbind(rowSums(S[,which(theta.idx==i),drop=FALSE]))
      }
      S <- newS
    }
    attributes(H)$score <- S
    ## LogLik
    attributes(H)$logL <- sum(status*log(lambda*p) + status*(p-1)*loglambdaT + status*eta - lambdaTp*expeta)
  }
  return(H)
}

weibull.maximize <- function(T,status,X,
                             theta.idx=NULL,theta0,niter=100,tol1=1e-9,tol2=1e-9,lambda1=0.5,lambda2=1,trace=0,...) {

  if (missing(X)) {
    theta0 <- c(0,0)
  } else {
    if (missing(theta0))
      theta0 <- rep(0,ifelse(is.null(theta.idx),
                             2+NCOL(X),length(unique(na.omit(theta.idx)))))
  }
  
  thetas <- theta0; logL <- c()
  for (i in 1:niter) {
    H <- hessian.weibull(theta0,all=TRUE,T=T,status=status,X=X,theta.idx=theta.idx)
    S <- colSums(attributes(H)$score)
    gamma <- lambda2*sqrt((t(S)%*%S)[1])*diag(NROW(H))
    theta0 <- theta0 - lambda1*solve(H-gamma)%*%S
    thetas <- rbind(thetas,as.vector(theta0))
    logL <- c(logL, attributes(H)$logL)
    if(trace>0)
      if (i%%trace==0) {
        cat("Iter=",i, ", logLik=",attributes(H)$logL,"\n",sep="")
        cat("theta=(",paste(formatC(theta0),collapse=";"),")\n",sep="")
      }
    if (i>1)
      if (sum(abs(S^2))<tol1 & abs(logL[i]-logL[i-1])<tol2) break;
  }
  
  coefs <- cbind(theta0,diag(solve(-H))^0.5); colnames(coefs) <- c("Estimate","Std.Err");
  attributes(coefs)$score <- S
  attributes(coefs)$logLik <- tail(logL,1)
  mynames <- c("log(scale)","log(shape)")
  if (!missing(X)) {
    xnames <- colnames(X); if (is.null(xnames)) xnames <- paste("x",1:NCOL(X),sep="")
    mynames <- c(mynames,xnames)
  }
  if (!is.null(theta.idx)) {
    mynames <- mynames[na.omit(unique(theta.idx))]
  }
  rownames(coefs) <- mynames
  coefs
}
