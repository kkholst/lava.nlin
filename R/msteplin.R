h <- function(eta,data,clustersize=1,modelpar,...) {
  theta <- modelpar$theta
  ny <- with(modelpar, nvar-nlatent)
  ny1 <- modelpar$nvar1
  ny2 <- ny-ny1
  npred <- modelpar$npred
  mu1 <- theta[1:ny1]
  mu2 <- theta[(ny1+1):ny]
  lambda1 <- c(1,theta[1:(ny1-1) + ny])
  lambda2 <- c(1,theta[ny1:(ny-2) +ny])
  beta <- theta[1:npred + (2*ny-2)]
  gamma <- tail(theta,1)
  res <- c()
  for (i in 1:ny1) {
    res <- c(res, data[i]-mu1[i]-lambda1[i]*eta[2])
  }
  for (i in 1:ny2) {
    res <- c(res, data[i+ny1]-mu2[i]-lambda2[i]*eta[3])
  }
  zeta2 <- eta[3] - gamma[1]*eta[1]
  zeta1 <- eta[2] - eta[1]
  zeta0 <- eta[1]
  for (i in 1:npred)
    zeta0 <- zeta0 - beta[i]*data[i+ny]  
  unlist(c(zeta0=zeta0,zeta1=zeta1,zeta2=zeta2,res))
}


Mstep_lin1 <- function(cdata,modelpar) {
  nvar <- modelpar$nvar
  npred <- modelpar$npred
  ny1 <- modelpar$nvar1
  ny <- nvar-3
  ny2 <- ny-ny1;


  
  f <- paste("eta0 ~ -1")
  if (npred>0) {
    f <- paste(f, paste("x",1:npred, sep="", collapse="+"), sep="+")
  }
  l.eta0 <- lm(f, cdata)
  l.eta1 <- lm(eta1 ~ -1 + eta0, cdata)
  l.eta2 <- lm(eta2 ~ -1 + eta0, cdata)
  
  ## cdata$eta00 <- cdata$eta0 - coef(l.eta0)[1]
  ## cdata$eta10 <- cdata$eta1 - coef(l.eta1)[1]
  ## cdata$eta20 <- cdata$eta2 - coef(l.eta2)[1]

  cdata$eta00 <- cdata$eta0
  cdata$eta10 <- cdata$eta1
  cdata$eta20 <- cdata$eta2

#  gamma0 <- coef(l.eta0)[-1]
#  gamma1 <- coef(l.eta1)[-1]
#  gamma2 <- coef(l.eta2)[-1]
  gamma0 <- coef(l.eta0)[]
  gamma1 <- coef(l.eta1)[]
  gamma2 <- coef(l.eta2)[]


  l1 <- list()
  for (i in 1:ny1) {
    f <- paste("y",i," ~ eta10",sep="")
    l1 <- c(l1, list(lm(f, cdata)))
  }
  l2 <- list()
  for (i in (1:ny2 + ny1)) {
    f <- paste("y",i," ~ eta20",sep="")
    l2 <- c(l2, list(lm(f, cdata)))    
  }
  c1 <- matrix(unlist(lapply(l1,coef)),ncol=2,byrow=TRUE)
  c2 <- matrix(unlist(lapply(l2,coef)),ncol=2,byrow=TRUE)

  
  mu <- c(c1[,1], c2[,1])

  ## Adjusting for fixed parameter
  b1 <- c1[1,2] 
  b2 <- c2[1,2]
###  theta.y <- c(c1[-1,2]/b1,c2[-1,2]/b2)
  theta.y <- c(c1[-1,2],c2[-1,2])
###  theta.eta <- c(gamma0*b1*gamma1, gamma2*b2/(b1*gamma1))
  theta.eta <- c(gamma0, gamma2)
  theta <- c(mu, theta.y, theta.eta)
  ## Variance parameters

##  cat(".......\n")
  N <- nrow(cdata)
  var.y <- unlist(lapply(c(l1,l2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- unlist(lapply(list(l.eta0,l.eta1,l.eta2),function(x) t(residuals(x))%*%residuals(x)/N))
##  sd.y <- unlist(c(lapply(c(l1,l2), function(x) summary(x)$sigma)))  
##  sd.eta <- unlist(c(lapply(list(l.eta0,l.eta1,l.eta2), function(x) summary(x)$sigma)))
##   sd.eta <- sd.eta*c(gamma1*b1,b1,b2)
##  sd. <- c(sd.eta, sd.y)
  ##  Sigma <- diag(sd.^2)
  ##  var.eta <- var.eta*(c(gamma1*b1,b1,b2)^2)
  Sigma <- diag(c(var.eta,var.y))
##  cat("ooost\n")
  
  ## cdata2 <- cdata
  ## cdata2$eta1 <- cdata$eta1*b1
  ## lm(y1 ~ eta1, cdata2)
  ## cdata2$eta0 <- cdata$eta0*gamma1*b1
  ## lm(eta1 ~ -1 + eta0, cdata2)
  ## cdata2$eta2 <- cdata$eta2*b2
  ## lm(y4 ~ eta2, cdata2)
  ## lm(eta2 ~ eta0-1, cdata2)
  ## diag(Sigma)[1:3]
  ## sd.eta^2*c(gamma1^2*b1^2,b1^2,b2^2)
  
  newmodelpar <- modelpar
  newmodelpar$theta <- theta; newmodelpar$theta.var <- as.double(Sigma)
  return(newmodelpar)
}

Mstep_lin1a <- function(cdata,modelpar) {
  nvar <- modelpar$nvar
  npred <- modelpar$npred
  ny1 <- modelpar$nvar1
  ny <- nvar-3
  ny2 <- ny-ny1;
  
  f <- paste("eta0 ~ 1")
  if (npred>0) {
    f <- paste(f, paste("x",1:npred, sep="", collapse="+"), sep="+")
  }
  l.eta0 <- lm(f, cdata)
  l.eta1 <- lm(eta1 ~ 1 + eta0, cdata)
  l.eta2 <- lm(eta2 ~ 1 + eta0, cdata)
  
  cdata$eta00 <- cdata$eta0 - coef(l.eta0)[1]
  cdata$eta10 <- cdata$eta1 - coef(l.eta1)[1]
  cdata$eta20 <- cdata$eta2 - coef(l.eta2)[1]

  gamma0 <- coef(l.eta0)[-1]
  gamma1 <- coef(l.eta1)[-1]
  gamma2 <- coef(l.eta2)[-1]
#  gamma0 <- coef(l.eta0)[]
#  gamma1 <- coef(l.eta1)[]
#  gamma2 <- coef(l.eta2)[]


  l1 <- list()
  for (i in 1:ny1) {
    f <- paste("y",i," ~ eta10",sep="")
    l1 <- c(l1, list(lm(f, cdata)))
  }
  l2 <- list()
  for (i in (1:ny2 + ny1)) {
    f <- paste("y",i," ~ eta20",sep="")
    l2 <- c(l2, list(lm(f, cdata)))    
  }
  c1 <- matrix(unlist(lapply(l1,coef)),ncol=2,byrow=TRUE)
  c2 <- matrix(unlist(lapply(l2,coef)),ncol=2,byrow=TRUE)

  
  mu <- c(c1[,1], c2[,1])

  ## Adjusting for fixed parameter
  b1 <- c1[1,2] 
  b2 <- c2[1,2]
###  theta.y <- c(c1[-1,2]/b1,c2[-1,2]/b2)
  theta.y <- c(c1[-1,2],c2[-1,2])
###  theta.eta <- c(gamma0*b1*gamma1, gamma2*b2/(b1*gamma1))
  theta.eta <- c(gamma0, gamma2)
  theta <- c(mu, theta.y, theta.eta)
  ## Variance parameters

##  cat(".......\n")
  N <- nrow(cdata)
  var.y <- unlist(lapply(c(l1,l2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- unlist(lapply(list(l.eta0,l.eta1,l.eta2),function(x) t(residuals(x))%*%residuals(x)/N))
##  sd.y <- unlist(c(lapply(c(l1,l2), function(x) summary(x)$sigma)))  
##  sd.eta <- unlist(c(lapply(list(l.eta0,l.eta1,l.eta2), function(x) summary(x)$sigma)))
##   sd.eta <- sd.eta*c(gamma1*b1,b1,b2)
##  sd. <- c(sd.eta, sd.y)
  ##  Sigma <- diag(sd.^2)
  ##  var.eta <- var.eta*(c(gamma1*b1,b1,b2)^2)
  Sigma <- diag(c(var.eta,var.y))
##  cat("ooost\n")
  
  ## cdata2 <- cdata
  ## cdata2$eta1 <- cdata$eta1*b1
  ## lm(y1 ~ eta1, cdata2)
  ## cdata2$eta0 <- cdata$eta0*gamma1*b1
  ## lm(eta1 ~ -1 + eta0, cdata2)
  ## cdata2$eta2 <- cdata$eta2*b2
  ## lm(y4 ~ eta2, cdata2)
  ## lm(eta2 ~ eta0-1, cdata2)
  ## diag(Sigma)[1:3]
  ## sd.eta^2*c(gamma1^2*b1^2,b1^2,b2^2)
  
  newmodelpar <- modelpar
  newmodelpar$theta <- theta; newmodelpar$theta.var <- as.double(Sigma)
  return(newmodelpar)
}

Mstep_lin1b <- function(cdata,modelpar) {
  nvar <- modelpar$nvar
  npred <- modelpar$npred
  ny1 <- modelpar$nvar1
  ny <- nvar-3
  ny2 <- ny-ny1;
  
  f <- paste("eta0 ~ 1")
  if (npred>0) {
    f <- paste(f, paste("x",1:npred, sep="", collapse="+"), sep="+")
  }
  l.eta0 <- lm(f, cdata)
  l.eta1 <- lm(eta1 ~ 1 + eta0, cdata)
  l.eta2 <- lm(eta2 ~ 1 + eta0, cdata)
  
  cdata$eta00 <- cdata$eta0 - coef(l.eta0)[1]
  cdata$eta10 <- cdata$eta1 - coef(l.eta1)[1]
  cdata$eta20 <- cdata$eta2 - coef(l.eta2)[1]

  gamma0 <- coef(l.eta0)[-1]
  gamma1 <- coef(l.eta1)[-1]
  gamma2 <- coef(l.eta2)[-1]


  l1 <- list()
  for (i in 1:ny1) {
    f <- paste("y",i," ~ eta10",sep="")
    l1 <- c(l1, list(lm(f, cdata)))
  }
  l2 <- list()
  for (i in (1:ny2 + ny1)) {
    f <- paste("y",i," ~ eta20",sep="")
    l2 <- c(l2, list(lm(f, cdata)))    
  }
  c1 <- matrix(unlist(lapply(l1,coef)),ncol=2,byrow=TRUE)
  c2 <- matrix(unlist(lapply(l2,coef)),ncol=2,byrow=TRUE)

  
  mu <- c(c1[,1], c2[,1])

  ## Adjusting for fixed parameter
  b1 <- c1[1,2] 
  b2 <- c2[1,2]
  theta.y <- c(c1[-1,2]/b1,c2[-1,2]/b2)
##  theta.y <- c(c1[-1,2],c2[-1,2])
  theta.eta <- c(gamma0*b1*gamma1, gamma2*b2/(b1*gamma1))
##  theta.eta <- c(gamma0, gamma2)
  theta <- c(mu, theta.y, theta.eta)
  ## Variance parameters

##  cat(".......\n")
  N <- nrow(cdata)
  var.y <- unlist(lapply(c(l1,l2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- unlist(lapply(list(l.eta0,l.eta1,l.eta2),function(x) t(residuals(x))%*%residuals(x)/N))
##  sd.y <- unlist(c(lapply(c(l1,l2), function(x) summary(x)$sigma)))  
##  sd.eta <- unlist(c(lapply(list(l.eta0,l.eta1,l.eta2), function(x) summary(x)$sigma)))
##   sd.eta <- sd.eta*c(gamma1*b1,b1,b2)
##  sd. <- c(sd.eta, sd.y)
  ##  Sigma <- diag(sd.^2)
  var.eta <- var.eta*(c(gamma1*b1,b1,b2)^2)
  Sigma <- diag(c(var.eta,var.y))
##  cat("ooost\n")
  
  ## cdata2 <- cdata
  ## cdata2$eta1 <- cdata$eta1*b1
  ## lm(y1 ~ eta1, cdata2)
  ## cdata2$eta0 <- cdata$eta0*gamma1*b1
  ## lm(eta1 ~ -1 + eta0, cdata2)
  ## cdata2$eta2 <- cdata$eta2*b2
  ## lm(y4 ~ eta2, cdata2)
  ## lm(eta2 ~ eta0-1, cdata2)
  ## diag(Sigma)[1:3]
  ## sd.eta^2*c(gamma1^2*b1^2,b1^2,b2^2)
  
  newmodelpar <- modelpar
  newmodelpar$theta <- theta; newmodelpar$theta.var <- as.double(Sigma)
  return(newmodelpar)
}
Mstep_lin1c <- function(cdata,modelpar) {
  nvar <- modelpar$nvar
  npred <- modelpar$npred
  ny1 <- modelpar$nvar1
  ny <- nvar-3
  ny2 <- ny-ny1;
  
  f <- paste("eta0 ~ -1")
  if (npred>0) {
    f <- paste(f, paste("x",1:npred, sep="", collapse="+"), sep="+")
  }
  l.eta0 <- lm(f, cdata)
  l.eta1 <- lm(eta1 ~ -1 + eta0, cdata)
  l.eta2 <- lm(eta2 ~ -1 + eta0, cdata)
  
##  cdata$eta00 <- cdata$eta0 - coef(l.eta0)[1]
##  cdata$eta10 <- cdata$eta1 - coef(l.eta1)[1]
##  cdata$eta20 <- cdata$eta2 - coef(l.eta2)[1]

  gamma0 <- coef(l.eta0)[]
  gamma1 <- coef(l.eta1)[]
  gamma2 <- coef(l.eta2)[]


  l1 <- list()
  for (i in 1:ny1) {
    f <- paste("y",i," ~ eta1",sep="")
    l1 <- c(l1, list(lm(f, cdata)))
  }
  l2 <- list()
  for (i in (1:ny2 + ny1)) {
    f <- paste("y",i," ~ eta2",sep="")
    l2 <- c(l2, list(lm(f, cdata)))    
  }
  c1 <- matrix(unlist(lapply(l1,coef)),ncol=2,byrow=TRUE)
  c2 <- matrix(unlist(lapply(l2,coef)),ncol=2,byrow=TRUE)

  
  mu <- c(c1[,1], c2[,1])

  ## Adjusting for fixed parameter
  b1 <- c1[1,2] 
  b2 <- c2[1,2]
  theta.y <- c(c1[-1,2]/b1,c2[-1,2]/b2)
##  theta.y <- c(c1[-1,2],c2[-1,2])
  theta.eta <- c(gamma0*b1*gamma1, gamma2*b2/(b1*gamma1))
##  theta.eta <- c(gamma0, gamma2)
  theta <- c(mu, theta.y, theta.eta)
  ## Variance parameters

##  cat(".......\n")
  N <- nrow(cdata)
  var.y <- unlist(lapply(c(l1,l2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- unlist(lapply(list(l.eta0,l.eta1,l.eta2),function(x) t(residuals(x))%*%residuals(x)/N))
##  sd.y <- unlist(c(lapply(c(l1,l2), function(x) summary(x)$sigma)))  
##  sd.eta <- unlist(c(lapply(list(l.eta0,l.eta1,l.eta2), function(x) summary(x)$sigma)))
##   sd.eta <- sd.eta*c(gamma1*b1,b1,b2)
##  sd. <- c(sd.eta, sd.y)
  ##  Sigma <- diag(sd.^2)
  var.eta <- var.eta*(c(gamma1*b1,b1,b2)^2)
  Sigma <- diag(c(var.eta,var.y))
##  cat("ooost\n")
  
  ## cdata2 <- cdata
  ## cdata2$eta1 <- cdata$eta1*b1
  ## lm(y1 ~ eta1, cdata2)
  ## cdata2$eta0 <- cdata$eta0*gamma1*b1
  ## lm(eta1 ~ -1 + eta0, cdata2)
  ## cdata2$eta2 <- cdata$eta2*b2
  ## lm(y4 ~ eta2, cdata2)
  ## lm(eta2 ~ eta0-1, cdata2)
  ## diag(Sigma)[1:3]
  ## sd.eta^2*c(gamma1^2*b1^2,b1^2,b2^2)
  
  newmodelpar <- modelpar
  newmodelpar$theta <- theta; newmodelpar$theta.var <- as.double(Sigma)
  return(newmodelpar)
}


############################

Mstep_lin0 <- function(cdata,modelpar) {
  nvar <- modelpar$nvar
  npred <- modelpar$npred
  ny1 <- modelpar$nvar1
  ny <- nvar-3
  ny2 <- ny-ny1;
  
  require(structural)
  m <- sem()
  regression(m,"y1") <- "eta1$1"
  for (i in 2:ny1) 
    regression(m, paste("y",i,sep="")) <- "eta1"
  regression(m,paste("y",ny1+1,sep="")) <- "eta1$1"
  for (i in (2:ny2 + ny1))
    regression(m, paste("y",i,sep="")) <- "eta2"
  regression(m, "eta1") <- "eta0$1"
  regression(m, "eta2") <- "eta0"
  regression(m, "eta0") <- c("x1","x2")

##  e <- estimate(m,dd)
  
  return(list(par=theta,var=Sigma))
}

