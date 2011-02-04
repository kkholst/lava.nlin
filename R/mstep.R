
Mstep_nsem1 <- function(cdata,modelpar,...) {
  ny1 <- modelpar$ny1
  ny2 <- modelpar$ny2;  
  ny <- ny1+ny2
  nlatent <- modelpar$nlatent
  nvar <- ny+nlatent
  npred <- modelpar$npred
  mynames <- c(paste("y",1:ny,sep=""),paste("x",1:npred,sep=""),paste("eta",0:2,sep=""))
  cdata <- data.frame(cdata)
  names(cdata) <- mynames
  

  f <- paste("eta0 ~ 1")
  if (npred>0) {
    f <- paste(f, paste("x",1:npred, sep="", collapse="+"), sep="+")
  }
  l.eta0 <- lm(f, cdata)
  l.eta1 <- lm(eta1 ~ 1 + eta0, cdata)
  l.eta2 <- lm(eta2 ~ 1 + eta0 + I(eta0^2), cdata)

  ##  browser()
  
  
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
  theta.eta <- c(gamma0*b1*gamma1, gamma2*b2/(b1*gamma1))
  theta <- c(mu, theta.y, theta.eta)
  ## Variance parameters
##  sd.y <- unlist(c(lapply(c(l1,l2), function(x) summary(x)$sigma)))  
##  sd.eta <- unlist(c(lapply(list(l.eta0,l.eta1,l.eta2), function(x) summary(x)$sigma)))
    ##  sd.eta <- sd.eta*c(gamma1*b1,b1,b2)
##  sd. <- c(sd.eta, sd.y)
##  Sigma <- diag(sd.^2)

  N <- nrow(cdata)
  var.y <- unlist(lapply(c(l1,l2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- unlist(lapply(list(l.eta0,l.eta1,l.eta2),function(x) t(residuals(x))%*%residuals(x)/N))
  var.eta <- var.eta*(c(gamma1*b1,b1,b2)^2)
  Sigma <- diag(c(var.y,var.eta))
  res <- list(theta=c(theta,diag(Sigma)),theta.var=as.double(Sigma)) 
  return(res)  
}
