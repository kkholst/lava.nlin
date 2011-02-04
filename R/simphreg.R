sim.phwreg <- function(X,lambda=0.8,p=2,beta=rep(1,NCOL(X)),
                        t=seq(0,10,length.out=NROW(X)),breakties=0,cens=rexp(NROW(X))) {
  a0 <- function(t) lambda*p*(lambda*t)^(p-1)
  A0 <- function(t) (lambda*t)^p
  A0i <- function(eta) eta^(1/p)/lambda
  A <- A0(t)
  n <- NROW(X)
  U <- rexp(n, 1) #give everyone a random death time, on the CH scale
  Z <- U*exp(-X%*%beta)
  T <- A0i(Z)
  if (breakties!=0)
    T <- T+runif(n,0,breakties)  
  Delta <- (T<cens)
  T[!Delta] <- cens[!Delta]
  data.frame(t=T,status=Delta*1)
}
