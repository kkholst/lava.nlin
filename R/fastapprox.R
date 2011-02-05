fastapprox <- function(x1,x2,y) {
  if (is.matrix(x1)) {
    y <- x1[,2]; x1 <- x1[,1]
  }    
  arglist <- list(name="FastApprox",
                  a=x1,
                  t=y,
                  z=x2,
                  DUP=FALSE,PACKAGE="lava.nlin")
  res <- do.call(".Call",arglist)
  return(res)
}

