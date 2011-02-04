fastapprox <- function(a,t,z) {
  arglist <- list(name="FastApprox",
                  a=a,
                  t=t,
                  z=z,
                  DUP=FALSE,PACKAGE="lava.nlin")
  res <- do.call(".Call",arglist)
  return(res)
}

