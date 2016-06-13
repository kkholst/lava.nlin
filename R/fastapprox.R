#' Fast approximation
#' 
#' Interpolation function
#' 
#' 
#' @param x1 Vector or matrix of observations
#' @param x2 Observation times
#' @param y Response paired with x1 (can optionally be given as the second
#' column in x1)
#' @return New list of responses y(x2) interpolated (ceiling) from x1
#' @author Klaus K. Holst
#' @keywords utilities
#' @examples
#' 
#' fastapprox(seq(-10,10,length.out=10),seq(-10,10,length.out=100),1:10)
#' 
#' @export fastapprox
fastapprox <- function(x1,x2,y) {
  if (is.matrix(x1)) {
    y <- x1[,2]; x1 <- x1[,1]
  }    
  arglist <- list("FastApprox",
                  a=x1,
                  t=y,
                  z=x2,
                  DUP=FALSE,PACKAGE="lava.nlin")
  res <- do.call(".Call",arglist)
  return(res)
}

