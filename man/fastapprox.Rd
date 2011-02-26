\name{fastapprox}
\alias{fastapprox}
\title{Fast approximation}
\description{
  Interpolation function
}
\usage{
fastapprox(x1,x2,y) 
}
\arguments{
  \item{x1}{Vector or matrix of observations}
  \item{x2}{Observation times}
  \item{y}{Response paired with x1 (can optionally be given as the second column in x1)}
}
\value{
  New list of responses y(x2) interpolated (ceiling) from x1
}
\examples{
fastapprox(seq(-10,10,length.out=10),seq(-10,10,length.out=100),1:10)
}
\author{Klaus K. Holst}
\keyword{ utilities }

