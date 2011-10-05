\name{rtorus}
\alias{rtorus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate points in the torus
}
\description{
This function generates \eqn{n} random points within the torus whose minor radius is \eqn{r}, major radius is \eqn{R} and center is \eqn{ct}.
}
\usage{
rtorus(n, r, R, ct = c(0, 0, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations.
}
  \item{r}{
Minor radius (radius of the tube).
}
  \item{R}{
Major radius (distance from the center of the tube to the center of the torus).
}
  \item{ct}{
A vector with the coordinates of the center of the torus.
}
}
\examples{
T1 <- rtorus(2000, 0.5, 1)
rgl.bbox()
rgl.points(T1, col = 4)
}