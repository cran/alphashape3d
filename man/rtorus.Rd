% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtorus.R
\name{rtorus}
\alias{rtorus}
\title{Generate points in the torus}
\usage{
rtorus(n, r, R, ct = c(0, 0, 0), rotx = NULL)
}
\arguments{
\item{n}{Number of observations.}

\item{r}{Minor radius (radius of the tube).}

\item{R}{Major radius (distance from the center of the tube to the center of
the torus).}

\item{ct}{A vector with the coordinates of the center of the torus.}

\item{rotx}{If not NULL, a rotation through an angle \code{rotx} (in
radians) about the \eqn{x}-axis is performed.}
}
\description{
This function generates \eqn{n} random points within the torus whose minor
radius is \eqn{r}, major radius is \eqn{R} and center is \eqn{ct}.
}
\examples{

T1 <- rtorus(2000, 0.5, 2.5)
bbox3d(color = c("white", "black"))
points3d(T1, col = 4)


T2 <- rtorus(2000, 0.5, 2.5, ct = c(2, 0, 0.5), rotx = pi/2)
points3d(T2, col = 2)

}
