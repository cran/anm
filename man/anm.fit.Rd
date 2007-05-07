\name{anmFit}
\alias{anmFit}
\title{Support function for anm}
\description{Basic computing engine called by \link{anm} to implement the analog method. This should usually not be used directly. 
}
\usage{
anmFit(x, y, tol = 1e-07, ...)
}
\arguments{
    \item{x}{design matrix of dimension n * p.}
    \item{y}{vector of observations of length n.}
    \item{tol}{if equal to 1, information is printed during the running of step.}
    \item{...}{currently disregarded.}
          }
\value{
 A list with components
 \item{coefficients}{vector containing the highest value among observations and the values of the predictors at this date.}
 \item{residuals}{\emph{n} vector.}
 \item{fitted.values}{\emph{n} vector.}
 \item{effects}{\emph{n} vector.}
 \item{rank}{integer, giving the rank.}
 \item{df.residual}{degrees of freedom of residuals.}
 \item{qr}{the QR decomposition, see \link{qr}.}
}
\references{URL \url{http://www.R-project.org/}}
\author{Alexandra Imbert}
\seealso{\link{anm}, \link{predictAnm}}
\keyword{models}



