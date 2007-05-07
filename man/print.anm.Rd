\name{printAnm}
\alias{printAnm}
\title{Print some components of an anm object.}
\description{Prints the coefficients of an \link{anm} object.
}
\usage{
printAnm(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
    \item{x}{the \link{anm} object.}
    \item{digits}{the vector defining the format of printing.}
    \item{...}{currently disregarded.}
          }
\references{URL \url{http://www.R-project.org/}}
\author{Alexandra Imbert}
\seealso{\link{anm}, \link{predictAnm}}
\examples{
library(survival)
library(clim.pact)
data(susendal)
data(temp.era)
y<-susendal$V6 # temperatures
X<- eof$PC[,c(1,2)]
calibration <- c(susendal$V4>1979 & susendal$V4<1990 & (susendal$V3==1 | susendal$V3==2 | susendal$V3==12))
evaluation <- c((susendal$V4>1990 & susendal$V4<1993 | susendal$V4==1990) & (susendal$V3==1 | susendal$V3==2 | susendal$V3==12))
y.calib <- y[calibration]
y.eval <- y[evaluation]
eof.calib <- c(eof$yy>1979 & eof$yy<1990)
eof.eval <- c(eof$yy> 1990 & eof$yy<1993| eof$yy==1990)
period <- c(calibration, evaluation)
y.period <- y[(susendal$V4>1979 & susendal$V4<1993) & (susendal$V3==1 | susendal$V3==2 | susendal$V3==12)]
test.data <-data.frame(y=y.period,
                       X1=X[eof$yy< 1993 & eof$yy> 1979,1],
                       X2=X[eof$yy< 1993 & eof$yy> 1979,2],
                       yy=eof$yy[eof.calib | eof.eval],
                       mm=eof$mm[eof.calib | eof.eval],
                       dd=eof$dd[eof.calib | eof.eval])
test.anm<-anm(y ~ X1 + X2,data=test.data)
printAnm(test.anm)
}
\keyword{models}
