\name{predictAnm}
\alias{predictAnm}
\title{Predict method for anm objects.}
\description{
Returns the predicted values based on the \link{anm} object.
}
\usage{
predictAnm(object,newdata=NULL,se.fit=FALSE, ...)
}
\arguments{
    \item{object}{the anm object inheriting from \link{anm} routine.}
    \item{newdata}{an optional independant data. If specified, only the vector of predictions is returned.}
    \item{se.fit}{if false, only the vector of predictions is returned.}
    \item{...}{further arguments passed to or from other methods.}
          }
\value{
 A list with components
 \item{problem.dimension}{the number of predictor variables.} 
 \item{period.length}{the time period.} 
 \item{d.min}{the vector of minimum distances.}
 \item{date.min}{the vector containing the dates corresponding to the minimum distances.}
 \item{analog}{the vector of predictions.}
 \item{maxi.anlg}{monthly maxima values of predictions.}
 \item{mini.anlg}{monthly minima values of predictions.}
 \item{error}{vector of errors between predictions and observations at each date.}
 \item{correlation}{correlation coefficient between predictions and observations.}
 \item{rmse}{root mean square errors between predictions and observations.}
}
\references{URL \url{http://www.R-project.org/}}
\author{R.E. Benestad and Alexandra Imbert}
\seealso{\link{anm}, \link{stepANM}}
\examples{
library(survival)
library(clim.pact)
data(temp.era)
data(susendal)
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
res <- predictAnm(test.anm)
}
\keyword{models}


