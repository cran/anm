\name{stepANM}
\alias{stepANM}
\title{Choose a model by the analog method in a stepwise algorithm}
\description{Performs the analog method step by step to select a model and plots on the same graph both correlation and rmse at each step.
}
\usage{
stepANM(anm.obj,trace=1,steps=8)
}
\arguments{
    \item{anm.obj}{the anm object inheriting from \link{anm} routine.}
    \item{trace}{if equal to 1, information is printed during the running of the stepwise algorithm.}
    \item{steps}{maximum number of steps, forced to the number of predictor variables if \emph{steps} exceeds it.}
          }
\value{
 A list with components
 \item{Call}{the matched call.}
 \item{PC}{the predictor variables selected.}
 \item{anm.obj}{the \link{anm} object selected.}
 \item{coefficients}{the coeffecients of the \link{anm} object.}
 \item{step.min}{the number of steps which returns the minimum rmse.}
 \item{model}{the model corresponding to the minimum rmse.}
 \item{Rmse}{the minimum root mean square error.}
 \item{correlation}{the correlation between predictions and observations for the selected model.}
}
\note{
The running of the stepwise algorithm can be quite slow especially if the number of steps specified in the \emph{steps} argument is high.
}  
\author{Alexandra Imbert}
\seealso{\link{anm}, \link{predictAnm}}
\examples{
library(survival)
library(clim.pact)
data(susendal)
data(temp.era)
y<-susendal$V6 # temperatures
X<- eof$PC[,c(1,2,3)]
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
                       X3=X[eof$yy< 1993 & eof$yy> 1979,3],
                       yy=eof$yy[eof.calib | eof.eval],
                       mm=eof$mm[eof.calib | eof.eval],
                       dd=eof$dd[eof.calib | eof.eval])
test.anm<-anm(formula=y ~ X1 + X2 + X3,data=test.data)
stepANM(test.anm,steps=3)
}
\keyword{models}
