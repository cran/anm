\name{plotANM}
\alias{plotANM}
\title{Plot Diagnostics for an anm Object.}
\description{Three plots are provided: a plot of the minimum distances versus time, a plot comparing the analogs from \link[clim.pact]{EOF} and observations and a plot of errors versus time.}
\usage{
plotANM(x,tmp,station,eof.file,leps)
}
\arguments{
    \item{x}{the anm object inheriting from \link{anm} routine and for which prediction is desired.}
    \item{tmp}{True if the analysis is on temperature, False if on precipitation.}
    \item{station}{the name of the station.}
    \item{eof.file}{string giving the name of the eof file used for the study.}
    \item{leps}{if true, postscripts are created for the plots.}
          }
\author{Alexandra Imbert}
\seealso{\link{anm}, \link{stepANM}, \link{predict.anm}, \link{print.anm}}
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
plotANM(test.anm,TRUE,"Susendal","eof_ERA-15_TEM_16E31E-64N73N_DJF",FALSE)
}
\keyword{models}

