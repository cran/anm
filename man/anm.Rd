\name{anm}
\alias{anm}
\title{The Analog method}
\description{
\emph{anm} is used to compute the analog method.
}
\usage{
anm(formula,data,weights=NULL,cross.valid=NULL)
}

\arguments{
    \item{formula}{a symbolic description of the model to be fit.}
    \item{data}{the data.frame containing the variables in the model.}
    \item{weights}{an optional matrix of weights to be used in the fitting process.}
    \item{cross.valid}{an optional matrix of booleans. If not specified, a cross validation is used in the fitting process.}
          }
\details{Models for \emph{anm} are specified symbolically. A typical model has the form \emph{predictand ~ terms} where \emph{terms} is a series of predictors whose specification can be of the for \emph{first + second}.
\link{anm} calls the lower level function \link[clim.pact]{anm.fit}.
}
\value{An object of class "anm".
An object of class "anm" is a list containing the following components: 
\item{coefficients}{a vector containing the values for the principal components corresponding to the maximum among observations.}
\item{contrasts}{(not used).}
\item{call}{the matched call.}
\item{terms}{the terms object used.}
\item{model}{the model frame used.}
\item{x}{the matrix used for predictors.}
\item{y}{the predictand.}
\item{weights}{the matrix of weights.}
\item{cross.valid}{equals to True if the cross.validation will be used for the fitting process.}
\item{data}{the input data.frame.}
}
\references{URL \url{http://www.R-project.org/}}
\author{Alexandra Imbert}
\seealso{link{predict.anm}, link{stepANM}}
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
anm(y ~ X1 + X2,data=test.data)
}
\keyword{models}




