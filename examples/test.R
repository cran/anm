# program to test the analog method
# source("test.R")

rm(list=ls())   # clear the memory
source("~/analog/analog.model/R/anm.R")
library(survival)
library(clim.pact)

# To test precipitations
load("~/analog/data/eof_ERA-15_TEM_16E31E-64N73N_DJF_day.Rdata")
#load("data/eof_ERA-15_SLP_5E12E-57N66N_DJF_day.Rdata")
#load("data/eof_ERA-15_TEM_16E31E-64N73N_DJF_day.Rdata")
rr<-read.table("~/alexi/scenario/obs80-99-txt/77750-obs2.ok")
y<-rr$V6
#print(summary(y))

iPCs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
X<- eof$PC[,iPCs]
calibration <- c(rr$V4>1979 & rr$V4<1990 & (rr$V3==1 | rr$V3==2 | rr$V3==12))
evaluation <- c((rr$V4>1990 & rr$V4<1993 | rr$V4==1990) & (rr$V3==1 | rr$V3==2 | rr$V3==12))
y.calib <- y[calibration]
y.eval <- y[evaluation]
eof.calib <- c(eof$yy>1979 & eof$yy<1990)
eof.eval <- c(eof$yy> 1990 & eof$yy<1993| eof$yy==1990)
period <- c(calibration, evaluation)
y.period <- y[(rr$V4>1979 & rr$V4<1993) & (rr$V3==1 | rr$V3==2 | rr$V3==12)]

test.data <-data.frame(y=y.period,
                       X1=X[eof$yy< 1993 & eof$yy> 1979,1],
                       X2=X[eof$yy< 1993 & eof$yy> 1979,2],
                       X3=X[eof$yy< 1993 & eof$yy> 1979,3],
                       X4=X[eof$yy< 1993 & eof$yy> 1979,4],
                       X5=X[eof$yy< 1993 & eof$yy> 1979,5],
                       X6=X[eof$yy< 1993 & eof$yy> 1979,6],
                       X7=X[eof$yy< 1993 & eof$yy> 1979,7],
                       X8=X[eof$yy< 1993 & eof$yy> 1979,8],
                       X9=X[eof$yy< 1993 & eof$yy> 1979,9],
                       X10=X[eof$yy< 1993 & eof$yy> 1979,10],
                       X11=X[eof$yy< 1993 & eof$yy> 1979,11],
                       X12=X[eof$yy< 1993 & eof$yy> 1979,12],
                       X13=X[eof$yy< 1993 & eof$yy> 1979,13],
                       X14=X[eof$yy< 1993 & eof$yy> 1979,14],
                       X15=X[eof$yy< 1993 & eof$yy> 1979,15],
                       X16=X[eof$yy< 1993 & eof$yy> 1979,16],
                       X17=X[eof$yy< 1993 & eof$yy> 1979,17],
                       X18=X[eof$yy< 1993 & eof$yy> 1979,18],
                       X19=X[eof$yy< 1993 & eof$yy> 1979,19],
                       X20=X[eof$yy< 1993 & eof$yy> 1979,20],
                       yy=eof$yy[eof.calib | eof.eval],
                       mm=eof$mm[eof.calib | eof.eval],
                       dd=eof$dd[eof.calib | eof.eval],
                       cal=c(rep(TRUE,length(y.calib)),rep(FALSE,length(y.eval))))

yy.length <- length(test.data$yy)
wgt <- cbind(rep(eof$W[1],yy.length),rep(eof$W[2],yy.length))
cpt <- 3
for (i.w in 3:dim(test.data)[2]){
 wgt <- cbind(wgt,rep(eof$W[cpt],yy.length))
 cpt <- cpt+1
 if (cpt==8) cpt <-1
}
crossval <- matrix(rep(T,yy.length*yy.length),yy.length,yy.length)

#test.anm<-anm(formula=y ~ X1 + X2+ X3+ X4+ X5+ X6+ X7+ X8+ X9+ X10+ X11+ X12+ X13+ X14+ X15+ X16+ X17+ X18+ X19+ X20,data=test.data,weights=NULL,cross.valid=crossval)
#test.anm<-anm(formula=y ~ X1 + X2+ X3+ X4+ X5+ X6+ X7,data=test.data,cross.valid=NULL,weights=NULL)
test.anm<-anm(formula=y ~ X1 + X2+ X3+ X4+ X5+ X6+ X7,data=test.data,weights=wgt)
#predict.anm(test.anm)
res.step <- stepANM(test.anm,steps=3)












