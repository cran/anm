.packageName <- "anm"
# Code to implement the analog method
# It is presented in the typical way of defining models in R.

# This routine gets the formula of the model to be fit: eg formula= "y ~ X"
# and a data frame containing the variables in the model.
# It returns an object of class "anm".

anm <- function(formula,data,weights=NULL,cross.valid=NULL) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$x <- mf$y <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms") # allow model.frame to update it
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf, contrasts=NULL)
  z<-anm.fit(x, y)
  date.max<- seq(1,dim(mf[1])[1],by=1)[is.element(mf[,1],max(mf[,1],na.rm=TRUE))]
  no.pc <- length(mf)
  if (!is.null(cross.valid)){
      cross.valid <- FALSE
      no.pc <- no.pc-1
  }else cross.valid <- TRUE
  if (!is.null(weights)) no.pc <- no.pc-1
  coeff <- rep(NA,length(mf))
  #coeff[1] <- max(mf[,1],na.rm=TRUE)
  coeff[1] <- 0 
  for (i.pc in 2:length(mf))  coeff[i.pc] <- mf[date.max,i.pc]
  z <- list(coefficients=coeff,contrasts=attr(x, "contrasts"),call=cl,
            terms= mt,model=mf,x=mf[2:no.pc],y=mf[1],
            weights=weights,cross.valid=cross.valid,data=data)
  class(z) <- "anm"
  invisible(z)
}

############################################
anm.fit <- function (x, y, tol = 1e-07, ...)
{
    if (is.null(n <- nrow(x))) stop("`x' must be a matrix")
    if(n == 0) stop("0 (non-NA) cases")
    p <- ncol(x)
    if (p == 0) {
        ## oops, null model
        cc <- match.call()
        cc[[1]] <- as.name("anm.fit.null")
        return(eval(cc, parent.frame()))
    }
    ny <- NCOL(y)
    ## treat one-col matrix as vector
    if(is.matrix(y) && ny == 1)
        y <- drop(y)
    if (NROW(y) != n)
	stop("incompatible dimensions")
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    z <- .Fortran("dqrls",
		  qr = x, n = n, p = p,
		  y = y, ny = ny,
                  tol = as.double(tol),
		  coefficients = mat.or.vec(p, ny),
		  residuals = y, effects = y, rank = integer(1),
		  pivot = 1:p, qraux = double(p), work = double(2*p),
                  PACKAGE="base")
    coef <- z$coefficients
    pivot <- z$pivot
    r1 <- 1:z$rank
    dn <- colnames(x); if(is.null(dn)) dn <- paste("x", 1:p, sep="")
    nmeffects <- c(dn[pivot[r1]], rep("", n - z$rank))
    if (is.matrix(y)) {
	coef[-r1, ] <- NA
	coef[pivot, ] <- coef
	dimnames(coef) <- list(dn, colnames(y))
	dimnames(z$effects) <- list(nmeffects,colnames(y))
    } else {
	coef[-r1] <- NA
	coef[pivot] <- coef
	names(coef) <- dn
	names(z$effects) <- nmeffects
    }
    z$coefficients <- coef
    r1 <- y - z$residuals ; 
    c(z[c("coefficients", "residuals", "effects", "rank")],
      list(fitted.values = r1, assign = attr(x, "assign"),
	   qr = z[c("qr", "qraux", "pivot", "tol", "rank")],
	   df.residual = n - z$rank))
}

####################################################
print.anm <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    invisible(x)
}

#####################################################
# Predicted values based on the model object. 
# "object" is the "anm" object inheriting from "anm" routine.

predict.anm<- function(object,newdata=NULL,se.fit=FALSE,...) {
nx <- length(object$x) # dimension of the problem
if (is.null(newdata)) {
   nt.1 <- length(object$y[,1])
   nt.2 <- nt.1
   object.y <- object$y[,1]
   newdata.y <- object.y
   mat.eval <- matrix(NA, ncol=nx, nrow=nt.2)
   mat.calib <- matrix(NA, ncol=nx, nrow=nt.1)
   for (i.npc in 1:nx){
      mat.calib[,i.npc] <- object$x[1:nt.1,i.npc]
      mat.eval[,i.npc] <- mat.calib[,i.npc]  
      }
} else {
   nt.1 <- length(object$y[,1])
   nt.2 <- length(newdata$X1)
   object.y <- object$y[,1]
   newdata.y <- object.y
   mat.eval <- matrix(NA, ncol=nx, nrow=nt.2)
   mat.calib <- matrix(NA, ncol=nx, nrow=nt.1)
   for (i.npc in 1:nx){
     mat.eval[,i.npc] <- eval(parse(text=paste("newdata$X",i.npc,sep="")))
     mat.calib[,i.npc] <- object$x[,i.npc]
     }
   }
date.min<-rep(NA,nt.2)
dist <- rep(NA,nt.2)
analog<-rep(NA,nt.2)
error<-rep(NA,nt.2)
same <- rep(FALSE,nt.1)
anlg.eof <- rep(NA,nx)
ref.eof <- rep(NA,nx)
cor.eof <- rep(NA,nt.2)
cpt <- 0
i.max <- 1
maxi.anlg <- rep(NA, floor((nt.2)/30))
mini.anlg <- rep(NA, floor((nt.2)/30))
# each day is considered iteratively
for (iit in 1:nt.2) {
  if (object$cross.valid) {
     same <- (object.y == newdata.y[iit])
      for (ix in 1:nx) {
       same <- same & (mat.calib[,ix] == mat.eval[iit,ix]) 
     }
  }
  d <- rep(0,nt.1)  # distance in phase space
  if (is.null(object$weights)){
  for (ix in 1:nx) {
      d <- d + (mat.calib[,ix]-mat.eval[iit,ix])^2
    }
  }
  if (!is.null(object$weights)){
    for (ix in 1:nx) {
      d <- d + object$weights[1,ix]*(mat.calib[,ix]-mat.eval[iit,ix])^2
     }
  }
  if ((object$cross.valid) & (sum(same)>0)) d[same] <- NA
  if (length(seq(1,length(d),by=1)[is.element(d,min(d,na.rm=TRUE))])>1)
  date.min[iit]<-seq(1,length(d),by=1)[is.element(d,min(d,na.rm=TRUE))][1]
  else
  date.min[iit]<- seq(1,length(d),by=1)[is.element(d,min(d,na.rm=TRUE))]
  dist[iit]<- min(d,na.rm=TRUE)
  analog[iit] <- object$y[date.min[iit],1]
  error[iit] <- abs(analog[iit]-newdata.y[iit])
  cpt <- cpt+1
  if(cpt==30){
     maxi.anlg[i.max] <- max(analog[(iit-29):iit],na.rm=TRUE)
     mini.anlg[i.max] <- max(-1*(analog[(iit-29):iit]),na.rm=TRUE)
     i.max <- i.max+1
     cpt <- 0
  }
}

# verifications with observations
if (is.null(newdata)){
  correlation <- cor(analog,newdata.y)
  rmse <- (1/nt.2)*sum((analog-newdata.y)^2)
  rmse <- sqrt(rmse)
  if (se.fit) predict <-list(problem.dimension=nx,
                             period.length=nt.1,
                             d.min=dist, date.min=date.min,
                             analog=analog,
                             maxi.anlg=maxi.anlg,
                             mini.anlg=mini.anlg, 
                             error=error,
                             correlation=correlation,
                             rmse=rmse)
  else predict <- analog
 } else predict <- analog
predict
}

######################################################
# ex. of call:
# plotANM(test.anm,tmp=TRUE,"Tromsoe","eof_ERA-15_TEM_16E31E-64N73N",FALSE)

plotANM <- function(x,tmp,station,eof.file,leps){
res.predict <- predict.anm(x,se.fit=TRUE)
nt <- res.predict$period.length
newdata.y <- x$y[1:nt,1]
vect.time <- c(1:nt) 
dates <- x$data$yy + x$data$mm/12 + x$data$dd/365.25
if (tmp) var <- "temperature"
else var <- "rainfall"
plot.new()
if (!leps) par(ask=TRUE)
 if (leps) { 
   figname<- paste("distmin_",var,"_",station,"_",eof.file,".eps",sep="")
   postscript(file = figname,onefile=TRUE,horizontal=FALSE)
   }
 par(ps=16,cex.sub=0.7,cex.main=0.9)
 plot(dates,res.predict$d.min,type="l",col="darkblue",lwd=2,xlab="time", ylab="dist")
 title(main = "Minimum distance for predictions",font.main = 4)
 text(quantile(dates,0.01),
     max(res.predict$d.min),pos=4,cex=0.6,
     paste("Daily winter",var,"in",station,"- EOF=",eof.file, sep=" "))
 grid()
 if (leps) { 
    dev.off()
    figname<- paste("obspred_",station,"_",var,"_",eof.file,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE)
   } 
 par(ps=16,cex.sub=0.7,cex.main=0.9)
 plot(dates,newdata.y,type="l",col="darkred",lwd=2,lty=1,cex=0.75,
       xlab="time",ylab=var)
 lines(dates,res.predict$analog, lty=2,col="grey40",lwd=2,cex=0.75)
 title(main = "Observations and predictions",font.main = 4)
 #mtext("(a)",line= 2,side=1,at=-30)
 legend(quantile(dates,0.01),
                 max(c(newdata.y,res.predict$analog)),
                 c("Obs.","Predictions"),cex=0.75,
                 col=c("darkred","grey40"),
                 lwd=c(2,2),lty=c(1,2),
                 merge=TRUE,bg="grey95")
 text(quantile(dates,0.01),
     min(c(newdata.y,res.predict$analog)),pos=4,cex=0.6,
     paste("Daily winter",var, "in",station,"eof=",eof.file,sep=" "))
 grid()

  if (leps) {
    dev.off()
    figname<- paste("error_",var,"_",station,"_",eof.file,".eps",sep="")
    postscript(file = figname,onefile=TRUE,horizontal=FALSE)
    }
  par(ps=16,cex.sub=0.7,cex.main=0.9)
  plot(dates,res.predict$error,type="l",lwd=2,col="darkblue",
        xlab="time",ylab=var)
  title(main ="Absolute value of Error between observations and predictions",
		font.main = 4)
 grid()
  #mtext("(b)",line= 2,side=1,at=1990)
  text(quantile(dates,0.01),
     max(res.predict$error),pos=4,cex=0.6,
     paste("Daily winter",var,"in",station,"- EOF=",eof.file,sep=" "))
  grid()

  if (leps) { 
    dev.off()
   }
}

###############################################################
# stepwise algorithm 
# anm.obj is the "anm" object returned by the "anm" routine.

stepANM <- function(anm.obj,trace=1,steps=8) {
library(xtable)
cor <- rep(NA,steps-1) # MODIF
rmse <- rep(NA,steps-1) # MODIF
incl <- "X1 + X2"
formula <- "y ~ X1 + X2"
formula <- as.formula(formula)
if (steps==2) stop("The number of steps and of predictor variables should be higher than 2.")
  if (steps > length(anm.obj$x)) steps <- length(anm.obj$x)
  # Determination of rmse and correlation at each step
  for (i.pc in 3:(steps+1)){
     test.anm<-anm(formula,data=anm.obj$data,weights=anm.obj$weights)
     result.predict <- predict.anm(test.anm,se.fit=TRUE)
     cor[i.pc-1] <- result.predict$correlation
     rmse[i.pc-1] <- result.predict$rmse
     if (trace==1){
       print(formula)
       cat("Rmse:\n",rmse[i.pc-1],"\n\n",sep="")
       cat("Correlation:\n",cor[i.pc-1],"\n\n",sep="")
     }
     incl <- paste(incl, " + X",i.pc, sep="")
     formula <- paste("y ~ ",incl,sep="")
     formula <- as.formula(formula)
   }
 # number of steps which returns the lowest rmse.
  min.rmse <- min(rmse,na.rm=TRUE) 
  step.min <- seq(1,steps,by=1)[is.element(rmse,min(rmse,na.rm=TRUE))]
  correlation <- cor[step.min]
    if (step.min == 2) model <- as.formula("y ~ X1 + X2")
    else{
       pred <- "X1 + X2"
       for (i in 3:step.min){
       pred <- paste(pred, " + X",i, sep="")
       }
       model <- paste("y ~ ",pred,sep="")
       model <- as.formula(model)
     }

   tt<- terms(model)
   tmp<-attr(tt,"term.labels")
  
# plotting of rmse and correlation.
  n.x <- c(1:steps)
  cor <- 10*cor
 #plot(c(1,steps),c(0,8),type="n",xlab="step",ylab="y")
 #lines(n.x,cor,type="l", col="darkblue")
 #lines(n.x,rmse, col="red",lwd=2,lty=2)
 #axis(4,seq(0,1,by=0.1)*10,labels=as.character(seq(0,1,by=0.1)), col="darkblue")
 #title(main = "Correlation and RMSE between analogs based on EOF and obs",
 #      font.main = 4) 
 #legend(6,3,c("cor","RMSE"), col=c("darkblue","red"),lwd=c(2,2),
 #       lty=c(2,2),bg="grey95",cex=0.8)

cat("Model:\n",deparse(model),"\n\n",sep="")
cat("Step.min:\n",step.min,"\n\n",sep="")
cat("Rmse.min:\n",min.rmse,"\n\n",sep="")
cat("Correlation:\n",correlation,"\n\n",sep="")

step.res <- list(Call=test.anm$call,PC=tmp,anm.obj=test.anm,coefficients=test.anm$coefficients, step.min = step.min,model = model,Rmse = min.rmse, correlation = correlation)
step.res 
}
