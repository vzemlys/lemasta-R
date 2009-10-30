rm(list=ls())

library(tseries)
library(lattice)
library(nleqslv)
library(ggplot2)

#################################################################
#Read data

lad <- read.csv("data/lemam_data.csv")
ladt <- ts(lad[,-1],start=c(1995,1),freq=4)
colnames(ladt)  <- tolower(colnames(ladt))

###Paskutinis stulpelis kartojasi du kartus, ismetu ji.
ladt <- ladt[,-ncol(ladt)]

source("10code.R")

################################################################
###Convert eviews formulas to R
eqstr <- read.eviews("eviews/LEMAM0922.txt")
eqstrm <- sub("(=)(.*)","-(\\2)",eqstr)
eqR <- lapply(eqstrm,function(l)eviewstoR(l,varnames=colnames(ladt)))

################################################################
##Insert exogenous forecasts

exof <- read.csv("data/scenarioDC.csv")
exof <- ts(exof[,-1],start=c(1995,1),freq=4)
colnames(exof)  <- tolower(colnames(exof))

###This line should be unnecessary if everything is ok
exof <- exof[,intersect(colnames(exof),colnames(ladt))]
###

###Future starts at 2009 q2
window(ladt,start=c(2009,2)) <- NA

for(nm in colnames(exof)) {
    window(ladt[,nm],start=c(2009,2)) <- window(exof[,nm],start=c(2009,2))
}

ee <- rbind(data.frame(name=colnames(exof),exo="Exog"),
            data.frame(name=setdiff(colnames(ladt),colnames(exof)),exo="Endog")
            )

################################################################
##Test the model

system.time(ftry <- eqforecast(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt,leave=TRUE,use.jacobian=TRUE))

#system.time(ftry <- eqforecast(start=c(2000,1),end=c(2008,4),eqR,ee,data=ladt,leave=FALSE,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1)))


endol <- ladt[,colnames(ftry$res)]


gofd <- cbind(endol,ftry$res)
colnames(gofd) <- c(colnames(endol),paste("f_",colnames(ftry$res),sep=""))

no <- dim(endol)[2]

xyplot(gofd,screens=rep(1:no,2),col=c(rep("blue",no),rep("red",no)))


