rm(list=ls())

library(tseries)
library(lattice)
library(nleqslv)
library(ggplot2)
library(foreach)
library(plyr)
library(reshape)

source("10code.R")

#################################################################
#Read data

lad <- read.csv("data/lemam_data.csv")
ladt <- ts(lad[,-1],start=c(1995,1),freq=4)
colnames(ladt)  <- tolower(colnames(ladt))

###Paskutinis stulpelis kartojasi du kartus, ismetu ji.
ladt <- ladt[,-ncol(ladt)]



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


