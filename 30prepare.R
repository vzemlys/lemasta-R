rm(list=ls())

library(tseries)
library(lattice)
library(nleqslv)
library(ggplot2)
library(foreach)
library(plyr)
library(reshape)
library(xtable)

source("10code.R")

#################################################################
#Read data

lad <- read.csv("data/VK0412_data.csv")
ladt <- ts(lad[,-1],start=c(1995,1),freq=4)
colnames(ladt)  <- tolower(colnames(ladt))



################################################################
###Convert eviews formulas to R
eqstr <- read.eviews("eviews/lygtys_01q1-09q4.txt")

eqstrm <- sub("(=)(.*)","-(\\2)",eqstr)

eqR <- lapply(eqstrm,function(l){print(l);eviewstoR(l,varnames=colnames(ladt))})

eqstr2 <- read.eviews("eviews/lygtys_95q1-09q4.txt")
eqstrm2 <- sub("(=)(.*)","-(\\2)",eqstr2)
eqR2 <- lapply(eqstrm2,function(l){print(l);eviewstoR(l,varnames=colnames(ladt))})

eqstr3 <- read.eviews("eviews/lygtys_97q1-09q4.txt")
eqstrm3 <- sub("(=)(.*)","-(\\2)",eqstr3)
eqR3 <- lapply(eqstrm3,function(l){print(l);eviewstoR(l,varnames=colnames(ladt))})

eqstr4 <- read.eviews("eviews/lygtys_99q1-09q4.txt")
eqstrm4 <- sub("(=)(.*)","-(\\2)",eqstr4)
eqR4 <- lapply(eqstrm3,function(l){print(l);eviewstoR(l,varnames=colnames(ladt))})


# "lygtys_97q1-09q4.txt" "lygtys_99q1-09q4.txt"
################################################################
##Insert exogenous forecasts

exof <- read.csv("data/VK0412_scen.csv")
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

##Last three of exof can be endogenous
ene <- c("x_r_sa"   ,   "bp_im_sa" ,   "g_n_sa" )

ee <- rbind(data.frame(name=colnames(exof),exo="Exog"),
            data.frame(name=setdiff(colnames(ladt),colnames(exof)),exo="Endog")
            )

ee$exo[ee$name%in% ene] <- "Endog"

