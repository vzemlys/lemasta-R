##Main file
##Dependencies:
##Data: data/all_data_all_data.csv,
##      data/data_sa_Feuille1.csv,
##      data/LEMAMmod_lemmdata.csv,
##      eviews/endoexo.csv
##Code: code.R,
##      eviews/eviewseq.R
##      


rm(list=ls())
require(tseries)

ad <- read.csv("data/all_data_all_data.csv")
sad <-read.csv("data/data_sa_Feuille1.csv")
lad <- read.csv("data/LEMAMmod_lemmdata.csv")

adt <- ts(ad[,-1],start=c(1995,1),end=c(2009,4),freq=4)
sadt <- ts(sad[,-1],start=c(1995,1),end=c(2009,4),freq=4)

ladt <- ts(lad[, -1],start=c(1995,1),end=c(2009,1),freq=4)

adt <- window(adt,start=c(1995,1),end=c(2009,1))
sadt <- window(sadt,start=c(1995,1),end=c(2009,1))

  
colnames(sadt) <- tolower(colnames(sadt))
colnames(adt)  <- tolower(colnames(adt))
colnames(ladt) <- tolower(colnames(ladt))
  
source("eviews/eviewseq.R")
source("code.R")

eqstr <- gsub("\n","",eqstr)

eqstrm <- gsub("=","-",eqstr)


eqR <- lapply(eqstrm,function(l)eviewstoR(l,varnames=colnames(adt)))

endoexo <- read.csv("eviews/endoexo.csv")

endoexo <- data.frame(name=tolower(as.character(endoexo[,1])),exo=endoexo[,2])

lower <- rep( -Inf,26)
upper <- rep(Inf,26)

upper[19] <- log(100)


ftry <- eqforecast(start=c(2008,4),end=c(2008,4),eqR,endoexo,data=ladt,lower=lower,upper=upper)

endol <- ladt[,colnames(ftry)]

colnames(ftry) <- paste("f_",colnames(ftry))


gofd <- cbind(endol,ftry)
colnames(gofd) <- c(colnames(endol),colnames(ftry))

no <- dim(endol)[2]

xyplot(gofd,screens=rep(1:no,2),col=c(rep("blue",no),rep("red",no)))



