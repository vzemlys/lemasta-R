##Main file
##Dependencies:
##Data: data/all_data_all_data.csv,
##      data/data_sa_Feuille1.csv,
##      data/LEMAMmod_lemmdata.csv,
##      data/LEMAMmod_aprasymas.csv
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


ftry <- eqforecast(start=c(2008,1),end=c(2008,4),eqR,endoexo,data=ladt,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1))

endol <- ladt[,colnames(ftry)]


gofd <- cbind(endol,ftry)
colnames(gofd) <- c(colnames(endol),paste("f_",colnames(ftry),sep=""))

no <- dim(endol)[2]

xyplot(gofd,screens=rep(1:no,2),col=c(rep("blue",no),rep("red",no)))

nmi <- read.csv("data/LEMAMmod_aprasymas.csv")

nmi[,1] <- tolower(as.character(nmi[,1]))

nmi <- data.frame(name=nmi[,1],saname=paste(nmi[,1],"_sa",sep=""),nicename=nmi[,2])

nmi <- nmi[order(nmi[,1]),]

res <- t(ftry)
rownames(res) <- NULL

res <- data.frame(as.character(nmi[match(colnames(endol),as.character(nmi[,2])),3]),res)


names(res) <- c("Rodiklis",paste("2008 K",1:4,sep=""))

require(xtable)


print(xtable(res),type="html",include.rows=FALSE,file="ftable.html")
