##Main file
##Dependencies:
##Data: data/all_data_all_data.csv,
##      data/data_sa_Feuille1.csv,
##      data/LEMAMmod_lemmdata.csv,
##      data/LEMAMmod_aprasymas.csv
##      eviews/endoexo.csv
##      eviews/LEMAM.txt
##Code: code.R,
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
  

source("code.R")

eqstr <- read.eviews("eviews/LEMAM.txt")

eqstrm <- sub("(=)(.*)","-(\\2)",eqstr)


eqR <- lapply(eqstrm,function(l)eviewstoR(l,varnames=colnames(ladt)))

ee <- read.csv("data/LEMAMmod_endoexo.csv")

ee <- data.frame(name=tolower(as.character(ee[,1])),exo=ee[,2],FT=ee[,3])

endl <- end(ladt)

fch <- length(ts(1,start=endl,end=c(2012,4),freq=4))-1

future <- ts(matrix(NA,ncol=ncol(ladt),nrow=fch),start=end(ladt)+c(0,1),end=c(2012,4),freq=4)


ladt <- ts(rbind(ladt,future),start=start(ladt),end=end(future),freq=4)

require(plyr)
require(forecast)
exnm <- as.character(ee$name[ee$exo=="Exog"])

ff <- list()
for(nm in exnm) {
    print(nm)
    x <- na.omit(window(ladt[,nm],end=endl))
    fch <- dim(qpadd(end(x),end(ladt)))[1]-1
    fc <- try(forecast(auto.arima(x),h=fch))
    ff <- c(ff,list(fc))
    window(ladt[,nm],start=end(x)+c(0,1)) <- fc$mean
}

#exf <- sapply(ff,function(x)x$mean)
#colnames(exf) <- exnm

#exff <- simplefc(window(ladt[,exnm],end=endl),ee[ee$exo=="Exog",],fch)


#for(nm in exnm) {
#    window(ladt[,nm],start=endl+c(0,1)) <- exf[,nm]
#}

lower <- rep(-Inf,26)
upper <- rep(Inf,26)

upper[19] <- log(100)


ftry <- eqforecast(start=c(2009,1),end=c(2012,4),eqR,ee,data=ladt,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1))

endol <- ladt[,colnames(ftry$res)]


gofd <- cbind(endol,ftry$res)
colnames(gofd) <- c(colnames(endol),paste("f_",colnames(ftry$res),sep=""))

no <- dim(endol)[2]

xyplot(gofd,screens=rep(1:no,2),col=c(rep("blue",no),rep("red",no)))

xyplot(ftry$data[,colnames(endol)])


#nmi <- read.csv("data/LEMAMmod_aprasymas.csv")

#nmi[,1] <- tolower(as.character(nmi[,1]))

#nmi <- data.frame(name=nmi[,1],saname=paste(nmi[,1],"_sa",sep=""),nicename=nmi[,2])

#nmi <- nmi[order(nmi[,1]),]

#res <- t(ftry)
#rownames(res) <- NULL

#res <- data.frame(as.character(nmi[match(colnames(endol),as.character(nmi[,2])),3]),res)


#names(res) <- c("Rodiklis",paste("2008 K",1:4,sep=""))

require(xtable)
require(reshape)
res <- produce.tb1(ftry$data)

write.table(res,file="data.csv",row.names=FALSE,sep="\t")

butt <- paste("<input type='submit' value='Vaizduoti' id='s",0:(nrow(res)-1),"'/>",sep="")

hres <- data.frame(res,butt)
colnames(hres) <- c("Rodiklis",2008:2012,"")

print(xtable(hres),type="html",include.rows=FALSE,file="ftable.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

#aldt <- adt

#ialnm <- intersect(colnames(adt),colnames(ladt))
#sdalnm <- setdiff(colnames(ladt),colnames(adt))

#aldt[,ialnm] <- ladt[,ialnm]

#alnm <- colnames(aldt)

#aldt <- cbind(aldt,ladt[,sdalnm])

#colnames(aldt) <- c(alnm,sdalnm)

#save(eqR,endoexo,aldt,lower,upper,nmi,ftry,file="lemasta.RData")

save(eqR,ee,adt,lower,upper,ftry,file="lemasta.RData")
