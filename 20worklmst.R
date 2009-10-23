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
ladt <- ts(lad[,-1],start=c(1995,1),end=c(2009,4),freq=4)


adt <- window(adt,start=c(1995,1),end=c(2009,1))
sadt <- window(sadt,start=c(1995,1),end=c(2009,1))

  
colnames(sadt) <- tolower(colnames(sadt))
colnames(adt)  <- tolower(colnames(adt))
colnames(ladt)  <- tolower(colnames(ladt))
  

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


require(xtable)
require(reshape)

res <- produce.tb1(ftry$data)



sc2 <- ladt

window(sc2[,"i_l"],start=c(2009,1)) <- 1.5*window(sc2[,"i_l"],start=c(2009,1))

ftry2 <- eqforecast(start=c(2009,1),end=c(2012,4),eqR,ee,data=sc2,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1))

res2 <- produce.tb1(ftry2$data)

sc3 <- ladt

window(sc3[,"i_l"],start=c(2009,1)) <- 0.5*window(sc3[,"i_l"],start=c(2009,1))

ftry3 <- eqforecast(start=c(2009,1),end=c(2012,4),eqR,ee,data=sc3,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1,maxit=50,factr=1e-3,pgtol=1e-8))

res3 <- produce.tb1(ftry3$data)

csvhtpair(res,1,"varno")
csvhtpair(res2,2,"varno")
csvhtpair(res3,3,"varno")


dd <- t(ladt[,ee$name[ee$exo=="Exog"]])[,11:14]

dd <- data.frame(rownames(dd),dd)

dd <- dd[-c(1,2,11),]

names(dd) <- c("Rodiklis",2009:2012)

a1 <- produce.form(dd,pref="scen1")

dd[1,3] <- 20
a2 <- produce.form(dd,pref="scen2")

dd[1,3] <- 30
a3 <- produce.form(dd,pref="scen3")

print(xtable(a1),type="html",include.rows=FALSE,file="form1.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

print(xtable(a2),type="html",include.rows=FALSE,file="form2.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

print(xtable(a3),type="html",include.rows=FALSE,file="form3.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

save(eqR,ee,ladt,lower,upper,ftry,file="lemasta.RData")
