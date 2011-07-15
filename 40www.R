source("30prepare.R")

################################################################
##Prepare exogenous variables (convert to yearly and back to quarterly)
exo2y <- read.csv("tables/exo2y.csv")
exo2yrest <- read.csv2("tables/exo2yrest.csv")
colnames(exo2yrest)[-2:-1] <- 2006:2011

q2y <- read.csv("tables/q2y.csv")

lydt <- q2y.meta(ladt,q2y)

exotb <-  produce.tb(lydt,exo2y)

exotb$level[,6:8] <- round(exotb$level[,6:8],2)

scy <- inverse.tb(exotb$level,exo2y)

scq <- y2q.meta(scy,exo2y)

ladt <- introduce.exo(scq,ladt,exo2y)

######################################################
##Prepare additional exogenous tables
exoadd2y <- read.csv("tables/exoadd2y.csv")
exoatb <- produce.tb(lydt,exoadd2y)
exoameta <- read.csv("tables/eqtoremove.csv")


###############################################################
###Prepare initial tables for www

##Do forecast
Rprof()
ftry <- eqforecast(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt,leave=TRUE,use.jacobian=TRUE,control=list(ftol=1e-3))
Rprof(NULL)

fp.prep <- eq.prepare(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt)
fp.prep2 <- eq.prepare(start=c(2009,1),end=c(2011,4),eqR2,ee,data=ladt)
fp.prep3 <- eq.prepare(start=c(2009,1),end=c(2011,4),eqR3,ee,data=ladt)
fp.prep4 <- eq.prepare(start=c(2009,1),end=c(2011,4),eqR4,ee,data=ladt)

metamodel <- list(list(eqlist=eqR,prep=fp.prep),
                  list(eqlist=eqR2,prep=fp.prep2),
                  list(eqlist=eqR3,prep=fp.prep3),
                  list(eqlist=eqR4,prep=fp.prep4)
                  )
Rprof()
fp <- eqprof(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt,leave=TRUE,use.jacobian=FALSE,control=list(ftol=1e-3),prep=fp.prep)
Rprof(NULL)


flydt <- q2y.meta(ftry$data,q2y)

###Read meta data
mreal <- read.csv("tables/tb1real.csv")
mnom <-  read.csv("tables/tb1nom.csv")

mprice <- read.csv("tables/tb2.csv")
mwf <- read.csv("tables/tb3.csv")

tbreal1 <- produce.tb(flydt,mreal,years=2006:2011,gdpshare=as.character(mreal[1,2]))

tbnom1 <- produce.tb(flydt,mnom,years=2006:2011,as.character(mnom[1,2]))

tb2 <- produce.tb(flydt,mprice,years=2006:2011)

tb3 <- produce.tb(flydt,mwf,years=2006:2011)

dd <- exotb$level

tb1 <- list(real=tbreal1,nominal=tbnom1)
tb2 <- list(real=tb2)
tb3 <- list(real=tb3)

exotb$level <- exotb$level[,-2]
exotb$growth <- exotb$growth[,-2]

rest <- list(table=exotb,rest=list(upper=exo2yrest[exo2yrest$Bound=="upper",-2],lower=exo2yrest[exo2yrest$Bound=="lower",-2]))
  
exoadd <- list(table=exoatb)

scen <- list(table=list(tb.conform(tb1),tb2,tb3),form=list(data=dd,start=2009,pref="scen"),rest=rest,exoadd=exoadd)


tbnames <- c("BVP ir jo dalys","Kainos","Darbo rinkos rodikliai")

scen1 <- scen.2.xml(scen,1,"Scenarijus 1",tbnames)
scen2 <- scen.2.xml(scen,2,"Scenarijus 2",tbnames)
scen3 <- scen.2.xml(scen,3,"Scenarijus 3",tbnames)
rest <- rest.2.xml(scen)
exoadd <- exoadd.2.xml(scen)

xml <- "<lemasta>"
xml <- paste(xml,scen1,sep="")
xml <- paste(xml,scen2,sep="")
xml <- paste(xml,scen3,sep="")
xml <- paste(xml,rest,sep="")
xml <- paste(xml,exoadd,sep="")
xml <- paste(xml,"</lemasta>",sep="")
write(xml,file="output/initial.xml")

save(ee,eqR,eqR2,ladt,exo2y,q2y,exoadd2y,exoameta,mreal,mnom,mprice,mwf,tbnames,fp.prep,fp.prep2,metamodel,file="output/lemasta.RData")

system("cp 10code.R output/code.R")


