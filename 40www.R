source("30prepare.R")

################################################################
##Prepare exogenous variables (convert to yearly and back to quarterly)
exo2y <- read.csv("tables/exo2y.csv")
q2y <- read.csv("tables/q2y.csv")

lydt <- q2y.meta(ladt,q2y)

exotb <-  produce.tb(lydt,exo2y)
scy <- inverse.tb(exotb$level,exo2y)

scq <- y2q.meta(scy,exo2y)

ladt <- introduce.exo(scq,ladt,exo2y)

###############################################################
###Prepare initial tables for www

##Do forecast
system.time(ftry <- eqforecast(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt,leave=TRUE,use.jacobian=TRUE))

flydt <- q2y.meta(ftry$data,q2y)

###Read meta data
mreal <- read.csv("tables/tb1real.csv")
mnom <-  read.csv("tables/tb1nom.csv")

tbreal1 <- produce.tb(flydt,mreal,years=2006:2011,gdpshare=as.character(mreal[1,2]))

tbnom1 <- produce.tb(flydt,mnom,years=2006:2011,as.character(mnom[1,2]))


csvhtpair(tbreal1$growth,1,"varno",catalogue="output/")
csvhtpair(tbreal1$growth,2,"varno",catalogue="output/")
csvhtpair(tbreal1$growth,3,"varno",catalogue="output/")

############################################################
##Produce forms

dd <- exotb$level

a1 <- produce.form(dd,start=2009,pref="scen1")
a2 <- produce.form(dd,start=2009,pref="scen2")
a3 <- produce.form(dd,start=2009,pref="scen3")

print(xtable(a1),type="html",include.rows=FALSE,file="output/form1.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

print(xtable(a2),type="html",include.rows=FALSE,file="output/form2.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

print(xtable(a3),type="html",include.rows=FALSE,file="output/form3.html",html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)


save(ee,eqR,ladt,exo2y,q2y,mreal,mnom,file="output/lemasta.RData")

system("cp 10code.R output/code.R")


