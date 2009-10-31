source("30prepare.R")


system.time(ftry <- eqforecast(start=c(2009,1),end=c(2011,4),eqR,ee,data=ladt,leave=TRUE,use.jacobian=TRUE))

#system.time(ftry <- eqforecast(start=c(2000,1),end=c(2008,4),eqR,ee,data=ladt,leave=FALSE,lower=lower,upper=upper,method="L-BFGS-B",control=list(trace=1)))



endol <- ladt[,colnames(ftry$res)]
gofd <- cbind(endol,ftry$res)
colnames(gofd) <- c(colnames(endol),paste("f_",colnames(ftry$res),sep=""))

no <- dim(endol)[2]

xyplot(gofd,screens=rep(1:no,2),col=c(rep("blue",no),rep("red",no)))

####Produce tables for comparison

q2y <- read.csv("tables/q2y.csv")

flydt <- q2y.meta(ftry$data,q2y)

mreal <- read.csv("tables/tb1real.csv")
mnom <-  read.csv("tables/tb1nom.csv")

tbreal1 <- produce.tb(flydt,mreal,years=2006:2011,gdpshare=as.character(mreal[1,2]))

tbnom1 <- produce.tb(flydt,mnom,years=2006:2011,as.character(mnom[1,2]))
