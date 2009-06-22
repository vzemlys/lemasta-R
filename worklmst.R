

ad <- read.csv("all_data_all_data.csv")

sad <-read.csv("data_sa_Feuille1.csv")

adt <- ts(ad[,-1],start=c(1995,1),end=c(2009,4),freq=4)
sadt <- ts(sad[,-1],start=c(1995,1),end=c(2009,4),freq=4)

adt <- window(adt,start=c(1995,1),end=c(2009,1))
sadt <- window(sadt,start=c(1995,1),end=c(2009,1))

colnames(sadt) <- tolower(colnames(sadt))
colnames(adt)  <- tolower(colnames(adt))

  
source("eviewseq.R")
source("code.R")

eqstr <- gsub("\n","",eqstr)

eqstrm <- gsub("=","-",eqstr)


eqR <- lapply(eqstrm,function(l)eviewstoR(l,varnames=colnames(adt)))

endoexo <- read.csv("endoexo.csv")

endoexo <- data.frame(name=tolower(as.character(endoexo[,1])),exo=endoexo[,2])

exonames <- as.character(endoexo$name[endoexo$exo=="Exog"])

it0 <- c(2008,4)

eqit <- lapply(eqR,function(l)edlagv(l,start=it0,end=it0,exonames=exonames))

eqbq <- lapply(eqit,function(l) {
    parse(text=paste("bquote(",paste(deparse(l,width=500),collapse=""),")",sep=""))
})

eqmod <- lapply(eqbq,function(l)eval(l,as.list(adt)))

noendog <- table(endoexo$exo)["Endog"]

subtb <- cbind(as.character(endoexo$name[endoexo$exo=="Endog"]),paste("y[",1:noendog,"]",sep=""))

eqdfsane <- lapply(eqmod,function(l)subvars(l,subtb,make.exp=TRUE))

modelis <- function(y,eqs) {
    res <- sapply(eqs,function(l)as.numeric(eval(l,list(y=y))))
    res
}

x0 <- as.numeric(window(adt[,subtb[,1]],start=it0,end=it0))

require(BB)

firstry <- dfsane(par = log(x0), fn = modelis, method=3,control=list(M=10,tol = 1.e-10),eqs=eqdfsane)
