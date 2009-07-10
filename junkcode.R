str1 <- "DLOG(C_R_SA) = 2.41488661147209 - 0.73340228709126*LOG(C_R_SA( - 1)) + 0.23722375445228*LOG(Y_DU_SA/P_C_SA) + 0.132143509191045*LOG(Y_MISR_SA( - 1)/P_C_SA( - 1)) + 0.0685163815802085*LOG(B_L_SA) + 0.0522760107429614*LOG(QM_SA) + 0.21150950597883*DLOG(Y_MISR_SA/P_C_SA) - 0.00422405406529888*D(I_L,0,2)"

eq1 <- parse(text=tolower(str1))


eq <- parse(text="dlog(y)=1+x(-1)")

eq <- parse(text="bquote(y-.(l(y,-1))+x)")

eq <- bquote(y-.(l(y,-1))+x)

txt <- "zz(-1)"

gsub("(zz)([(])","l\\2\\1,",txt)

exdl <-parse(text="dlog(y)-1-dlog(x)")

exd <- parse(text="dlog(y)-1-dlog(x(-1),0,2)-d(z,0,2)")

exl <- parse(text="log(y(-1)/y(-2))-1-x(-1)+log(z)")

lx <- edlogv(exl[[1]],varnames=c("y","x","z"))

blx <- edlag(lx,start=3,end=3)

blz <- edlagv(lx,start=3,end=3,exonames=c("x","z"))

bq <- parse(text=paste("bquote(",deparse(blz,width=500),")",sep=""))


exonames <- as.character(endoexo$name[endoexo$exo=="Exog"])

it0 <- c(2008,4)

eqit <- lapply(eqR,function(l)edlagv(l,start=it0,end=it0,exonames=exonames))

eqbq <- lapply(eqit,function(l) {
    parse(text=paste("bquote(",paste(deparse(l,width=500),collapse=""),")",sep=""))
})

eqmod <- lapply(eqbq,function(l)eval(l,as.list(ladt)))

noendog <- table(endoexo$exo)["Endog"]

subtb <- cbind(as.character(endoexo$name[endoexo$exo=="Endog"]),paste("y[",1:noendog,"]",sep=""))

eqdfs <- lapply(eqmod,function(l)subvars(l,subtb,make.exp=TRUE))

modelis <- function(y,eqs) {
    res <- sapply(eqs,function(l)as.numeric(eval(l,list(y=y))))
    res
}

mod.optim <- function(y,eqs) {
        res <- sapply(eqs,function(l)as.numeric(eval(l,list(y=y))))
        sum(res^2)
}

eqo <- sapply(eqmod,function(l) {
    aa <- paste(deparse(l,wid=500),collapse="")
    paste("(",aa,")^2",sep="")
})

eqo <- parse(text=paste(eqo,collapse="+"))

eqo <- subvars(eqo[[1]],cbind(subtb[,1],subtb[,1]),make.exp=TRUE)

eqog <- lapply(subtb[,1],function(l)D(eqo,name=l))

eqoy <- subvars(eqo,subtb)
eqogy <-lapply(eqog,function(l)subvars(l,subtb))

   
mod.o <- function(y) {
    as.numeric(eval(eqoy,list(y=y)))
}

mod.ograd <- function(y) {
    res <- sapply(eqogy,function(l)as.numeric(eval(l,list(y=y))))
    res
}


  

x0 <- as.numeric(window(ladt[,subtb[,1]],start=it0,end=it0))

require(BB)

#firstry <- dfsane(par = log(x0), fn = modelis, method=3,control=list(M=10,tol = 1.e-10),eqs=eqdfsane)

fo <- optim(log(x0),mod.optim,eqs=eqdfs,method="BFGS",control=list(trace=1))

fog <- optim(par=log(x0),fn=mod.o,gr=mod.ograd,method="BFGS",control=list(trace=1))

lower <- rep( -Inf,26)
upper <- rep(Inf,26)

upper[19] <- log(100)

fogs <- optim(par=log(x0),fn=mod.o,gr=mod.ograd,method="L-BFGS-B",control=list(trace=1))



tld <- ladt
tld <- cbind(qpadd(start(ladt),end(ladt)),ladt)

colnames(tld) <- c("year","quarter",colnames(ladt))
