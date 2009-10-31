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


dd <- t(ladt[,ee$name[ee$exo=="Exog"]])[,11:14]

dd <- data.frame(rownames(dd),dd)

str <- 'd971014;1;1;1;1;
d984994;0;0;0;0;
i_ilg;15.21333;11.97;13.45333;11.51333;
i_l;14.66333;12.59667;13.21;12.36667;
l_sa;1820.787;1810.955;1761.165;1717.401;
m_n_sa;6115.145;6772.301;6557.701;6566.46;
p_de_sa;0.9896704;0.9887405;0.9888156;0.9944268;
p_mino_sa;1.194639;1.187765;1.198492;1.169425;
p_m_sa;1.046018;1.051855;1.046247;0.9990417;
p_x_sa;0.972795;0.9756456;0.9400556;0.9166854;
t;11;12;13;14;
y_r_add;-99.67067;337.1502;-317.0353;70.86014;
y_r_de_sa;482027.9;484760.6;489101.7;487592.9;
'

todf <- function(x) {
    data.frame(row=2:length(x)-1,variable=x[1],value=as.numeric(x[-1]))

}


compare.fit <- function(fit,data) {
    require(ggplot2)
    ##fit and data are timeseries matrixes. Data may be have more
    ##columns than fit

    data <- data[,colnames(fit)]

    dt <- melt(ts.2.arr(data))
    ft <- melt(ts.2.arr(fit))
    
    colnames(dt)[1:2] <- c("time","variable") ->colnames(ft)[1:2]

    pp <- qplot(x=time,y=value,data=dt,geom="line")+geom_line(data=ft,colour="red")+facet_wrap(~variable,scales="free_y")
    pp
    
}

ts.2.arr <- function(x) {
    array(x,dim(x),dimnames=list(time(x),colnames(x)))
}

ch <- ladt[,c("c_r_sa_0","i_r_sa_0","y_nonx_sa_0")]
colnames(ch) <- c("c_r_sa","i_r_sa","y_nonx_sa")


gofdc <- cbind(gofd,ch)

colnames(gofdc) <- c(colnames(gofd),paste("x_",colnames(ftry$res),sep=""))

no <- dim(endol)[2]

xyplot(gofdc,screens=rep(1:no,3),col=c(rep("blue",no),rep("red",no),rep("green",no)))


cc <- cbind(ch,ladt[,setdiff(colnames(ladt),colnames(ch))])
colnames(cc) <- c(colnames(ch),setdiff(colnames(ladt),colnames(ch)))

gg <- ts(rbind(window(ladt[,colnames(ftry$res)],end=c(2006,4)),ftry$res),start=start(ladt),freq=frequency(ladt))

dd <- cbind(gg,ladt[,setdiff(colnames(ladt),colnames(gg))])
colnames(dd) <- c(colnames(gg),setdiff(colnames(ladt),colnames(gg)))

window(ladt[,c("c_r_sa_0","i_r_sa_0","y_nonx_sa_0")],start=c(2007,1))
