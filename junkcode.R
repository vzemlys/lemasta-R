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

