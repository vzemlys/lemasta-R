

eqforecast <- function(start,end,eq,endoexo,data,leave=TRUE,use.jacobian=TRUE,...) {
    exonames <- as.character(endoexo$name[endoexo$exo=="Exog"])

    noendog <- table(endoexo$exo)["Endog"]
    subtb <- cbind(as.character(endoexo$name[endoexo$exo=="Endog"]),paste("y[",1:noendog,"]",sep=""))
    
    mod.fn <- function(x) {
        sapply(eqosy,function(l)as.numeric(eval(l,list(y=x))))
    }

    mod.jac <- function(x) {
        sapply(eqosgy,function(l)sapply(l,function(ll)as.numeric(eval(ll,list(y=x)))))
    }
        
        
    timem <- qpadd(start,end)
    indt <- ts(1:dim(data)[1],start=start(data),end=end(data),freq=4)
    res <- numeric()
    for (i in 1:dim(timem)[1])     {

        it0 <- timem[i,]
       
        
        eqit <- lapply(eq,function(l)edlagv(l,start=it0,end=it0,exonames=exonames))
        
        eqit <- lapply(eqit,function(l) {
            parse(text=paste("bquote(",paste(deparse(l,width=500),collapse=""),")",sep=""))
        })
        
        
        eqmod <- lapply(eqit,function(l)eval(l,as.list(data)))

        eqq <- eqs.nleqslv(eqmod,subtb,use.jacobian=use.jacobian)
        
        eqosy <- eqq$fn
        eqosgy <- eqq$grad

        x0 <- as.numeric(window(lag(data[,subtb[,1]],-1),start=it0,end=it0))


        if(use.jacobian)    
          fogs <- nleqslv(x0,mod.fn,jac=mod.jac,...)
        else
          fogs <- nleqslv(x0,mod.fn,...)
     
        ind <- as.numeric(window(indt,start=it0,end=it0))
        ###If some values were not NA, leave them intact
      
        fc <- fogs$x
        if(leave) {
            cd <- data[ind,subtb[,1]]
            if(sum(!is.na(cd))>0) {
                fc[!is.na(cd)] <- cd[!is.na(cd)]
            }
        }
        data[ind,subtb[,1]] <- fc
        res <- rbind(res,fc)

    }
    res <- ts(res,start=start,end=end,frequency=4)
    colnames(res) <- subtb[,1]
    list(res=res,data=data)
}

eqs.nleqslv <- function(eqmod,subtable,use.jacobian=TRUE) {

    eqoy <-lapply(eqmod,function(l)subvars(l,subtable))

    if(use.jacobian) {
        eqog <- lapply(subtable[,1],function(l)lapply(eqmod,function(ll,nm=l)D(ll,name=nm)))
    

        eqogy <- lapply(eqog,function(l)lapply(l,function(ll)subvars(ll,subtable)))
    }
    else
      eqogy <- NULL
    
    list(fn=eqoy,grad=eqogy)
}

eviewstoR <-function(str,varnames) {
    varnames <- tolower(varnames)
    str <- tolower(str)
    
    eq <- parse(text=str)
    lnames <- intersect(varnames,all.names(eq))
    res <- edlogv(eq[[1]],varnames=lnames)
    res
}

read.eviews <- function(file) {
    ###Reads eviews equations output, strips
    ###lines starting with @INNOV, and removes @IDENTITY
    ss <-scan(file,what=character(),sep='\n',blank.lines.skip=TRUE,comment.char="'")
    ss <- gsub("@INNOV.*","",ss)
    ss <- gsub("@IDENTITY +","",ss)
    ss <- ss[ss!=""]
    ss <- gsub("[[].*[]]","",ss) ###Remove error terms in brackets
   
    ss <- gsub("[+-] *$","",ss) ### Remove trailing arithmetic signs at the end of the line
    ss
}

qpadd <- function(start,end) {
##Padd in quarters given start and end in terms c(year,quarter)
    q <- start[2]:4
    
    y <- rep(start[1],length(q))

    if(end[1]>start[1]+1) {

        my <- rep((start[1]+1):(end[1]-1),each=4)
        q <- c(q,rep(1:4,length.out=length(my)))
        y <- c(y,my)
    }
    if(end[1]>start[1]) {
        eq <- 1:end[2]
        y <- c(y,rep(end[1],length(eq)))
        q <- c(q,eq)
    }
    res <- cbind(y,q)
    colnames(res) <- NULL
    rownames(res) <- NULL
    res
}

simplefc <- function(data,tb,ahead) {
    require(tseries)
    require(forecast)
    
    ff <- numeric()
    for(nm in colnames(data)) {
        trans <- tb$FT[tb$name==nm]
        if(trans=="skip") {
            res <- rep(data[dim(data)[1],nm],ahead)
        }
        else {
            if(trans=="log") {
                res <- forecast(auto.arima(log(na.omit(data[,nm]))),h=ahead)
                res <- exp(res$mean)
            }
            else {
                res <- forecast(auto.arima(na.omit(data[,nm])),h=ahead)
                res <- res$mean
            }
        }
        ff <- cbind(ff,res)
    }
    colnames(ff) <- colnames(data)
    ff
}
plot.forecast <- function(x,fc,varn,labels) {

    ##Find the variables who do not have seasonality removed.
   
    svarn <- paste(varn,"_sa",sep="")
    wsn <- setdiff(svarn,colnames(x))

    if(length(wsn)>0) {
        wsnr <- sapply(wsn,function(s)substring(s,1,nchar(s)-3))
        svarn <- c(setdiff(svarn,wsn),wsnr)
    }
        
    fcn <- intersect(colnames(fc),svarn)
    
    if(length(fcn)>0) {
        endol <- x[,svarn,drop=FALSE]
        gofd <- cbind(endol,fc[,fcn])

        colnames(gofd) <- c(svarn,paste("f_",fcn,sep=""))
        no <- dim(endol)[2]
        fscr <- match(fcn,svarn)
        ylab=as.character(labels[match(varn,as.character(labels[,1])),3])
        ylab[is.na(ylab)] <- "a"
        
        print(xyplot(gofd,screens=c(1:no,fscr),col=c(rep("blue",no),rep("red",length(fcn)))))
    }
    else {
        print(xyplot(x[,svarn,drop=FALSE]))
    }


}

edlogv <- function(expr,varnames) {
    if(length(expr)==1) {
        return(expr)
    }
    if(length(expr)==2) {
        if(expr[[1]]==as.name("dlog")) {
            ie <- deparse(edlogv(expr[[2]],varnames))
            e <- paste("log(",ie,")-lag(log(",ie,"),-1)",sep="")
            return(parse(text=e)[[1]])
        }
        else {
            nm <- deparse(expr[[1]])
            if(nm %in% varnames) {
                ##Goodie, we have a lag!
                nol <- deparse(expr[[2]])
                e <- paste("lag(",nm,",",nol,")",sep="")
                return(parse(text=e)[[1]])
            }
            else {
                expr[[2]] <- edlogv(expr[[2]],varnames)
            }
        }
    }
    if(length(expr)>=3) {
        if(expr[[1]]==as.name("d")) {
                ie <- deparse(edlogv(expr[[2]],varnames))
                nol <- deparse(expr[[4]])
                e <- paste(ie,"-lag(",ie,",-",nol,")",sep="")
                return(parse(text=e)[[1]])
            }
        else {
            if(expr[[1]]==as.name("dlog")) {
                ie <- deparse(edlogv(expr[[2]],varnames))
                nol <- deparse(expr[[4]])
                e <- paste("log(",ie,")-lag(log(",ie,"),-",nol,")",sep="")
                nn <- as.numeric(deparse(expr[[3]]))
                if(nn>0) {
                    for(i in 1:nn) {
                        ie <- e
                        e <- paste(ie,"-lag(",ie,",-1)",sep="")
                    }
                }
                  
                return(parse(text=e)[[1]])
            }
            else {
                expr[[2]] <- edlogv(expr[[2]],varnames)
                expr[[3]] <- edlogv(expr[[3]],varnames)
            }
        }
    }
        
    
    
    return(expr)
}

subvars <- function(expr,tb,make.exp=FALSE) {
    if(length(expr)==1) {
        nm <- deparse(expr)
        if(nm %in% tb[,1]) {
            if(make.exp) {
                ie <- tb[tb[,1]==nm,2]
                e <- paste("exp(",ie,")",sep="")
            }
            else {
                e <- tb[tb[,1]==nm,2]
            }
            return(parse(text=e)[[1]])
        }
        else
          return(expr)
    }
    else {
        for (i in 2:length(expr)) {
            expr[[i]] <- subvars(expr[[i]],tb,make.exp)
        }
          
    }
        
    return(expr)
}



edlagv <- function(expr,start,end,exonames) {
    if(length(expr)==1) {
        nm <- deparse(expr)
        if(nm %in% exonames) {
            e <- paste(".(as.numeric(window(",nm,",start=",deparse(start),",end=",deparse(end),")))",sep="")
            return(parse(text=e)[[1]])
        }
        else         return(expr)
    }
    if(length(expr)==2) {
        expr[[2]] <- edlagv(expr[[2]],start,end,exonames)
    }
    if(length(expr)==3) {
        if(expr[[1]]==as.name("lag")) {
            ie <- deparse(expr)
            e <- paste(".(as.numeric(window(",ie,",start=",deparse(start),",end=",deparse(end),")))",sep="")
            return(parse(text=e)[[1]])
        }
        else {
            expr[[2]] <- edlagv(expr[[2]],start,end,exonames)
            expr[[3]] <- edlagv(expr[[3]],start,end,exonames)
        }
    }
    return(expr)
}


edlag <- function(expr,start,end) {
    if(length(expr)==1) {
        return(expr)
    }
    if(length(expr)==2) {
        expr[[2]] <- edlag(expr[[2]],start,end)
    }
    if(length(expr)==3) {
        if(expr[[1]]==as.name("lag")) {
            ie <- deparse(expr)
            e <- paste(".(window(",ie,",start=",deparse(start),",end=",deparse(end),"))",sep="")
            return(parse(text=e)[[1]])
        }
        else {
            expr[[2]] <- edlag(expr[[2]],start,end)
            expr[[3]] <- edlag(expr[[3]],start,end)
        }
    }
    return(expr)
}

edlog <- function(expr) {
    if(length(expr)==1) {
        return(expr)
    }
    if(length(expr)==2) {
        if(expr[[1]]==as.name("dlog")) {
            ie <- deparse(edlog(expr[[2]]))
            e <- paste("log(",ie,")-lag(log(",ie,"),-1)",sep="")
            return(parse(text=e)[[1]])
        }
        else return(expr)

        
    }
    if(length(expr)>=3) {
        if(expr[[1]]==as.name("d")) {
                ie <- deparse(edlog(expr[[2]]))
                nol <- deparse(expr[[4]])
                e <- paste(ie,"-lag(",ie,",-",nol,")",sep="")
                return(parse(text=e)[[1]])
            }
        else {
            if(expr[[1]]==as.name("dlog")) {
                ie <- deparse(edlog(expr[[2]]))
                nol <- deparse(expr[[4]])
                e <- paste("log(",ie,")-lag(log(",ie,"),-",nol,")",sep="")
                return(parse(text=e)[[1]])
            }
            else {
                expr[[2]] <- edlog(expr[[2]])
                expr[[3]] <- edlog(expr[[3]])
            }
        }
    }
        
    
    
    return(expr)
}

q2y.fun <- function(data,fun.aggregate=sum)  {
     tld <- data
     tld <- cbind(qpadd(start(data),end(data)),tld)
     colnames(tld) <- c("year","quarter",colnames(data))

     ag <- recast(as.data.frame(tld),year~variable,fun.aggregate=fun.aggregate,id.var=c("year","quarter"))
     ts(ag[,-1],start=start(data)[1])
}

q2y.stream <- function(data) {
     q2y.fun(data)
}

q2y.index <- function(data) {
    q2y.fun(data,function(x)x[length(x)])
}



q2y.meta <- function(data,meta) {
    
    q2s <- q2y.stream(data[,as.character(meta$Name[meta$Type=="stream"]),drop=FALSE])
    q2i <- q2y.index(data[,as.character(meta$Name[meta$Type=="index"]),drop=FALSE])

    lydt <- cbind(q2s,q2i)
    colnames(lydt) <- c(colnames(q2s),colnames(q2i))
    lydt
}

y2q <- function(data,fun.expand) {
    res <- apply(data,2,function(x)sapply(x,fun.expand))
    if(ncol(data)==1) {
        res <- matrix(res,ncol=1)
        colnames(res) <- colnames(data)
    }
    ts(res,start=c(start(data),1),frequency=4)
}

y2q.stream <- function(data) {
    y2q(data,function(x)rep(x/4,4))
}

y2q.index <- function(data) { 
    y2q(data,function(x)c(NA,NA,NA,x))
}

y2q.meta <- function(data,meta) {
    
    y2s <- y2q.stream(data[,as.character(meta$Name[meta$Type=="stream"]),drop=FALSE])
    y2i <- y2q.index(data[,as.character(meta$Name[meta$Type=="index"]),drop=FALSE])

    lydt <- cbind(y2s,y2i)
    colnames(lydt) <- c(colnames(y2s),colnames(y2i))
    lydt
}

introduce.exo <- function(x,data,meta,start=c(2009,2)) {
    
    ee <- end(x)
    if(start[2]!=1) {
        si <- start
        ss <- c(start[1],1)
    }
    st <- si
    tmp <- foreach(l=as.list(x),nm=colnames(x)) %do% {
        if(meta[meta$Name==nm,"Type"]=="stream")
          st <- ss
        else
          st <- si
        
        window(data[,nm],start=st,end=ee) <- window(l,start=st,end=ee)
        fillna <- na.approx(data[,nm])
        window(data[,nm],start=start(fillna)) <- fillna
       
    }
    data
}

inverse.tb <- function(x,meta) {
    #x turi būti lentelė su stulpeliu rodiklis, kuriame yra pavadinimai
    #atitinkantys meta stulpelio rodiklis pavadinimus.
    #kituose stulpeliuose turi būti metiniai duomenys.
    #meta turi turėti stulpelius rodiklis, pavadinimas ir inverse
    #rezultatas bus metinė laiko eilutė kurioje vienetai bus
    #tokie patys kaip ir originaliuose duomenyse.

    years <- as.numeric(colnames(x)[-2:-1])

    x <- x[match(as.character(meta$Rodiklis),as.character(x$Rodiklis)),]
    nm <- as.character(meta$Name)

    tt <- ts(t(x[,-2:-1]),start=years[1])
    colnames(tt) <- nm
           
    formulas <- as.character(meta$Inverse)
    
    level <- foreach(l=formulas) %do% eval(parse(text=l),envir=as.list(tt))

    res <- ts(sapply(level,function(x)x),start=start(tt))
    colnames(res) <- nm
    res
}

produce.tb <- function(x,meta,years=2006:2011,gdpshare=NULL){
   
    formulas <-as.character(meta[,2])
    fnames <- as.character(meta[,1])
    
      
    years <- min(years):max(years)
  
    level <- foreach(l=formulas) %do% eval(parse(text=l),envir=as.list(x))

    growth <- foreach(el=level) %do% { (el/lag(el,-1)-1)*100}

    if(!is.null(gdpshare)) {
        no <- which(formulas==gdpshare)
        gdp <- level[[no]]
        gdpshare <- foreach(el=level) %do% {(el/gdp)*100 }
    }
    

    res <- list(level=level,growth=growth)
    res$gdpshare <- gdpshare

    nmr <- names(res)
    
    res <- foreach(el=res) %do% {
        ll <- t(sapply(el,window,start=min(years),end=max(years)))
        dt <- data.frame(Rodiklis=fnames,Vienetai=as.character(meta$Vienetai),ll)
        names(dt)[-2:-1] <- years
        dt
    }
    
    names(res) <- nmr

    res
    
}

csvht <- function(tb,scenno,startvar=0,var="varno",scen="scenno",catalogue="") {
    require(xtable)
    dtln <- paste(catalogue,"datalev",suffix,".csv",sep="")
    dtgr <- paste(catalogue,"datagro",suffix,".csv",sep="")
    htname <- paste(catalogue,"ftable",suffix,".html",sep="")

    dml <- dim(level)
    dmg
}

csvhtpair <- function(res,suffix,cssattr="varno",scenattr="scenno",catalogue="") {
##cssattr is very important, since it is used in javascript code
##to determine which variable is being compared
    
    require(xtable)
    dtname <- paste(catalogue,"data",suffix,".csv",sep="")
    htname <- paste(catalogue,"ftable",suffix,".html",sep="")
    write.table(res,file=dtname,row.names=FALSE,quote=FALSE,sep="\t")

    ids <- 1:dim(res)[1]-1
    butt <- paste("<input type='submit' value='Lyginti' ",cssattr,"='",ids,"'/>",sep="")

    radiol <- paste("<input type='radio' name='show",suffix,"type",ids,"' value='Level' ",cssattr,"='",ids,"' ",scenattr,"='",suffix, "' checked'/>",sep="")
    radiog <- paste("<input type='radio' name='show",suffix,"type",ids,"' value='Growth' ",cssattr,"='",ids,"' ",scenattr,"='",suffix, "'/>",sep="")
    #radio <- paste(radiol,radiog,sep="")
    
    hres <- data.frame(res,radiol,radiog,butt)
    colnames(hres) <- c(names(res),"Lygis","Augimas","")

    print(xtable(hres),type="html",include.rows=FALSE,file=htname,html.table.attributes='border="1" id="table1",cellpading="2"',sanitize.text.function=function(x)x)

}

produce.form <- function(etb,start=2009,prefix="") {

    no <- dim(etb)[1]
    nms <- paste(prefix,"egzo",1:no,"[]",sep="")
    years <- as.numeric(colnames(etb)[-1])

    shy <- as.character(years[years<start])
    wry <- as.character(years[years>=start])
    
    write.col.fun <- function(col,nm) {
        col <- prettyNum(round(col,2))
        vals <- paste("<input name='",nm,"' value=",col," type='text' size='5'/>",sep="")
        vals
    }
    show.col.fun <- function(col,nm) {
        col <- prettyNum(round(col,2))
        vals <- paste(col,"<input name='",nm,"' value=",col," type='hidden'/>",sep="")
        vals
    }

    scols <- sapply(etb[,shy],show.col.fun,nm=nms)
    wcols <- sapply(etb[,wry],write.col.fun,nm=nms)
    
    first <- paste(etb[,1],"<input name='",nms,"' value='",etb[,1],"' type='hidden'/>",sep="")
    res <- data.frame(first,scols,wcols)

    names(res) <- c("Rodiklis",names(etb)[-1])
    res

}

todf <- function(x) {
    data.frame(row=2:length(x)-1,variable=x[1],value=as.numeric(x[-1]))

}

doforecast <- function(x,sceno,years=2006:2011) {
    require(tseries)
    require(nleqslv)
    require(plyr)
    require(reshape)
    require(foreach)

    colnames(x) <- c("Rodiklis",years)
    
    scy <- inverse.tb(x,exo2y)
    scq <- y2q.meta(scy,exo2y)
    dd <- introduce.exo(scq,ladt,exo2y)
    
    ftry <- eqforecast(start=c(2009,1),end=c(2011,4),eqR,ee,data=dd,leave=TRUE,use.jacobian=TRUE,control=list(ftol=1e-3))

    flydt <- q2y.meta(ftry$data,q2y)
    tbreal1 <- produce.tb(flydt,mreal,years=years,gdpshare=as.character(mreal[1,2]))
    tbnom1 <- produce.tb(flydt,mnom,years=years,as.character(mnom[1,2]))
    
    csvhtpair(tbreal1$growth,sceno);
}