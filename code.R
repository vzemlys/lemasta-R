require(tseries)
require(lattice)

eqforecast <- function(start,end,eq,endoexo,data,...) {
    exonames <- as.character(endoexo$name[endoexo$exo=="Exog"])

    noendog <- table(endoexo$exo)["Endog"]
    subtb <- cbind(as.character(endoexo$name[endoexo$exo=="Endog"]),paste("y[",1:noendog,"]",sep=""))
    
    mod.o <- function(y) {
        as.numeric(eval(eqoy,list(y=y)))
    }

    mod.ograd <- function(y) {
        res <- sapply(eqogy,function(l)as.numeric(eval(l,list(y=y))))
        res
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
       # browser()
        
        eqs <- eqs.optim(eqmod,subtb)
        eqoy <- eqs$fn
        eqogy <- eqs$grad
        
        x0 <- as.numeric(window(lag(data[,subtb[,1]],-1),start=it0,end=it0))
        ##Submit  lag as a starting value, since at current time the values might not exist

        
        fogs <- optim(par=log(x0),fn=mod.o,gr=mod.ograd,...)
        
        ind <- as.numeric(window(indt,start=it0,end=it0))
        ###If some values were not NA, leave them intact
        fc <- exp(fogs$par)
        cd <- data[ind,subtb[,1]]
        if(sum(!is.na(cd))>0) {
            fc[!is.na(cd)] <- cd[!is.na(cd)]
        }
        
        data[ind,subtb[,1]] <- fc
        res <- rbind(res,fc)
    }
    res <- ts(res,start=start,end=end,frequency=4)
    colnames(res) <- subtb[,1]
    list(res=res,data=data)
}

eqs.optim <- function(eqmod,subtable) {
    eqo <- sapply(eqmod,function(l) {
        aa <- paste(deparse(l,wid=500),collapse="")
        paste("(",aa,")^2",sep="")
    })

    eqo <- parse(text=paste(eqo,collapse="+"))

    eqo <- subvars(eqo[[1]],cbind(subtable[,1],subtable[,1]),make.exp=TRUE)

    eqog <- lapply(subtable[,1],function(l)D(eqo,name=l))

    eqoy <- subvars(eqo,subtable)
    eqogy <-lapply(eqog,function(l)subvars(l,subtable))

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


produce.tb1 <- function(data) {
    f <- list(expression((y_r_sa/lag(y_r_sa,-1)-1)*100),
              expression(w_b_sa/lag(w_b_sa,-1)*100),
              expression(w_b_sa),
              expression((x_r_sa-m_r_sa)/y_r_sa),
              expression((c_r_sa/lag(c_r_sa,-1)-1)*100),
              expression((i_r_sa/lag(i_r_sa,-1)-1)*100),
              expression((y_n_sa/lag(y_n_sa,-1)-1)*100)
              )
    nf <- c("BVP augimas, grandine susietos apimties augimas, proc.",
            "Vidutinio mėnesinio bruto darbo užmokesčio indeksai, ankstesnis laikotarpis = 100",
            "Vidutinis mėnesinis bruto darbo užmokestis, Lt",
            "Prekių ir paslaugų balansas, proc. BVP",
            "Vartojimo augimas, grandine susietos apimties augimas, proc.",
            "Bendrojo pagrindinio kapitalo formavimo augimas, grandine susietos apimties augimas, proc.",
            "BVP augimas to meto kainomis, proc.")
    
    tld <- data
    tld <- cbind(qpadd(start(data),end(data)),tld)

    colnames(tld) <- c("year","quarter",colnames(data))

    agm <- recast(as.data.frame(tld),year~variable,mean,id.var=c("year","quarter"))
    ags <- recast(as.data.frame(tld),year~variable,sum,id.var=c("year","quarter"))

    agm <- ts(agm[,-1],start=start(data)[1])
    ags <- ts(ags[,-1],start=start(data)[1]) ## Will not be used for the mooment
    
    rl <- lapply(f,eval,envir=as.list(agm))
    res <- t(sapply(rl,window,start=2008,end=2012))
    res <- data.frame(Rodiklis=nf,res)
    names(res)[-1] <- 2008:2012
    res
}
