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
   
    timem <- fillstartend(start,end)
    res <- numeric()
    for (i in 1:dim(timem)[1])     {
        it0 <- timem[i,]
        
        eqit <- lapply(eq,function(l)edlagv(l,start=it0,end=it0,exonames=exonames))
        eqit <- lapply(eqit,function(l) {
            parse(text=paste("bquote(",paste(deparse(l,width=500),collapse=""),")",sep=""))
        })
        eqmod <- lapply(eqit,function(l)eval(l,as.list(data)))

        eqs <- eqs.optim(eqmod,subtb)
        eqoy <- eqs$fn
        eqogy <- eqs$grad
        
        x0 <- as.numeric(window(data[,subtb[,1]],start=it0,end=it0))

        
        fogs <- optim(par=log(x0),fn=mod.o,gr=mod.ograd,...)
        res <- rbind(res,exp(fogs$par))
    }
    res <- ts(res,start=start,end=end,frequency=4)
    colnames(res) <- subtb[,1]
    res
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
  
fillstartend <- function(start,end) {
    
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
    cbind(y,q)
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
