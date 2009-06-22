l <- function(x,lag=-1){
    n <- length(x)
    x[n+lag]
}

eviewstoR.old <- function(str,varnames) {
    ###substitute a(-1) to l(a,-1) and dlog(a) to log(a)-l(log(a),-1),
    ###where a is the variable a

    varnames <- tolower(varnames)
    str <- tolower(str)
    
    eq <- parse(text=str)

    lnames <- intersect(varnames,all.names(eq))
    lnames
    for(l in lnames) {
        pattern <- paste("(",l,")([(])",sep="")
        str <- gsub(pattern,"l\\2\\1,",str)
        str
    }

    str
}

eviewstoR <-function(str,varnames) {
    varnames <- tolower(varnames)
    str <- tolower(str)
    
    eq <- parse(text=str)
    lnames <- intersect(varnames,all.names(eq))
    res <- edlogv(eq[[1]],varnames=lnames)
    res
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
