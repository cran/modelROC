#' @method auc cph
#' @rdname auc
#'
#' @export
auc.cph <- function(...,model=NULL,x=NULL,method=c('NNE','KM')){
    auc.coxph(...,model=model,x=x,method=method)
}


#' @param x can be logical or characters. TRUE means all x variable in regression
#'     will be calculated. One or more characters will be calculated only.
#' @param model can be logical or characters. FALSE means no model TP and FP,
#'     characters mean model names.
#' @param method NNE or KM
#' @rdname auc
#'
#' @return one auc_coxph for cox regression. model means model names,
#' @export
#' @method auc coxph
auc.coxph <- function(...,model=NULL,x=NULL,method=c('NNE','KM')){
    method=match.arg(method)
    fitname <- do::get_names(...)
    if (isFALSE(model)) model=NULL
    if (isFALSE(x)) x= NULL
    if (isTRUE(model)) model= rep(TRUE,length(fitname))
    if (!is.null(model) &length(fitname) != length(model)) stop(tmcn::toUTF8("\u6709"),length(fitname),tmcn::toUTF8("\u4E2A\u6A21\u578B,\u4F46\u6709"),length(model),tmcn::toUTF8("\u4E2Amodel\u540D\u79F0"))
    lp <- lapply(fitname, function(i) auci(fiti=i,
                                                     modeli=model[fitname==i],
                                                     x=x,
                                                     method=method))
    pp <- do.call(rbind,lp)
    class(pp) <- c('auc_coxph','data.frame')
    pp
}
auci <- function(fiti,modeli=NULL,x=NULL,method=c('NNE','KM')){
    method=match.arg(method)
    fitg <- get(fiti,envir = .GlobalEnv)
    data <- eval(fitg$call$data)
    vtime <- data[,do::model.y(fitg)[1]]
    vstatus <- data[,do::model.y(fitg)[2]]
    linerpredictor <- data.frame(model=exp(fitg$linear.predictors))
    if (is.logical(x[1])){
        if (x[1]){
            x <- do::model.x(fitg)
        }else{
            x <- NULL
        }
    }
    x <- x[ x %in% do::model.x(fitg)]
    if (is.logical(modeli)){
        if (modeli){
            if (!is.null(x) & (fiti %in% x)) stop(tmcn::toUTF8("model\u548Cx\u4E0D\u80FD\u6709\u540C\u540D:"),fiti)
            vx <- c(fiti,x)
            xmt <- cbind(linerpredictor,data[,x,drop=FALSE])
            colnames(xmt) <- vx
        }else{
            if (is.null(x)) stop(tmcn::toUTF8("x\u548Cmodel\u4E0D\u80FD\u540C\u65F6\u4E3ANULL"))
            vx <- x
            xmt <- data[,x,drop=FALSE]
            colnames(xmt) <- vx
        }
    }else{
        if (is.null(modeli)){
            if (is.null(x)){
                stop(tmcn::toUTF8("model\u548Cx\u4E0D\u80FD\u540C\u65F6\u4E3ANULL"))
            }else{
                vx <- x
                xmt <- data[,x,drop=FALSE]
                colnames(xmt) <- vx
            }
        }else{
            if (!is.null(x) & (modeli %in% x)) stop(tmcn::toUTF8("model\u548Cx\u4E0D\u80FD\u6709\u540C\u540D:"),modeli)
            vx <- c(modeli,x)
            xmt <- cbind(linerpredictor,data[,x,drop=FALSE])
            colnames(xmt) <- vx
        }

    }
    head(xmt)
    # x is not null
    lp <- lapply(vx, function(j){
        pb <- txtProgressBar(max = length(unique(vtime)),width = 30,style=3)
        cat('   ',j)
        lpx <- lapply(1:length(unique(vtime)), function(i){
            setTxtProgressBar(pb,value = i)
            r <- survivalROC::survivalROC(Stime=vtime,
                                          status=vstatus,
                                          marker = xmt[,j],
                                          predict.time =unique(vtime)[i],
                                          method=method,
                                          span = 0.25*NROW(data)^(-0.20))
            Yd <- ifelse(r$AUC >= 0.5,
                         paste0(round(r$cut.values[which.max(r$TP-r$FP)],3),collapse = ', '),
                         paste0(round(r$cut.values[which.min(r$TP-r$FP)],3),collapse = ', '))

            data.frame(model=fiti,
                       time=r$predict.time,
                       marker=j,
                       AUC=r$AUC,
                       Youden = Yd
            )
        })
        close(pb)
        do.call(rbind,lpx)
    })
    do.call(rbind,lp)
}









