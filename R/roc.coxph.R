#' @method roc cph
#' @rdname roc
#' @references Heagerty PJ, Lumley T, Pepe MS. \emph{Time-dependent ROC curves
#' for censored survival data and a diagnostic marker.} Biometrics, 2000.
#' @export
roc.cph <- function(...,times=NULL,model=NULL,x=NULL,method=c('NNE','KM')){
    roc.coxph(...,times=times,model=model,x=x,method=method)
}

#' @param times one or more times for cox regression
#' @param x can be logical or characters. TRUE means all x variable in regression
#'     will be calculated. One or more characters will be calculated only.
#' @param model can be logical or characters. FALSE means no model TP and FP,
#'     characters mean model names.
#' @param method NNE or KM
#' @rdname roc
#'
#' @return one roc_coxph for cox regression. model means model names,
#' @export
#' @method roc coxph
roc.coxph <- function(...,times=NULL,model=NULL,x=NULL,method=c('NNE','KM')){
    method=match.arg(method)
    fitname <- do::get_names(...)
    if (isFALSE(model)) model=NULL
    if (isFALSE(x)) x= NULL
    if (isTRUE(model)) model= rep(TRUE,length(fitname))
    if (!is.null(model) &length(fitname) != length(model)) stop(tmcn::toUTF8("\u6709"),length(fitname),tmcn::toUTF8("\u4E2A\u6A21\u578B,\u4F46\u6709"),length(model),tmcn::toUTF8("\u4E2Amodel\u540D\u79F0"))
    lp <- lapply(fitname, function(i) roci(fiti=i,times=times,
                                           modeli=model[fitname==i],
                                           x=x,
                                           method=method))
    pp <- do.call(rbind,lp)
    class(pp) <- c('roc_coxph','data.frame')
    pp
}
roci <- function(fiti,times=NULL,modeli=NULL,x=NULL,method=c('NNE','KM')){
    method=match.arg(method)
    fitg <- get(fiti,envir = .GlobalEnv)
    data <- eval(fitg$call$data)
    vtime <- data[,do::model.y(fitg)[1]]
    if (is.null(times)) times=median(vtime)
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
    lp <- lapply(times, function(i){
        lpx <- lapply(vx, function(j){
            r <- survivalROC::survivalROC(Stime=vtime,
                                       status=vstatus,
                                       marker = xmt[,j],
                                       predict.time =i,
                                       method=method,
                                       span = 0.25*NROW(data)^(-0.20))
            Yd <- r$TP-r$FP
            Yd.max <- ifelse(r$AUC >= 0.5,
                         paste0(round(r$cut.values[which.max(Yd)],3),collapse = ', '),
                         paste0(round(r$cut.values[which.min(Yd)],3),collapse = ', '))

            data.frame(model=fiti,
                       time=r$predict.time,
                       marker=j,
                       AUC=r$AUC,
                       FP=r$FP,
                       TP=r$TP,
                       cutoff=r$cut.values,
                       Youden = Yd,
                       Youden.max = Yd.max
                       )
        })
        do.call(rbind,lpx)
    })
    do.call(rbind,lp)
}









