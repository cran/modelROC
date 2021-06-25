#' @method roc glm
#' @rdname roc
#' @importFrom stats median
#' @export
#' @references Pepe, Margaret Sullivan. \emph{The statistical evaluation
#' of medical tests for classification and prediction.} Medicine, 2003.
#'
#' @references Zou, Kelly H., W. J. Hall, and David E. Shapiro.
#' \emph{Smooth non-parametric receiver operating characteristic (ROC)
#'  curves for continuous diagnostic tests.} Statistics in medicine
#'  16, no. 19 (1997): 2143-2156.
roc.glm <- function(...,negref=0,model=NULL,x=NULL,method=c("empirical", "binormal","nonparametric")){
    roc.lrm(...,negref=negref,model=model,x=x,method=method)
}

#' @param negref negative reference for each model
#' @method roc lrm
#' @rdname roc
#' @export
#'
roc.lrm <- function(...,negref=0,model=NULL,x=NULL,method=c("empirical", "binormal","nonparametric")){
    method=match.arg(method)
    fitname <- do::get_names(...)
    if (isFALSE(model)) model=NULL
    if (isFALSE(x)) x= NULL
    if (is.null(model) & is.null(x)) stop(tmcn::toUTF8("model\u548Cx\u4E0D\u80FD\u540C\u65F6\u4E3ANULL\u6216\u8005FALSE"))
    if (isTRUE(model)) model= rep(TRUE,length(fitname))
    if (!is.null(model) & length(fitname) != length(model)) stop(tmcn::toUTF8("\u6709"),length(fitname),tmcn::toUTF8("\u4E2A\u6A21\u578B,\u4F46\u6709"),length(model),tmcn::toUTF8("\u4E2Amodel\u540D\u79F0"))
    if (length(negref)==1) negref <- rep(negref,length(fitname))
    if (length(fitname) != length(negref)) stop(tmcn::toUTF8("\u6709"),length(fitname),tmcn::toUTF8("\u4E2A\u6A21\u578B,\u4F46\u6709"),length(negref),tmcn::toUTF8("\u4E2Anegref"))
    lp <- lapply(fitname, function(i) lrmi(fiti=i,
                                           negref=negref[fitname==i],
                                           modeli=model[fitname==i],
                                           x=x,
                                           method=method))
    pp <- do.call(rbind,lp)
    class(pp) <- c('roc_logit','data.frame')
    pp
}
lrmi <- function(fiti,negref=0,modeli=NULL,x=NULL,method=c("empirical", "binormal","nonparametric")){
    method=match.arg(method)
    fitg <- get(fiti,envir = .GlobalEnv)
    data <- eval(fitg$call$data)
    class <- data[,do::model.y(fitg)[1]]
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
    # x is not null
    lp <- lapply(vx, function(j){
        r <- ROCit::rocit(score = xmt[,j],
                          class = class,
                          negref = negref,
                          method = method)
        Yd <- r$TPR-r$FPR
        Yd.max <- ifelse(r$AUC >= 0.5,
                     paste0(round(r$Cutoff[which.max(Yd)],3),collapse = ', '),
                     paste0(round(r$Cutoff[which.min(Yd)],3),collapse = ', '))

        data.frame(model=fiti,
                   marker=j,
                   AUC=r$AUC,
                   FP=r$FPR,
                   TP=r$TPR,
                   cutoff=r$Cutoff,
                   Youden = Yd,
                   Youden.max = Yd.max)
    })
    do.call(rbind,lp)
}









