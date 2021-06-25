#' Extract data from roc() function
#'
#' @param x result of roc() function
#' @param incomparables ignore
#' @param ... ignore
#' @return a dataframe.
#' @name unique
#' @export
#' @method unique roc_coxph
unique.roc_coxph <- function(x, incomparables = FALSE, ...){
    x <- as.data.frame(x)
    unique(x[,c('model','time','marker','AUC','Youden.max')])
}
#' @rdname unique
#' @export
#' @method unique roc_logit
unique.roc_logit <- function(x, incomparables = FALSE, ...){
    x <- as.data.frame(x)
    unique(x[,c('model','marker','AUC','Youden.max')])
}
#' @rdname unique
#' @export
#' @method unique auc_coxph
unique.auc_coxph <- function(x, incomparables = FALSE, ...){
    x <- as.data.frame(x)
    unique(x[,c('model','marker')])
}
