#' roc for model
#'
#' @param ... one or more fit
#'
#' @return roc dataframe
#' @export
#'
roc <- function(...) UseMethod('roc')
