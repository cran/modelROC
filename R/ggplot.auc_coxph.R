#' @title Plot for ROC curve
#'
#' @param data resultes of roc() function
#' @param mapping ignore
#' @param color one or more colors
#' @param lwd logical or integers
#' @param ... ignore
#' @param grid.space space between grids, default is 2
#' @param ncol number of column for grid plot
#' @param environment ignore
#' @importFrom ggplot2 geom_line aes_string facet_wrap scale_color_manual xlab ylab scale_x_continuous scale_y_continuous geom_hline aes theme_bw theme unit element_blank geom_line aes_string facet_wrap facet_grid vars scale_color_manual xlab ylab scale_x_continuous scale_y_continuous coord_fixed theme_bw geom_abline theme unit geom_line aes_string facet_wrap scale_color_manual xlab ylab scale_x_continuous scale_y_continuous coord_fixed theme_bw geom_abline theme unit
#' @name ggplot
#' @method ggplot auc_coxph
#' @return a ggplot picture.
#' @export
#' @examples
#' \donttest{
#' library(ggDCA)
#' library(rms)
#' library(modelROC)
#' ###  COX ----------
#'
#' fit <- cph(Surv(time,status)~ANLN+CENPA+GPR182,LIRI)
#'
#' ####        one model, one time ####
#' #----            auc for model
#'
#' r <- auc(fit,
#'           model='This is model') # one model
#' unique(r)
#' ggplot(r)
#'
#' r <- auc(fit,
#'           model=TRUE)            # all model
#' unique(r)
#' ggplot(r)
#'
#'
#' #----            auc for x
#' r <- auc(fit,
#'           x='ANLN')              # one x
#' unique(r)
#' ggplot(r)
#'
#' r <- auc(fit,
#'           x=c('ANLN','CENPA'))   # more x
#' unique(r)
#' ggplot(r)
#'
#' r <- auc(fit,
#'           x=TRUE)                # all x
#' unique(r)
#' ggplot(r)
#'
#' r <- auc(fit,
#'           model=TRUE,            # one model
#'           x=TRUE)                # all x
#' unique(r)
#' ggplot(r)
#'
#' r <- auc(fit,
#'           model='Three Genes',   # specify model names
#'           x=TRUE)                # all X
#' unique(r)
#' ggplot(r)
#'
#'
#'
#' ####        more models ####
#' fit2 <- cph(Surv(time,status)~ANLN+CENPA,LIRI)
#'
#'
#' r <- auc(fit,fit2,
#'           model=c('Three Genes','Two Genes'))            #
#' unique(r)
#' ggplot(r)
#'
#'
#'
#' r <- auc(fit,fit2,
#'           model=TRUE,
#'           x=TRUE)
#' unique(r)
#' ggplot(r)
#' }
ggplot.auc_coxph <- function(data,
                             mapping,
                             color=NULL,
                             lwd=1.05,
                             grid.space=2,
                             ncol=NULL,
                             ...,
                             environment = parent.frame()){
    # data=pp
    pp <- as.data.frame(data)
    pp <- as.data.frame(pp)
    head(pp)
    mA <- unique(data)
    # length
    (lenmodel <- length(unique(pp$model)))
    (lenmarker <- length(unique(pp$marker)))

    # plot
    if (lenmodel == 1 & lenmarker == 1){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='time',y='AUC'))
    }else if(lenmodel == 1 & lenmarker > 1){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='time',y='AUC',color='marker'))
    }else if(lenmodel > 1 & lenmodel == lenmarker){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='time',y='AUC',color='marker'))
    }else if(lenmodel > 1 & lenmodel != lenmarker){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='time',y='AUC',color='marker')) +
            facet_wrap(~model,ncol = ifelse(is.null(ncol),length(unique(pp$model)),ncol))
    }
    if (!is.null(color)){
        p <- p+scale_color_manual(values = color)
    }
    p +
        xlab('Time')+
        ylab('AUC') +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        geom_hline(aes(yintercept=0.85), color="gray",
                   linetype="dashed", size=0.75)+
        geom_hline(aes(yintercept=0.70), color="gray",
                   linetype="dashed", size=0.75)+
        geom_hline(aes(yintercept=0.50), color="gray",
                   linetype="dashed", size=0.75)+
        theme_bw() +
        theme(legend.position = 'right',
              panel.spacing = unit(grid.space, "lines"),
              panel.grid = element_blank())
    # legend.title = element_blank()
}


