#' @importFrom ggplot2 ggplot
#' @rdname ggplot
#' @method ggplot roc_logit
#' @importFrom utils head setTxtProgressBar txtProgressBar
#' @export
#'
ggplot.roc_logit <- function(data,
                             mapping,
                             color=NULL,
                             lwd=1.05,
                             grid.space=2,
                             rank=FALSE,
                             ...,
                             environment = parent.frame()){
    # data=pp
    pp <- as.data.frame(data)
    pp <- as.data.frame(pp)
    mA <- unique(data)
    if (rank){
        level <- unique(mA$marker[order(mA$AUC,decreasing = TRUE)])
    }else{
        level <- unique(mA$marker)
    }
    pp$marker <-factor(pp$marker,levels = level)
    # length
    (lenmodel <- length(unique(pp$model)))
    (lenmarker <- length(unique(pp$marker)))

    # plot
    if (lenmodel == 1 & lenmarker == 1){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='FP',y='TP'))
    }else if(lenmodel == 1 & lenmarker > 1){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='FP',y='TP',color='marker'))

    }else if(lenmodel > 1 & lenmodel == lenmarker){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='FP',y='TP',color='marker'))

    }else if(lenmodel > 1 & lenmodel != lenmarker){
        p <- ggplot(pp)
        p <- p + geom_line(aes_string(x='FP',y='TP',color='marker')) +
            facet_wrap(~model,ncol = length(unique(pp$model)))
    }
    if (!is.null(color)){
        p <- p+scale_color_manual(values = color)
    }
    p +
        xlab('1 - Specificity')+
        ylab('Sensitivity') +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
        coord_fixed() +
        theme_bw() +
        geom_abline(intercept = 0, slope = 1, color="gray",
                    linetype="dashed", size=0.75)+
        theme(legend.position = 'right',
              panel.spacing = unit(grid.space, "lines"))
    # legend.title = element_blank()
}


