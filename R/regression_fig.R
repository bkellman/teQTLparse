#' Generate Regression Figure
#'
#' Takes regression objects and visualizes them
#' @param obj regression object output by $gen_regression$ function
#' @param genes gene annotation data.frame from doi: 10.5281/zendo.838733
#' @return grid.arrange output of a list of ggplot grobs
#' @export
#'
regression_fig<-function(obj,genes){
  grid.arrange(grobs=lapply(names(obj$regression),function(xn){
    print(xn)
    x = obj$regression[[xn]]
    dat = as.data.frame(cbind(coef(summary(x)),confint(x)))[-1,]
    colnames(dat) = c('beta','stdErr','t','pr_t','CI_low','CI_high')
    dat$var = rownames(dat)
    dat$var = sapply(dat$var,function(s) paste(genes$Symbol[as.character(genes$Name)%in%strsplit(s,':')[[1]]],collapse=':') )
    ggplot(data=dat,aes(x=var,y=beta,col=-log(pr_t,10)))+ geom_pointrange(aes(ymin = CI_low, ymax = CI_high))+
      geom_hline(yintercept=0)+ggtitle(xn)+theme_classic()+ #scale_color_gradientn(colours = heat.colors(100))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }))
}
