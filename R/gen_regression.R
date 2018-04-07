#' Generate Regressions
#'
#' Takes predicted trans-eQTLs and fits regressions to describe the relation between predicted TFs and the target
#' @param obj an teQTL class object that has been initialized using "init_teQTL"
#' @param gl named character vector of ensembl gene IDs. Names should be human readable gene-names
#' @param tissues character vector of tissues (from the samples data.frame in doi: 10.5281/zendo.838733) specifying those tissues of where the regression should be performed (not yet implimented). If no tissues are specified, the default is all tissues.
#' @param genes gene annotation data.frame from doi: 10.5281/zendo.838733
#' @param exp gene expression matrix from doi: 10.5281/zendo.838733
#' @param samples samples data.frame in doi: 10.5281/zendo.838733
#' @param filePrefix string, filename for the output
#' @param dosave boolean, if true, the regression visualization will be saved in the location specified by the filePrefix
#' @param min number, minimum -log10(cummulative probability of TF occcurances), default is 2
#' @return object of class "init_teQTL" extended to contain regression objects for each target stored in obj$regressions and visualizations of those regressions stored in obj$regression_fig
#' @export
#'
gen_regression <- function(obj, gl, tissues=NULL , genes , exp,samples ,filePrefix='file',dosave=T,min=2){
  if(sum(names(obj)%in%c('TF_v_targets.cdf_spec','TF_v_targets.cdf_all'))<2){stop('must run gen_heatmap before gen_regression')}
  assc = obj$TF_v_targets.cdf_spec
  obj$regression = lapply(names(gl),function(g){
    print(g)
    TFi = rownames(assc)[assc[,g]>min]
    #if(length(TFi)>10){TFi = (rownames(assc)[order(-assc[,g])])[1:10]}
    if(length(TFi)==0){return(NA)}
    TFid=as.character(genes$Name)[as.character(genes$Symbol) %in% TFi]
    gid=as.character( gl[names(gl)==g] )
    TSidx = samples$Tissue %in% tissues
    if(is.null(tissues)){
      dat=as.data.frame(scale(log(t(exp[ c(TFid,gid) , ])+1)))
    }else{
      dat=as.data.frame(scale(log(t(exp[ c(TFid,gid) , TSidx ])+1)))
    }
    # univariate screen
    uni_screen = p.adjust( sapply( TFid , function(tfi){ anova(glm(as.formula(paste(gid,'~',tfi)),data=dat),test='LRT')[2,5] }) , 'fdr')
    if(length(TFi)>3){
      TFid = TFid[uni_screen<=max(1e-20,median(uni_screen))]
      # cbx screen
      cbx_screen =  p.adjust( apply( cbx<-combn(TFid,2),2 , function(tfi){
        anova(glm(as.formula(paste(gid,'~',tfi[1],'+',tfi[2])),data=dat),
              glm(as.formula(paste(gid,'~',tfi[1],'*',tfi[2])),data=dat),test='LRT')[2,5]
      }) , 'fdr')
      if(dim(cbx)[2]>4){ cbx = cbx[,cbx_screen<=max(1e-20,median(cbx_screen))] }
      # multivariate model
      ret=stepAIC( glm( as.formula( paste(gid,'~',paste(TFid,collapse='+'))),data=dat) , direction='both',
                   scope=list(lower=~1,upper=as.formula(paste('~',paste(c(TFid,apply(cbx,2,paste,collapse=':')),collapse='+')))),k=length(TFi))
    }else{
      ret=stepAIC( glm( as.formula( paste(gid,'~',paste(TFid,collapse='+'))),data=dat) , direction='both',
                   scope=list(lower=~1,upper=~.^2),k=length(TFi))
    }
    if(dim(coef(summary(ret)))[1]==1){return(NA)}
    ret
  })
  names(obj$regression) = names(gl)
  obj$regression =  obj$regression[unlist(lapply(obj$regression,function(x) class(x)[1]!='logical'))]

  g = regression_fig(obj,genes)

  obj$regression_fig = g

  if(dosave){
    ggsave(g,filename=paste0(filePrefix,'regressions.pdf'),height=length(gl),width=length(gl))
  }

  obj
}
