#' Generate Heatmaps
#'
#' Generate heatmaps  of trans-eQTL likelihood: -log10( cummumaltive probability of number of hits for general or specific TF)
#' @param obj an teQTL class object that has been initialized using "init_teQTL"
#' @param gl named character vector of ensembl gene IDs. Names should be human readable gene-names
#' @param tissues character vector of tissues (from the samples data.frame in doi: 10.5281/zendo.838733) specifying those tissues of where the regression should be performed (not yet implimented). If no tissues are specified, the default is all tissues.
#' @param genes gene annotation data.frame from doi: 10.5281/zendo.838733
#' @param exp gene expression matrix from doi: 10.5281/zendo.838733
#' @param samples samples data.frame in doi: 10.5281/zendo.838733
#' @param filePrefix string, filename for the output
#' @param thresh number [0,1], minimum cummulative probability of TF occcurances. default is .9
#' @param dosave boolean, if true, the regression visualization will be saved in the location specified by the filePrefix
#' @param meancor number, minimum correlation to be included in the visualization
#' @return object of class "init_teQTL" extended to contain gene expression correlation between Targets and TFs in obj$cr_TSspec (only specified tissues) and obj$cr_all (all tissues) respectively
#' @export
#'
gen_cor <- function(obj, gl, tissues , genes , exp,samples ,filePrefix='file',thresh=.9,dosave=T,meancor=.1){
  if(sum(names(obj)%in%c('TF_v_targets.cdf_spec','TF_v_targets.cdf_all'))<2){stop('must run gen_heatmap before gen_cor')}
  cr_spec = matrix(0,nrow=length(unique(obj$TF_sel$names)),ncol=length(names(gl)),dimnames=list(unique(obj$TF_sel$names),names(gl)))
  cr_all = matrix(0,nrow=length(unique(obj$TF_sel$names)),ncol=length(names(gl)),dimnames=list(unique(obj$TF_sel$names),names(gl)))
  for(i in 1:nrow(obj$TF_sel)){
    if(obj$TF_sel$cdf_all[i]>thresh | obj$TF_sel$cdf_spec[i]>thresh){
      for(g in strsplit(as.character(obj$TF_sel$gene_assc[i]),',')[[1]]){
        TFi=obj$TF_sel$names[i]
        TFid=as.character(genes$Name)[as.character(genes$Symbol)==TFi]
        gid=gl[names(gl)==g]
        if(length(gid)>1){stop('ambiguous gene names. please confirm that official gene names are unique')}
        TSidx = samples$Tissue %in% tissues

        cr_all[TFi,g] = cor( exp[ TFid , ] , exp[ gid , ] , method='spearman')
        cr_spec[TFi,g] = cor( exp[ TFid , TSidx ] , exp[ gid , TSidx ] , method='spearman')
      }
    }
  }

  if(dosave){
    pdf(paste0(filePrefix,'.correlation.pdf'),height=30,width=20)
    heatmap.2(cr_spec[rowMeans(abs(cr_spec))>meancor,],trace='none',main=paste('spearman correlation\n',paste(tissues,collapse=',')),mar=c(8,8),col=cm.colors)
    heatmap.2(cr_all[rowMeans(abs(cr_all))>meancor,],trace='none',main='spearman correlation',mar=c(8,8),col=cm.colors)
    dev.off()
  }

  obj$cr_TSspec = cr_spec
  obj$cr_all = cr_all

  obj
}
