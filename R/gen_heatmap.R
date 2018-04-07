#' Generate Heatmaps
#'
#' Generate heatmaps  of trans-eQTL likelihood: -log10( cummumaltive probability of number of hits for general or specific TF)
#' @param obj an teQTL class object that has been initialized using "init_teQTL"
#' @param gl named character vector of ensembl gene IDs. Names should be human readable gene-names
#' @param filePrefix string, filename for the output
#' @param thresh number [0,1], minimum cummulative probability of TF occcurances. default is .9
#' @param dosave boolean, if true, the regression visualization will be saved in the location specified by the filePrefix
#' @param min number, minimum -log10(cummulative probability of TF occcurances), default is 1
#' @return object of class "init_teQTL" extended to contain TF specific and general aggregated -log10 of cummulative probabilities for TF hit counts for each target stored in obj$TF_v_targets.cdf_all and obj$TF_v_targets.cdf_spec respectively
#' @export
#'
gen_heatmap <- function( obj, gl ,filePrefix='file',thresh=.9,dosave=T,minNeglog10Perc=1){
  tmp1 = matrix(0,nrow=length(unique(obj$TF_sel$names)),ncol=length(names(gl)),dimnames=list(unique(obj$TF_sel$names),names(gl)))
  tmp2 = matrix(0,nrow=length(unique(obj$TF_sel$names)),ncol=length(names(gl)),dimnames=list(unique(obj$TF_sel$names),names(gl)))
  for(i in 1:nrow(obj$TF_sel)){
    if(obj$TF_sel$cdf_all[i]>thresh){
      tmp1[obj$TF_sel$names[i] ,  strsplit(as.character(obj$TF_sel$gene_assc[i]),',')[[1]] ] =  -log((1-obj$TF_sel$cdf_all[i])+1e-3,10)
    }
    if(obj$TF_sel$cdf_spec[i]>thresh){
      tmp2[obj$TF_sel$names[i] ,  strsplit(as.character(obj$TF_sel$gene_assc[i]),',')[[1]] ] =  -log((1-obj$TF_sel$cdf_spec[i])+1e-3,10)
    }
  }

  if(dosave){
    pdf(paste0(filePrefix,'.percentile.pdf'),height=30,width=20)
    heatmap.2(tmp1[rowSums(tmp1)>minNeglog10Perc,],trace='none',
              main='-log10(1-percentile(# of targets hit))\nrelative to ALL TFs with random targets',mar=c(8,8))
    heatmap.2(tmp2[rowSums(tmp2)>minNeglog10Perc,],trace='none',
              main='-log10(1-percentile(# of targets hit))\nrelative to THIS TFs with random targets',mar=c(8,8))
    dev.off()
  }

  obj$TF_v_targets.cdf_all = tmp1
  obj$TF_v_targets.cdf_spec = tmp2

  obj
}
