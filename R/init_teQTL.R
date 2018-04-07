#' Initialize the teQTL object
#'
#' Initilizes the fundamental teQTL object
#' @param gl named character vector of ensembl gene IDs. Names should be human readable gene-names
#' @param tissues character vector of tissues (from the samples data.frame in doi: 10.5281/zendo.838733) specifying those tissues of where the regression should be performed (not yet implimented). If no tissues are specified, the default is all tissues.
#' @param edges edges data.frame containing predicted trans-eQTL interactions from doi: 10.5281/zendo.838733
#' @param genes gene annotation data.frame from doi: 10.5281/zendo.838733
#' @param netTS netTS data.frame in doi: 10.5281/zendo.838733
#' @param expTS expTS data.frame in doi: 10.5281/zendo.838733
#' @param minTS integer, minimum number of specified tissues which must contain a teQTL for it to be considered. default is 0 which disables this parameter because it is too limiting.
#' @param rand inteber, number of random genelists to examine for the background
#' @param filePrefix string, filename for the output
#' @param thresh number [0,1], minimum cummulative probability of TF occcurances. default is .9
#' @param dosave boolean, if true, the regression visualization will be saved in the location specified by the filePrefix
#' @return object of class "init_teQTL" extended to contain TF specific and general cummulative probabilities for TF hit counts for each target stored in obj$TF_sel
#' @export
#'
init_teQTL <- function( gl, tissues, edges , genes , netTS , expTS  , minTS=0 , rand=50 , onlyCanonical=T,filePrefix='file',thresh=.9,dosave=T){
  # filter possible TFs by expression within those tissues
  keep_gen = rowSums( expTS[,colnames(expTS) %in% tissues] ) >= minTS
  g = genes[keep_gen,]
  g = genes

  # filter possible TFs by tissue of activity AND canonicalTF AND TS-expression
  keep_net = rowSums( netTS[,colnames(netTS) %in% tissues] ) >= minTS
  if(onlyCanonical){ keep_net = keep_net & edges$Prior==1 }
  keep_net = keep_net & edges$TF %in% g$Symbol
  # get regulators of gene list
  e = edges[keep_net & edges$Gene %in% gl,]
  # get regulators of n random gene lists
  e_rand = lapply(1:rand, function(x){
    gl_i = sample( g$Name , length(gl) ) # generate geen
    edges[keep_net & edges$Gene %in% gl_i,]
  } )

  # get occurances of regulators
  occ = table( droplevels(e)$TF)

  # get percentile for occurances of any regulator
  ecdf_all = ecdf( unlist(lapply(e_rand,function(x) table(droplevels(x)$TF))) )
  cdf_all = ecdf_all( occ )
  # get percentile for occurances of a specific regulator
  cdf_spec = sapply(names(occ),function(o_i){
    ecdf_spec = ecdf( unlist(lapply(e_rand,function(x) sum(x$TF==o_i))) )
    ecdf_spec( occ[o_i] )
  })
  gene_assc = sapply( names(occ) , function(o_i){
    paste( names(gl)[gl%in%e$Gene[e$TF==o_i]] , collapse=',' )
  })
  obj=list()
  obj$gl = gl
  obj$TF_sel = data.frame(names=names(occ),occurance=as.numeric(occ),cdf_all=cdf_all,cdf_spec=cdf_spec,gene_assc=gene_assc)
  if(dosave){ write.xlsx(obj$TF_sel,file=paste0(filePrefix,'.xlsx')) }
  attr(obj,'class')='teQTL'
  obj
}
