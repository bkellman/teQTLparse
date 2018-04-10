test<-function(){

  # download GTEx PANDA tissues Rdata from: https://doi.org/10.5281/zenodo.838734
  load('GTEx_PANDA_tissues.RData')

  out=list()
  #"Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_other","Brain_cerebellum","Brain_basal_ganglia","Breast","Lymphoblastoid_cell_line","Fibroblast_cell_line","Colon_sigmoid","Colon_transverse","Gastroesophageal_junction","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Skeletal_muscle","Tibial_nerve","Ovary","Pancreas","Pituitary", "Prostate","Skin","Intestine_terminal_ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_blood"
  tissues = c("Artery_aorta","Artery_coronary","Artery_tibial")


  ###########
  r = read.table('hs_genes.txt',sep='\t',header=T)
  gl = r$To
  names(gl) = r$From

  out[['general_HS_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                                                      minTS=0 , rand=n , onlyCanonical=T,filePrefix='general_HS_regulators')
  ###########
  r = read.table('lysosome_genes.txt',sep='\t',header=T)
  gl = r$To
  names(gl) = r$From

  out[['general_lysosome_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                                                            minTS=0 , rand=n , onlyCanonical=T,filePrefix='general_lysosome_regulators')
  ###########
  r = read.table('cholesterol_biosynthesis.txt',sep='\t',header=T)
  gl = r$To
  names(gl) = r$From

  out[['general_cholesterol_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                                                               minTS=0 , rand=n , onlyCanonical=T,filePrefix='general_cholesterol_regulators')

  ###########
  r = read.table('Cellular_iron_msigdb.txt',sep='\t',header=F)
  colnames(r) = c('From','To','Species','LongName')
  gl = r$To
  names(gl) = r$From

  out[['Cellular_iron_regulators']] = init_teQTL( gl, tissues=c('liver'), edges , genes , netTS , expTS  ,
                                                         minTS=0 , rand=n , onlyCanonical=T,filePrefix='Cellular_iron_regulators')

  ######

  ### gen heatmaps
  for(n in names(out)){
    out[[n]] = gen_heatmap( out[[n]], out[[n]]$gl ,filePrefix=n,thresh=.9,dosave=T)
  }
  ### gen correlation
  for(n in names(out)){
    out[[n]] = gen_cor( out[[n]], out[[n]]$gl, tissues , genes , exp,samples ,filePrefix=n,thresh=.9,dosave=T,meancor=.05)
  }
  ### gen regressions
  for(n in names(out)){
    out[[n]] = gen_regression(out[[n]], out[[n]]$gl, tissues , genes , exp,samples ,filePrefix=n,thresh=.9,dosave=T,min=2)
  }

  ### tissue expss
  for(n in names(out)){
    tissues_expss(obj=out[[n]],exp,genes,samples,filePrefix=n,dosave=T)
  }

  ######
  save(out,file='tfs.rda')
}
