test<-function(){

  load('GTEx_PANDA_tissues.RData')

  out=list()
  #"Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_other","Brain_cerebellum","Brain_basal_ganglia","Breast","Lymphoblastoid_cell_line","Fibroblast_cell_line","Colon_sigmoid","Colon_transverse","Gastroesophageal_junction","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Skeletal_muscle","Tibial_nerve","Ovary","Pancreas","Pituitary", "Prostate","Skin","Intestine_terminal_ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_blood"
  tissues = c("Artery_aorta","Artery_coronary","Artery_tibial")

  gl = c('ENSG00000166507','ENSG00000122863','ENSG00000147119','ENSG00000137573','ENSG00000070614',
         'ENSG00000138604','ENSG00000145681','ENSG00000038427','ENSG00000114378','ENSG00000068001')
  names(gl) = c('NDST2','CHST3','CHST7','SULF1','NDST1','GLCE','HAPLN1','VCAN','HYAL1','HYAL2')

  out[['ShearSpecific_glycocalyx_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                                                                    minTS=0 , rand=n , onlyCanonical=T,filePrefix='ShearSpecific_glycocalyx_regulators')

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
