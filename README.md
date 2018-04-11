# teQTLparse - trans-eQTL parse
A package which parses and ranks predicted trans-eQTLs predicted in [Sonawane 2017](https://doi.org/10.1016/j.celrep.2017.10.001). teQTLparse compares the number of times each predicted TF maps to a given gene list. It compares that number of hits to the expected hits of an arbitrary TF in a random gene list and the expected hits of that specific TF in a random gene list. This comparison is made simply using a cumulative distribution describing the results of multiple random gene list queries.

## Further Explanation
Please see the slideshow: RegulatorPredictor.pptx

# Example: Cholesterol, Lysosome, and Iron Metabolism KEGG Pathway Regulators

Input and Output for these examples can be found in the "examples" folder. Each KEGG pathway is formated in a ".txt" file: col1 is gene names, col2 is Ensembl gene ids. The init_teQTL function reads the number of times genes in the given gene list are predicted to be regulated and the identity of the predicted TF. n-random gene lists, with the same number of genes as the given gene list, are also queried to create a background frequency for any TF as well as specific TFs; the frequency with which they are found to regulate genes in a random gene list. The percentile of predicted regulation ("occurrences") for each TF is recorded and expanded into output files:
- <filename>.xlsx contains the core output in five columns: 
  - names: predicted TF
  - occurrence: number of times this TF was predicted to regulate a gene in the gene list
  - cdf_all: percentile of this TFs occurrences compared to the occurrence frequency for ANY TF
  - cdf_spec: percentile of this TFs occurrences compared to the occurrence frequency for THIS TF
  - gene_assc: genes from the gene list predicted to be regulated by this TF
- <filename>.percentile.pdf contains the xlsx data reformated as a heatmap to help highlight groups of coregulated and co-regulating genes and TFs respectively
- <filename>.correlation.pdf contains the Spearman correlation between gene list and TF genes from the GTEx expression data.
- <filename>.regressions.pdf is the first pass as most specifically describing the correlation between each target (generalist) and regulators (TFs) using multiple regression (univariate selection followed by backward model selection with interactions). These results are from a temporary implementation, later iteration with use the [RegressionPipeline](https://github.com/LewisLabUCSD/RegressionModelPipeline)

```R
# download GTEx PANDA tissues Rdata from: https://doi.org/10.5281/zenodo.838734
load('GTEx_PANDA_tissues.RData')

out=list()
#"Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_other","Brain_cerebellum","Brain_basal_ganglia","Breast","Lymphoblastoid_cell_line","Fibroblast_cell_line","Colon_sigmoid","Colon_transverse","Gastroesophageal_junction","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Skeletal_muscle","Tibial_nerve","Ovary","Pancreas","Pituitary", "Prostate","Skin","Intestine_terminal_ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_blood"
tissues = c("Artery_aorta","Artery_coronary","Artery_tibial")


########### Lysosome
r = read.table('lysosome_genes.txt',sep='\t',header=T)
gl = r$To
names(gl) = r$From

out[['general_lysosome_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                        minTS=0 , rand=n , onlyCanonical=T,filePrefix='general_lysosome_regulators')
########### Cholesterol
r = read.table('cholesterol_biosynthesis.txt',sep='\t',header=T)
gl = r$To
names(gl) = r$From

out[['general_cholesterol_regulators']] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
                        minTS=0 , rand=n , onlyCanonical=T,filePrefix='general_cholesterol_regulators')

########### Iron Metabolism
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
  out[[n]] = gen_cor( out[[n]], out[[n]]$gl, tissues , genes , exp,samples ,
                        filePrefix=n,thresh=.9,dosave=T,meancor=.05)
}
### gen regressions
for(n in names(out)){
  out[[n]] = gen_regression(out[[n]], out[[n]]$gl, tissues , genes , exp,samples ,
                        filePrefix=n,thresh=.9,dosave=T,min=2)
}
### save
save(out,file='tfs.rda')
```
