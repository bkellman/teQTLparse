
library(teQTLparse)
library(openxlsx)
library(ggplot2)



# download GTEx PANDA tissues Rdata from: https://doi.org/10.5281/zenodo.838734
load('GTEx_PANDA_tissues.RData')

#"Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Artery_aorta","Artery_coronary","Artery_tibial","Brain_other","Brain_cerebellum","Brain_basal_ganglia","Breast","Lymphoblastoid_cell_line","Fibroblast_cell_line","Colon_sigmoid","Colon_transverse","Gastroesophageal_junction","Esophagus_mucosa","Esophagus_muscularis","Heart_atrial_appendage","Heart_left_ventricle","Kidney_cortex","Liver","Lung","Minor_salivary_gland","Skeletal_muscle","Tibial_nerve","Ovary","Pancreas","Pituitary", "Prostate","Skin","Intestine_terminal_ileum","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_blood"
# tissues = c("Artery_aorta","Artery_coronary","Artery_tibial")
tissues = c("Adipose_subcutaneous","Adipose_visceral","Adrenal_gland","Liver","Pancreas","Pituitary","Spleen")

###### start benchmark
out_bench=list()

r = read.table('teQTLparse/examples/Cholesterol/cholesterol_biosynthesis.txt',sep='\t',header=T)
gl = r$To
names(gl) = r$From

## initialize at different sampling levels
if(!file.exists('teQTLparse/examples/benchmark.rda')){
	for(iter in c(50,100,500,1000,1500,2000)){
	  out_bench[[paste0('general_cholesterol_regulators.n_',iter)]] = init_teQTL( gl, tissues, edges , genes , netTS , expTS  ,
	                          minTS=0 , rand=iter , onlyCanonical=T,filePrefix='general_cholesterol_regulators')
	}
	save(out_bench,file='teQTLparse/examples/benchmark.rda')
}else{
	load('teQTLparse/examples/benchmark.rda')
}

## load benchmark
library(enrichR)

dbs <- listEnrichrDbs()
dbi=c('ChEA_2016','TRANSFAC_and_JASPAR_PWMs','ARCHS4_TFs_Coexp','Enrichr_Submissions_TF-Gene_Coocurrence','ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X','ENCODE_TF_ChIP-seq_2015')

enriched <- enrichr(names(gl), dbi)
enriched <- do.call(rbind,enriched)
enriched$study = unlist(lapply(strsplit(rownames(enriched),'\\.'),function(x) x[1] ))
enriched$TF = unlist(lapply(strsplit(enriched$Term,'_| '),function(x) x[1] ))

# agg = aggregate(enriched$Adjusted.P.value,by=list(enriched$TF),median)
# levels(enriched$TF) = factor( as.character(enriched$TF) , levels=as.character(agg[,1])[order(agg[,2])])
# ggplot(enriched[enriched$TF %in% agg[agg[,2]<.1,1],],aes(x=TF,y=-log(Adjusted.P.value)))+geom_boxplot()+
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1))


## literature
lit = c('PPARA','PPARG','PPARB','PPARD','CEBPA','CEBPB','CEBPD','CEBPG','CEBPZ','MLXIPL','SREBP1','SREBP2','SREBP1c','DLK1','ZNF423')

##############################
####### where known TFs show up intrinsically
# teQTL discovered
out_bench[[1]]$TF_sel[grepl(paste(lit,collapse='|'),out_bench[[1]]$TF_sel$names),]
# enrichR discovered
 enriched[grepl(paste(lit,collapse='|'),enriched$TF),]
###################################

###################################
# error rate for each method
library(ROCR)


# use F measure because True Negative (genes that are not TFs) is unknown and/or poorly defined

# threshs=c(.4,.3,.2,.1,.05)

all = unique(c(unique(unlist(lapply(out_bench,function(x) as.character(x$TF_sel$names)))),lit))

# teQTL
perf=list()
for(i in 1:length(out_bench)){
	inam = gsub('^n_','',strsplit(names(out_bench)[i],'\\.')[[1]])
	# get p-values for the selection
	p = p.adjust(1-out_bench[[i]]$TF_sel$cdf_spec,'fdr')
	# get TFs corresponding to those p-values
	TF = as.character(out_bench[[i]]$TF_sel$names)
	# set scores to 1-p for observed TFs, assign 0 to missed TFs
	predi = sapply(all,function(x) ifelse(x%in%TF,1-p[TF==x],0))
	# set labels to 1=T, 0=F
	labels = ifelse( all%in%lit,1,0)
	pred <- prediction( predi, labels)
	## precision/recall curve (x-axis: recall, y-axis: precision)
	perf1 = performance(pred, "prec", "rec")
	perf2 = performance(pred, 'f')@y.values[[1]]
	perf3 = performance(pred, 'phi')@y.values[[1]]
	perf4 = performance(pred,"tpr","fpr")
	perf5 = performance(pred,"auc")@y.values[[1]]


	perf[[names(out_bench)[i]]] = data.frame(precision=perf1@x.values[[1]],recall=perf1@y.values[[1]],
		alpha=perf1@alpha.values[[1]],f=perf2,mcc=perf3,auc=perf5,
		tpr=perf4@x.values[[1]],fpr=perf4@y.values[[1]],
		method=paste('teQTL',inam[2],sep='_') )
}

perf_n = do.call(rbind,perf)
perf_n$method = factor(as.character(perf_n$method),levels=unique(perf_n$method))
 ggplot(perf_n,aes(color=method))+
 	geom_line(aes(x=recall,y=precision,size=f))
 ggplot(perf_n,aes(color=method))+
 	geom_line(aes(x=tpr,y=fpr,size=mcc))
  ggplot(perf_n,aes(color=method))+
 	geom_line(aes(x=f,y=mcc))

# other
perf_compare=list()
perf_compare$teQTL_1000 = perf$general_cholesterol_regulators.n_1000
for(i in unique(enriched$study)){
	# get p-values for the selection
	p = enriched$Adjusted.P.value[ enriched$study==i ]
	# get TFs corresponding to those p-values
	TF = enriched$TF[ enriched$study==i ]
	# set scores to 1-p for observed TFs, assign 0 to missed TFs
	predi = sapply(all,function(x) ifelse(x%in%TF,1-p[TF==x],0))
	# set labels to 1=T, 0=F
	labels = ifelse( all%in%lit,1,0)
	pred <- prediction( predi, labels)
	## precision/recall curve (x-axis: recall, y-axis: precision)
	perf1 = performance(pred, "prec", "rec")
	perf2 = performance(pred, 'f')@y.values[[1]]
	perf3 = performance(pred, 'phi')@y.values[[1]]
	perf4 = performance(pred,"tpr","fpr")
	perf5 = performance(pred,"auc")@y.values[[1]]

	perf_compare[[i]] = data.frame(precision=perf1@x.values[[1]],recall=perf1@y.values[[1]],
		alpha=perf1@alpha.values[[1]],f=perf2,mcc=perf3,auc=perf5,
		tpr=perf4@x.values[[1]],fpr=perf4@y.values[[1]],
		method=i)
}

perf_compare_df = do.call(rbind,perf_compare)
 ggplot(perf_compare_df,aes(color=method))+
 	geom_line(aes(x=recall,y=precision,size=f))
 ggplot(perf_compare_df,aes(color=method))+
 	geom_line(aes(x=tpr,y=fpr,size=mcc))
  ggplot(perf_compare_df,aes(color=method))+
 	geom_line(aes(x=f,y=mcc))

###################################

##################################
# teQTL v other , highlight if in CholTF

