# function to load all data
load_data = function(){
  source('./utils/copynumber.R')
  annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',sep='\t',row.names=1)
  
  fhr = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
  fhr = fhr[fhr$risk != -1,]
  fhr$risk = factor(fhr$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))
  data.label = fhr['risk']
  
  data.cn = read.table('./matrices/broad_cn_matrix_integer.tsv',sep='\t')
  data.cn = as.matrix(data.cn)
  
  data.rna = read.delim('./matrices/gene_exp_logcpm_proteincoding_signif772.tsv',sep='\t')
  data.rna = as.matrix(data.rna)
  gene.name.rna = annot[colnames(data.rna),'Gene.name']
  colnames(data.rna) = ifelse(gene.name.rna=="",colnames(data.rna),gene.name.rna)
  
  data.mut = read.delim('./matrices/gene_mut_matrix_signif33.tsv',sep='\t')
  data.mut[data.mut > 1] = 1L
  data.mut = as.matrix(data.mut)
  colnames(data.mut) = annot[colnames(data.mut),'Gene.name']
  
  # canonical IG translocations - 3 significant ones
  data.ig = read.delim('./matrices/canonical_ig_translocations_signif3.tsv',sep='\t')
  data.ig = as.matrix(data.ig)
  
  # gene fusion data
  # data.fus = read.delim('./matrices/gene_fusion_matrix_signif6.tsv',sep='\t')
  # data.fus = as.matrix(data.fus)
  
  assign("data.cn",data.cn,envir=.GlobalEnv)
  assign("data.rna",data.rna,envir=.GlobalEnv)
  assign("data.mut",data.mut,envir=.GlobalEnv)
  assign("data.ig",data.ig,envir=.GlobalEnv)
  assign("data.label",data.label,envir=.GlobalEnv)
  
  return(0)
}

load_common_data = function(){
  load_data()
  common.id = Reduce(intersect,list(rownames(data.cn),rownames(data.rna),rownames(data.mut),rownames(data.ig)),rownames(data.label))
  
  common.data.cn = data.cn[common.id,,drop=FALSE]
  common.data.rna = data.rna[common.id,,drop=FALSE]
  common.data.ig = data.ig[common.id,,drop=FALSE]
  common.data.label = data.label[common.id,,drop=FALSE]
  common.data.mut = data.mut[common.id, colSums(data.mut[common.id,]) > 0, drop=FALSE] # filter genes without anymore mutations
  
  assign("common.data.cn", common.data.cn, envir=.GlobalEnv)
  assign("common.data.rna", common.data.rna, envir=.GlobalEnv)
  assign("common.data.ig", common.data.ig, envir=.GlobalEnv)
  assign("common.data.label", common.data.label, envir=.GlobalEnv)
  assign("common.data.mut", common.data.mut, envir=.GlobalEnv)
}
