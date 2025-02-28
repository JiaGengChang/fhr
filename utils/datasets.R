source('./utils/copynumber.R')

# function to load all data
load_data = function(){
  # clin = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
  # clin = clin[clin$risk != -1,]
  # clin$risk = factor(clin$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))
  # write.table(clin['risk'],'./matrices/risk_labels.tsv',sep='\t',quote=F)
  data.label = read.table('./matrices/risk_labels.tsv',sep='\t')
  
  data.cn = read.table('./matrices/broad_cn_matrix_integer.tsv',sep='\t')
  data.cn = as.matrix(data.cn)
  
  data.rna = read.delim('./matrices/gene_exp_logcpm_proteincoding_signif772.tsv',sep='\t')
  data.rna = as.matrix(data.rna)
  
  data.mut = read.delim('./matrices/gene_mut_matrix_signif33.tsv',sep='\t')
  data.mut = as.matrix(data.mut)
  
  # gene fusion data
  data.fus = read.delim('./matrices/gene_fusion_matrix_signif6.tsv',sep='\t')
  data.fus = as.matrix(data.fus)
  
  assign("data.cn",data.cn,envir=.GlobalEnv)
  assign("data.rna",data.rna,envir=.GlobalEnv)
  assign("data.mut",data.mut,envir=.GlobalEnv)
  assign("data.fus",data.fus,envir=.GlobalEnv)
  assign("data.label",data.label,envir=.GlobalEnv)
  
  return(0)
}
