
# function to load all data
load_all_data = function(){
  
  fhr = read.delim('./annotations/fhr-annotations-raw.2Mar25.tsv',row.names = 1)
  fhr = fhr[fhr$risk != -1,]
  fhr$risk = factor(fhr$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))
  all.label = fhr['risk']
  
  # arm level integer cn
  all.cn.arm = read.table('./matrices/broad_cn_matrix_integer.tsv',sep='\t')
  all.cn.arm = as.matrix(all.cn.arm)
  
  # segment level integer cn
  all.cn.seg = read.delim('./matrices/segment_cn_matrix_uncorrelated.tsv',sep='\t')
  all.cn.seg = as.matrix(all.cn.seg)
  
  all.rna = read.delim('./matrices/gene_exp_logcpm_proteincoding.tsv',sep='\t')
  all.rna = as.matrix(all.rna)
  annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',sep='\t',row.names=1) 
  # rename ENSG IDs to gene names where possible
  gene.name.rna = annot[colnames(all.rna),'Gene.name']
  # There are 8 pairs duplicated gene names
  duplicated.names = gene.name.rna[duplicated(gene.name.rna) & gene.name.rna!=""]
  duplicated = gene.name.rna %in% duplicated.names
  missing = gene.name.rna==""
  gene.name.id.rna = ifelse(missing | duplicated, colnames(all.rna),gene.name.rna)
  colnames(all.rna) = gene.name.id.rna
  
  all.mut = read.delim('./matrices/gene_mut_matrix.tsv',sep='\t',row.names=1)
  all.mut = as.matrix(all.mut)
  
  # canonical IG translocations
  all.ig = read.delim('./matrices/canonical_ig_translocations.tsv',sep='\t')
  all.ig = as.matrix(all.ig)
  
  assign("all.label", all.label, envir=.GlobalEnv)
  assign("all.rna", all.rna, envir=.GlobalEnv)
  assign("all.cn.arm", all.cn.arm, envir=.GlobalEnv)
  assign("all.cn.seg", all.cn.seg, envir=.GlobalEnv)
  assign("all.ig", all.ig, envir=.GlobalEnv)
  assign("all.mut", all.mut, envir=.GlobalEnv)
  
  return(invisible(list(label=all.label,rna=all.rna,cn.arm=all.cn.arm,cn.seg=all.cn.seg,ig=all.ig,mut=all.mut)))
}

# load data for train and validation sets
load_dev_data = function(){
  
  # load all data
  if (!(exists("dev.ig") & exists("dev.mut") & exists("dev.rna") & exists("dev.cn.arm") & exists("dev.cn.seg") & exists("dev.label"))) {
    load_all_data()
  }
  
  # public IDs of patients with all data,
  # excluding 150 patients in test set
  dev.id = read.csv('./data/splits/dev.txt',sep=',')$PUBLIC_ID
  
  dev.label = all.label[dev.id,,drop=FALSE]
  dev.rna = all.rna[dev.id,]
  dev.cn.arm = all.cn.arm[dev.id,]
  dev.cn.seg = all.cn.seg[dev.id,]
  dev.mut = all.mut[dev.id, colSums(all.mut[dev.id,]) > 0]
  dev.ig = all.ig[dev.id,]
  
  assign("dev.label", dev.label, envir=.GlobalEnv)
  assign("dev.rna", dev.rna, envir=.GlobalEnv)
  assign("dev.cn.arm", dev.cn.arm, envir=.GlobalEnv)
  assign("dev.cn.seg", dev.cn.seg, envir=.GlobalEnv)
  assign("dev.mut", dev.mut, envir=.GlobalEnv)
  assign("dev.ig", dev.ig, envir=.GlobalEnv)
  
  return(invisible(list(label=dev.label,rna=dev.rna,cn.arm=dev.cn.arm,cn.seg=dev.cn.seg,mut=dev.mut,ig=dev.ig)))
}

load_split_data = function(shuffle=10,fold=5){
  # load data for a given train test split
  train.id = read.csv(sprintf('./data/splits/%d/train_%d.txt',shuffle,fold))$PUBLIC_ID
  valid.id = read.csv(sprintf('./data/splits/%d/valid_%d.txt',shuffle,fold))$PUBLIC_ID
  
  # load data with dev information
  if (!(exists("dev.ig") & exists("dev.mut") & exists("dev.rna") & exists("dev.cn.arm") & exists("dev.cn.seg") & exists("dev.label"))) {
    load_dev_data()
  }
  
  # split into train and valid subsets
  assign("train.id", train.id, envir=.GlobalEnv)
  assign("valid.id", valid.id, envir=.GlobalEnv)
  assign("train.label", dev.label[train.id,,drop=F], envir=.GlobalEnv)
  assign("valid.label", dev.label[valid.id,,drop=F], envir=.GlobalEnv)
  assign("train.rna", dev.rna[train.id,], envir=.GlobalEnv)
  assign("valid.rna", dev.rna[valid.id,], envir=.GlobalEnv)
  assign("train.cn.arm", dev.cn.arm[train.id,], envir=.GlobalEnv)
  assign("valid.cn.arm", dev.cn.arm[valid.id,], envir=.GlobalEnv)
  assign("train.cn.seg", dev.cn.seg[train.id,], envir=.GlobalEnv)
  assign("valid.cn.seg", dev.cn.seg[valid.id,], envir=.GlobalEnv)
  assign("train.ig", dev.ig[train.id,], envir=.GlobalEnv)
  assign("valid.ig", dev.ig[valid.id,], envir=.GlobalEnv)
  assign("train.mut", dev.mut[train.id,], envir=.GlobalEnv)
  assign("valid.mut", dev.mut[valid.id,], envir=.GlobalEnv)
  
  return(invisible(list(train=train.id,valid=valid.id)))
}