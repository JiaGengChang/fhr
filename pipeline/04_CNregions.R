library(GenomicRanges)
library(iClusterPlus)
library(cluster)
source('./utils/copynumber.R')

mmseg = read.delim('./data/Copy Number Estimates_MMRF_CoMMpass_IA22_genome_gatk_cna.seg')
mmseg['Chromosome'] = substr(mmseg$Chromosome,4,6)
mmseg['Chromosome'] = as.integer(ifelse(mmseg$Chromosome=='X',23,ifelse(mmseg$Chromosome=='Y',24,mmseg$Chromosome)))

# there is no data of the germline copy number variations
# leave epsilon parameter as adaptive
# finds 1464 features
cn=CNregions(seg=mmseg, adaptive=TRUE,
               frac.overlap=0.5, rmSmallseg=TRUE)

# 1128 by 1464
# write.table(cnregions,'./matrices/gene_cn_matrix_CNregions.tsv',sep='\t')

# remove subsequent visits

# cn = read.delim('./matrices/segment_cn_matrix.tsv')
names_cn = substr(rownames(cn),1,9)
chroms_cn = sapply(colnames(cn), function(...) as.integer(strsplit(strsplit(...,"\\.")[[1]][[1]],"chr")[[1]][[2]]))
cn = cn[!duplicated(names_cn), colnames(cn)[order(chroms_cn)]]
rownames(cn) = substr(rownames(cn),1,9)

# locations of values < -9
# missing.ind = arrayInd(which(cn < -9),dim(cn))
# missing.row = unique(missing.ind[,1])

# Option 1: drop individuals with missing CN data
# cn = cn[-missing.row,]
# Option 2: clip the -30 values to -8.79
# cn[cn < -9] = min(cn[cn > -9])

# Convert segment mean to integer CN
cn.out = apply(apply(cn, 2, seg_mean_to_cn), 2, as.integer)
rownames(cn.out) = rownames(cn)

# write integer level segment cn data
write.table(cn.out,'./matrices/segment_cn_matrix.tsv',sep='\t',quote=F)
