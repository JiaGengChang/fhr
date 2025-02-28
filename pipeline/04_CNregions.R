library(GenomicRanges)
library(iClusterPlus)
library(cluster)

mmseg = read.delim('./data/Copy Number Estimates_MMRF_CoMMpass_IA22_genome_gatk_cna.seg')
mmseg['Chromosome'] = substr(mmseg$Chromosome,4,6)
mmseg['Chromosome'] = as.integer(ifelse(mmseg$Chromosome=='X',23,ifelse(mmseg$Chromosome=='Y',24,mmseg$Chromosome)))

# optional - convert segment mean to integer copy number
seg_mean_to_cn = function(x){
  ifelse(x > 0.66, 2,
         ifelse(x > 0.38, 1,
                ifelse(x > -0.5, 0,
                       ifelse(x > -1.5, -1, -2))))
}


# there is no data of the germline copy number variations
# leave epsilon parameter as adaptive
# finds 1464 features
cnregions=CNregions(seg=mmseg, adaptive=TRUE,
               frac.overlap=0.5, rmSmallseg=TRUE)

# 1128 by 1464
write.table(cnregions,'./matrices/gene_cn_matrix_CNregions.tsv',sep='\t')

