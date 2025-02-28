# input './matrices/gene_cn_matrix.tsv.gz'
# output './matrices/gene_cn_integer_matrix.tsv.gz'

seg = read.delim('~/Documents/fhr/data/Copy Number Estimates_MMRF_CoMMpass_IA22_genome_gatk_cna.seg')
plot = hist(seg[seg$Num_Probes>1000,'Segment_Mean'],
     breaks=10000,
     xlim=c(-2,2),
     main="Length >= 1000",
     xlab="Segment mean",)
lines(x=c(-1.5,-1.5),y=c(0,7000),col='red')
lines(x=c(-0.5,-0.5),y=c(0,7000),col='red')
lines(x=c(0.38,0.38),y=c(0,7000),col='red')
lines(x=c(0.66,0.66),y=c(0,7000),col='red')

data = read.delim('./matrices/gene_cn_matrix.tsv.gz',row.names = 1)
data = data[, 2:ncol(data)] # drop sample column

# convert to integer CN using
# segment mean = log2(copy number/2)
# cn = data
# cn[data = log2(4/2)] = 4
# cn[data = log2(3/2) ] = 3
# cn[data = log2(2/2)] = 2 # normal
# cn[data = log2(1/2)] = 1
# cn[data < log2(1/2)] = 0

seg_mean_to_cn = function(x){
  ifelse(x > 0.66, 2,
         ifelse(x > 0.38, 1,
                ifelse(x > -0.5, 0,
                       ifelse(x > -1.5, -1, -2))))
}

cn = apply(apply(data, 2, seg_mean_to_cn), 2, as.integer)
rownames(cn) = rownames(data)

write.table(cn,'./matrices/gene_cn_integer_matrix.tsv.gz',sep='\t',quote=F)
