library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')
source('./utils/copynumber.R')
source('./utils/common_fhr_data.R')

# option 1 - use segment mean
data = read.delim('./matrices/gene_cn_matrix.tsv.gz',row.names = 1)
data = data[, 2:ncol(data)] # drop sample column

# option 2 - use integer CNA -2, -1, 0, +1, +2
data = read.delim('./matrices/gene_cn_integer_matrix.tsv.gz')

fhr.0 = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
# option 1 - compare FHR, GHR, and SR. n=884
fhr = fhr.0[(fhr.0['risk'] != -1),]
# option 2 - compare FHR and SR. n=621 
fhr = fhr.0[(fhr.0['risk'] != -1) & (fhr.0['risk'] != 1),]

# get common ids, data, and labels
subset_to_labelled(data, fhr)

# Option 1 - Perform ANOVA on segment mean
p_values = rep(1,ncol(common_data))
for(i in 1:ncol(common_data)){
  # for SR - GHR - FHR comparison
  aov_result = aov(common_data[,i] ~ common_fhr$risk)
  p_values[i] <- summary(aov_result)[[1]]["Pr(>F)"][[1]][1]
  # for SR - FHR only comparison
  # p_values[i] = t.test(common_data[common_fhr$risk == 'SR',i], common_data[common_fhr$risk == 'FHR',i])$p.value
}

# Option 2 - Chi-square test on integer CN categories and risk categories
# not work for FHR-SR 
p_values = apply(common_data, 2, function(...) chisq.test(table(..., common_fhr$risk))$p.value)

# adjust for FDR
p_values_BH = p.adjust(p_values, method='BH')
length(p_values_BH[p_values_BH < 0.05])

# adjust for FWER
p_values_bf = p.adjust(p_values, method='bonferroni')
length(p_values_bf[p_values_bf < 0.05])

# pick genes to use for heatmap
genes = colnames(common_data)[p_values_bf < 0.05]
genes = colnames(common_data)[order(p_values_bf)[1:500]]

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

# rename gene id to gene name
chromosomes = as.factor(annot[genes,'Chromosome.scaffold.name'])
heatmap_data = tapply(1:length(genes), chromosomes, function(cols) {
  rowMeans(common_data[, cols])
})
heatmap_data = do.call(cbind, heatmap_data)

# option 1 - use all genes
heatmap_data = common_data[,genes]

# option 2 - use PCA
pr = prcomp(common_data[,genes])
nPCs = 15
heatmap_data = pr$x[,1:nPCs]

# samples appear as columns, genes appear as rows
fig = pheatmap(t(heatmap_data),
               annotation_col = common_fhr['risk'],
               annotation_colors = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222")),
               scale = "none",
               color = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100),
               clustering_distance_cols = 'euclidean',
               clustering_distance_rows = 'euclidean',
               clustering_method = 'complete',
               show_rownames=T, 
               show_colnames=F,
               cutree_cols = 1,
               cutree_rows = 1,
               treeheight_row = 0,
               cluster_rows=F,
               main="Copy number of 6460 genes (Chi-sq FWER-adjusted p < 0.05)",
)

save_pheatmap(fig, './assets/pheatmap-cn-by-gene-6460genes-pc.png',width=8,height=8)


barplot(table(chromosomes),main='6460 genes with significant association to risk',xlab='Chromosome')