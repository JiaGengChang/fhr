library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')

data = read.delim('./matrices/gene_cn_matrix.tsv.gz',row.names = 1)
data = data[, 2:ncol(data)] # drop sample column

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

# identify the names of labeled patients
common_ids = intersect(rownames(data), rownames(fhr))

# subset to common ids
common_data = data[rownames(data) %in% common_ids,]
common_fhr = fhr[rownames(fhr) %in% common_ids,]

# match the order of observations
common_data = common_data[match(common_ids, rownames(common_data)),]
common_fhr = common_fhr[match(common_ids, rownames(common_fhr)),]
# rename risk labels from integers to names
# common_fhr$risk = factor(common_fhr$risk,levels=c(0,1,2),labels=c("SR","GHR","FHR"))

# Perform ANOVA
p_values = numeric(ncol(common_data))
for(i in 1:ncol(common_data)){
  aov_result = aov(common_data[,i] ~ common_fhr$risk)
  p_values[i] <- summary(aov_result)[[1]]["Pr(>F)"][[1]][1]
}

# adjust for FDR
p_values_BH = p.adjust(p_values, method='BH')

length(p_values_BH[p_values_BH < 0.05])

# adjust for FWER
p_values_bf = p.adjust(p_values, method='bonferroni')

length(p_values_bf[p_values_bf < 0.05])

# pick genes to use for heatmap
genes = colnames(common_data)[p_values_bf < 0.05]

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

# rename gene id to gene name
chromosomes = as.factor(annot[genes,'Chromosome.scaffold.name'])
heatmap_data = tapply(1:length(genes), chromosomes, function(cols) {
  rowMeans(common_data[, cols])
})
heatmap_data = do.call(cbind, heatmap_data)

heatmap_data = common_data[,genes]

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))

# pick colors for heatmap
heatmap_colors = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

# samples appear as columns, genes appear as rows
fig = pheatmap(t(heatmap_data),
               annotation_col = common_fhr['risk'],
               annotation_colors = annotation_colours,
               scale = "none",
               color = heatmap_colors,
               clustering_distance_cols = 'euclidean',
               clustering_distance_rows = 'euclidean',
               clustering_method = 'complete',
               show_rownames=T, 
               show_colnames=F,
               cutree_cols = 3,
               cutree_rows = 2,
               treeheight_row = 0,
               cluster_rows=T,
               cex=0.9,
               main="Copy number of significant genes (ANOVA)",
)

save_pheatmap(fig, './assets/pheatmap-fhr-vs-rest-14genes.png')