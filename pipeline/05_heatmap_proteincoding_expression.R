library(pheatmap)
library(RColorBrewer)
source('./utils/common_fhr_data.R')
source('./utils/save_pheatmap.R')

# option 1 - use separate DGE scores
scores1 = read.delim('./scores/DGE_FHR_v_SR.tsv',row.names=1)
genes1 = rownames(scores1[scores1$adj.P.Val < 0.05,])
scores2 = read.delim('./scores/DGE_FHR_v_GHR.tsv',row.names=1)
genes2 = rownames(scores2[1:50,])
scores3 = read.delim('./scores/DGE_GHR_v_SR.tsv',row.names=1)
genes3 = rownames(scores3[1:50,])
genes = union(genes1, union(genes2, genes3))

# option 2 - use combined DGE scores


# read gene expression data and risk labels
data = read.delim('./matrices/gene_exp_logcpm_proteincoding.tsv')

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names=1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

# get common ids, data and fhr labels
subset_to_labelled(data, fhr)

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))

# pick colors for heatmap
heatmap_colors = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

# pick genes for heatmap
heatmap_data = common_data[,genes.signif]
# rename gene id to gene name
colnames(heatmap_data) = annot[colnames(heatmap_data),'Gene.name']

# samples appear as columns, genes appear as rows
fig = pheatmap(t(heatmap_data),
         annotation_col = common_fhr['risk'],
         annotation_colors = annotation_colours,
         scale = "row",
         color = heatmap_colors,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'complete',
         show_rownames=T, 
         show_colnames=F,
         cutree_cols = 1,
         cutree_rows = 1,
         treeheight_row = 0,
         cluster_rows=T,
         cex=1,
         main="Top 274 protein coding DEGs, logFC>1, p<0.05",
         )

save_pheatmap(fig, './assets/pheatmap-rna-274-genes.png',height=14,width=10)

# DEPRECATED - using linear model is a better approach
# perform t-test here
X = data[rownames(data) %in% common_ids, ]
X = X[match(common_ids, rownames(X)),]
p_values = rep(1, ncol(X))
# possible options are FHR vs rest, FHR vs GHR
for (i in 1:ncol(X)){
  x = X[common_fhr$risk != 'SR', i] # group 1
  y = X[common_fhr$risk == 'SR', i] # group 2
  # ignore columns with too low variance
  if (var(X[,i]) > 1){
    p_values[i] = t.test(x, y)$p.value
  }
}
p_values_BH = p.adjust(p_values,method='BH')
p_values_bf = p.adjust(p_values,method='bonferroni')

length(colnames(X)[p_values_BH < 0.05])
length(colnames(X)[p_values_bf < 0.05])

dge_genes_fhr = colnames(X)[p_values_bf < 0.05]
dge_genes_fhr = colnames(X)[order(p_values_bf)[1:100]]

pval.df = cbind(annot[colnames(X),],p_values,p_values_BH,p_values_bf)
rownames(pval.df) = colnames(X)
write.table(pval.df,'./scores/gene-exp-sr-p-vals.tsv')
