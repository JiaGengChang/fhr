library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')

scores1 = read.delim('./scores/DGE_FHR_v_SR.tsv',row.names=1)
genes1 = rownames(scores1[scores1$adj.P.Val < 0.05,])
scores2 = read.delim('./scores/DGE_FHR_v_GHR.tsv',row.names=1)
genes2 = rownames(scores2[1:50,])
scores3 = read.delim('./scores/DGE_GHR_v_SR.tsv',row.names=1)
genes3 = rownames(scores3[1:50,])

genes = union(genes1, union(genes2, genes3))

# read gene expression data and risk labels
data = read.delim('./matrices/gene_exp_logcpm.tsv',row.names = 1)

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

# subset to common ids
common_ids = intersect(rownames(data), rownames(fhr))
common_data = data[rownames(data) %in% common_ids, dge_genes_fhr]
common_fhr = fhr[rownames(fhr) %in% common_ids,]

# match the order of observations
common_data = common_data[match(common_ids, rownames(common_data)),]
common_fhr = common_fhr[match(common_ids, rownames(common_fhr)),]
# rename risk labels from integers to names
common_fhr$risk = factor(common_fhr$risk,levels=c(0,1,2),labels=c("SR","GHR","FHR"))

# rename gene id to gene name
colnames(common_data) = annot[colnames(common_data),'Gene.name']

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))
# annotation_colours = list(risk = c(GHR="#F1A340",SR="#F7F7F7",FHR="#998EC3"))

# pick colors for heatmap
heatmap_colors = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

# samples appear as columns, genes appear as rows
fig = pheatmap(t(common_data),
         annotation_col = common_fhr['risk'],
         annotation_colors = annotation_colours,
         scale = "row",
         color = heatmap_colors,
         clustering_distance_cols = 'euclidean',
         clustering_distance_rows = 'euclidean',
         clustering_method = 'complete',
         show_rownames=T, 
         show_colnames=F,
         cutree_cols = 3,
         cutree_rows = 3,
         treeheight_row = 0,
         cluster_rows=T,
         cex=1,
         main="Expression of top 100 genes for SR vs GHR/FHR",
         )

save_pheatmap(fig, './assets/pheatmap-sr-100-genes.png',height=18)

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
