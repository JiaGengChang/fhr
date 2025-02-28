library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')
source('./utils/common_fhr_data.R')

# Use combined DGE scores
scores = read.delim('./scores/DGE_totalRNA_all_contrasts.txt',sep=' ',row.names=1)
genes = rownames(scores)[scores$adj.P.Val < 0.05]
genes = rownames(scores)[1:100] # top 100 total RNA
genes = rownames(scores)[rownames(scores) %in% rownames(annot.mrna)][1:100] # top 100 protein coding autosomal

# read gene expression data and risk labels
data = read.delim('./matrices/gene_exp_logcpm_totalRNA.txt',sep=' ',row.names = 1)

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-all.txt',row.names=1)

# subset to common ids
subset_to_labelled(data, fhr)
common_fhr$risk = replace(as.character(common_fhr$risk),is.na(common_fhr$risk),"NR")

# pick colors for risk labels
annotation_colours = list(risk = c(NR="#999999",SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))

# pick colors for heatmap
heatmap_colors = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)

# subset to genes
heatmap_data = common_data[,genes]
gene.names = annot[colnames(heatmap_data),'Gene.name']
colnames(heatmap_data) = ifelse(is.na(gene.names) | gene.names=="",colnames(heatmap_data),gene.names)


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
               main=sprintf("All genes"),
)

save_pheatmap(fig, './assets/pheatmap-mRNA-proteincoding-100-genes.png',height=14,width=20)

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
