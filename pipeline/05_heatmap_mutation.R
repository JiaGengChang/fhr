library(pheatmap)
library(RColorBrewer)
library(gplots)
source('./utils/save_pheatmap.R')
source('./utils/common_fhr_data.R')
source('./utils/init.R')

# scores from scikit learn feature selection
scores = read.delim('./scores/gene_mut_scores.tsv',row.names=1)
genes = rownames(scores[(scores['chi2_p']<=0.05) & (scores['freq'] >= 10),])

data = read.delim('./matrices/gene_mut_matrix_gt1.tsv.gz',row.names=1)
data = data[,2:ncol(data)]

fhr = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
fhr = fhr[fhr$risk != -1,]

subset_to_labelled(data, fhr)

# Select genes to use for heatmap

# Prepare data for chi2 test
X = common_data
# X = data[rownames(data) %in% common_ids,]
# X = X[match(common_ids, rownames(X)),]

CONVERT_TO_BINARY=F

if (CONVERT_TO_BINARY){
  tmp=rownames(X)
  X = apply(X,2, function(...) as.integer(as.logical(...)))
  rownames(X)=tmp
}

y = common_fhr$risk
# y = fhr[rownames(fhr) %in% common_ids,]
# y = y[match(common_ids, rownames(y)),'risk'] # ordinal vars: 0, 1, 2

# perform chi2 test for the 3 risk categories
p_values = apply(X, 2, function(...){chisq.test(table(...,y))$p.value})

# perform chi2 test for the 3 risk categories
p_values_BH <- p.adjust(p_values, method = "BH")
p_values_Bonferroni <- p.adjust(p_values, method = "bonferroni")

length(p_values_BH[p_values_BH < 0.05]) # 80
length(p_values_Bonferroni[p_values_Bonferroni < 0.05]) # 52

# pick the genes to use for the heatmap
# 44 genes: FDR 0.25
# 33 genes: FDR 0.05
genes = colnames(common_data)[(p_values_BH < 0.05) & colSums(common_data)>=1]
genes = colnames(common_data)[p_values < 0.05]
length(genes) # ENSGID

# write dataset
write.table(data[,genes],'./matrices/gene_mut_matrix_signif33.tsv',sep='\t',quote=F)

# for DAVID
write.table(genes,'genes.txt',quote=F, row.names = F, col.names = F)

heatmap_data = common_data[,genes]
# rename heatmap columns to gene names 
gene.names=annot[colnames(heatmap_data),'Gene.name']
colnames(heatmap_data) = ifelse(gene.names=="",colnames(heatmap_data),gene.names)

# normalize by the mutation burden of each patient
# or by mutation burden over the genes of interest of each patient
heatmap_data_scaled = scale(heatmap_data,center=T,scale=T)

# samples appear as columns, genes appear as rows
fig = pheatmap(t(heatmap_data_scaled),
               annotation_col = common_fhr['risk'],
               annotation_colors = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222")),
               scale = "none",
               color = colorRampPalette(colorpanel(11,low="white",high="black"))(100),
               clustering_distance_cols = 'euclidean',
               clustering_distance_rows = 'euclidean',
               clustering_method = 'complete',
               show_rownames=T, 
               show_colnames=F,
               cutree_cols = 1,
               cutree_rows = 1,
               treeheight_row = 0,
               treeheight_col = 30,
               cluster_rows=T,
               main="Nonsyn. mutations in 36 significant genes (chi-squared FDR<0.05)",
)

save_pheatmap(fig, './assets/pheatmap-mut-top44-greyscale.png',height=8,width=8)

# save scores
pvals.df = cbind(annot[colnames(X),1:2],p_values,p_values_BH,p_values_Bonferroni)
rownames(pvals.df) = colnames(X)
write.table(pvals.df,'./scores/gene-mut-p-vals.tsv',quote=F)

# number of mutations per gene
final.gene.counts = sort(colSums(common_data[,genes]),decreasing = T)
final.gene.names = annot[names(final.gene.counts),'Gene.name']
names(final.gene.counts) = final.gene.names
table(final.gene.counts)


# number of FHR patients with mutation
common_data[common_fhr$risk=='FHR',genes]

#  ENSG00000213281,ENSG00000163517,ENSG00000113811,ENSG00000132382,ENSG00000161547,ENSG00000215695,ENSG00000113638,ENSG00000118496,ENSG00000135439,ENSG00000170275,ENSG00000179627,ENSG00000095970,ENSG00000136881,ENSG00000104415,ENSG00000173801,ENSG00000149485,ENSG00000106565,ENSG00000149474,ENSG00000163914,ENSG00000006210,ENSG00000068078,ENSG00000221954,ENSG00000236981,ENSG00000170153,ENSG00000170035,ENSG00000134317,ENSG00000181718,ENSG00000108061,ENSG00000162441,ENSG00000141510,ENSG00000158092,ENSG00000163961,ENSG00000137312
