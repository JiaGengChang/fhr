library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')

scores = read.delim('./scores/gene_mut_scores.tsv',row.names=1)

genes = rownames(scores[(scores['chi2_p']<=0.05) & (scores['freq'] >= 50),])

data = read.delim('./matrices/gene_mut_matrix.tsv.gz',row.names=1)
# drop sample column
data = data[,2:ncol(data)]

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

# identify the names of labeled patients
common_ids = intersect(rownames(data), rownames(fhr))

# subset to common ids
common_data = data[rownames(data) %in% common_ids, genes]
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

# normalize by the mutation burden of each patient
# or by mutation burden over the genes of interest of each patient
common_data_norm = common_data / rowSums(common_data + 0.5)

# samples appear as columns, genes appear as rows
fig = pheatmap(t(common_data_norm),
               annotation_col = common_fhr['risk'],
               annotation_colors = annotation_colours,
               scale = "row",
               color = heatmap_colors,
               clustering_distance_cols = 'maximum',
               clustering_distance_rows = 'maximum',
               clustering_method = 'complete',
               show_rownames=T, 
               show_colnames=F,
               cutree_cols = 1,
               cutree_rows = 1,
               treeheight_row = 0,
               treeheight_col = 30,
               cluster_rows=T,
               cex=1,
               main="NS Mutation of top 44 significant genes (chi-squared)",
)

save_pheatmap(fig, './assets/pheatmap-mut-top44.png',height=8)


# Prepare data for chi2 test
X = data[rownames(data) %in% common_ids,]
X = X[match(common_ids, rownames(X)),]

y = fhr[rownames(fhr) %in% common_ids,]
y = y[match(common_ids, rownames(y)),'risk']#==2

# perform chi2 test for the 3 risk categories
p_values = apply(X, 2, function(...){chisq.test(table(...,y))$p.value})

freq = colSums(X)

# number of significant genes before multiple testing correction
length(p_values[p_values < 0.05 ])

# BH adjustment to control the FDR
p_values_BH <- p.adjust(p_values, method = "BH")

# Bonferroni adjustment to control the FWER
p_values_Bonferroni <- p.adjust(p_values, method = "bonferroni")

length(p_values_BH[p_values_BH < 0.05]) # 80
length(p_values_Bonferroni[p_values_Bonferroni < 0.05]) # 52

# pick the genes to use for the heatmap
genes = colnames(X)[(p_values_BH < 0.05) & colSums(X) >= 2]
length(genes)

pvals.df = cbind(annot[colnames(X),1:2],p_values,p_values_BH,p_values_Bonferroni)
rownames(pvals.df) = colnames(X)
write.table(pvals.df,'./scores/gene-mut-p-vals.tsv',quote=F)
