library(pheatmap)
library(RColorBrewer)
library(gplots)
source('./utils/save_pheatmap.R')
source('./utils/common_fhr_data.R')

data = read.delim('./matrices/gene_fusion_matrix.tsv',row.names = 1)
data = data[, 2:ncol(data)] # drop sample column

fhr = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
fhr = fhr[(fhr['risk'] != -1),] # remove NA labels

subset_to_labelled(data, fhr) # 698
common_data = common_data[,colSums(common_data)>0]

# Chi-square test
p_values = apply(common_data,2,function(...){chisq.test(table(...,common_fhr$risk))$p.value})

# Perform two group rank sum test for FHR vs non-FHR
# p_values = numeric(ncol(common_data))
# for(i in 1:ncol(common_data)){
#   x = as.numeric(risk[common_data[,i] > 0] == 2) # are carriers FHR
#   y = as.numeric(risk[common_data[,i]==0]==2) # are non-carriers FHR
#   p_values[i] = wilcox.test(common_data[,i], as.numeric(common_fhr$risk==2))$p.value
# }

# adjust for FDR
p_values_BH = p.adjust(p_values, method='BH')
length(p_values_BH[p_values_BH < 0.05 & !is.na(p_values_BH)])

# adjust for FWER
p_values_bf = p.adjust(p_values, method='bonferroni')
length(p_values_bf[(p_values_bf < 0.05) & !is.na(p_values_bf)])

# pick genes to use for heatmap
# option 1 - use test results
genes = colnames(common_data)[(p_values_bf < 0.05) & !is.na(p_values_bf) & (colSums(common_data) >= 1)]
genes = colnames(common_data)[(p_values_BH < 0.05) & !is.na(p_values_BH)]

length(genes)

# write dataset 
write.table(data[,genes], './matrices/gene_fusion_matrix_signif6.tsv',sep='\t',quote=F)

freq = colSums(common_data)
pvals = cbind(p_values,p_values_bf,p_values_BH,freq)
rownames(pvals) = colnames(data)
write.table(pvals,'./scores/gene-fusion-p-vals.tsv',quote=F)

# option 2 - go by frequency
genes = names(sort(colSums(common_data),decreasing = T)[1:10])

# read gene id to gene name pairs
annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names=1)

heatmap_data = common_data[,genes]

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))

# pick colors for heatmap
heatmap_colors = brewer.pal(n=3, name="Blues")
heatmap_colors = colorpanel(2,low="#F7F7F7",high="black")

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
               cutree_cols = 1,
               cutree_rows = 1,
               treeheight_row = 0,treeheight_col = 30,
               cluster_rows=T,
               cex=1,
               main="Top gene fusion events by Chi-Squared test",
)

save_pheatmap(fig, './assets/pheatmap-gene-fusion-chisq6.png',height=3,width=8)
