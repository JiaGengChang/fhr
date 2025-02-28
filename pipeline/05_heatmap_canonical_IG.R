library(pheatmap)
library(RColorBrewer)
source('./utils/common_fhr_data.R')
source('./utils/save_pheatmap.R')

data = read.delim('./matrices/canonical_ig_translocations.tsv',sep='\t')
fhr = read.delim('./annotations/fhr-annotations.tsv',row.names=1)
fhr = fhr[(fhr['risk'] != -1),]

subset_to_labelled(data,fhr)

annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))
heatmap_colors = c("#EEEEEE","black")
heatmap_data = common_data

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
               treeheight_row = 0,
               cluster_rows=T,
               cex=1,
               main="Canonical IG translocations",
)

save_pheatmap(fig, './assets/pheatmap-canonical-IG.png',height=8,width=8)

p.values = apply(common_data, 2, function(x) fisher.test(table(x,common_fhr$risk))$p.value)
# > p.values
# WHSC1        CCND3          MYC         MAFA        CCND1        CCND2          MAF         MAFB 
# 1.200475e-70 9.032373e-02 9.825739e-02 6.939878e-01 2.447600e-11 1.605104e-01 1.277680e-20 8.992261e-01 
# > which(p.values < 0.05)
# WHSC1 CCND1   MAF 
# 1     5     7 

p.values.adj = p.adjust(p.values, method='BH')
# > p.values.adj
# WHSC1        CCND3          MYC         MAFA        CCND1        CCND2          MAF         MAFB 
# 9.603797e-70 1.572118e-01 1.572118e-01 7.931289e-01 6.526933e-11 2.140139e-01 5.110721e-20 8.992261e-01 
# > which(p.values.adj < 0.05)
# WHSC1 CCND1   MAF 
# 1     5     7 

# write significant features to dataset
res = data[names(which(p.values.adj<0.05))]
write.table(res,'./matrices/canonical_ig_translocations_signif3.tsv',sep='\t',quote=F)
