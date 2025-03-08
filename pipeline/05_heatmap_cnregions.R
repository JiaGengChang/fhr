library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')
source('./utils/common_fhr_data.R')

cn = read.delim('./matrices/segment_cn_matrix_uncorrelated.tsv',sep='\t')


# anova-test to identify significant regions
fhr.0 = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
# Option 1 - compare FHR, GHR, SR
fhr = fhr.0[(fhr.0['risk'] != -1),]

subset_to_labelled(data,fhr)

# test each feature
p_values = rep(1,ncol(common_data))
for(i in 1:ncol(common_data)){
  # ANOVA test for SR-GHR-FHR
  aov_result = aov(common_data[,i] ~ common_fhr$risk)
  p_values[i] <- summary(aov_result)[[1]]["Pr(>F)"][[1]][1]
  # t test for SR-FHR
  # p_values[i] = t.test(common_data[common_fhr$risk=='SR',i],
                       # common_data[common_fhr$risk=='FHR',i])$p.value
}
p_values = apply(common_data, 2, function(...)chisq.test(table(common_fhr$risk,...))$p.value)

p_values_BH = p.adjust(p_values, method='BH')
length(p_values_BH[p_values_BH < 0.05])
p_values_bf = p.adjust(p_values, method='bonferroni')
length(p_values_bf[p_values_bf < 0.05])

# pick genes to use for heatmap
features = colnames(common_data)[p_values_bf < 0.05]
# option 2 - use top 100
features = colnames(common_data)[order(p_values_bf,decreasing=F)[1:100]]

# option 1 - use raw features
heatmap_data = common_data[,features]
# option 2 - use PCs
pr = prcomp(common_data[,features],scale=T)
nPCs = 10
heatmap_data = pr$x[,1:nPCs]

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))

# pick colors for heatmap
heatmap_colors = colorRampPalette(rev(brewer.pal(n=5, name="RdBu")))(5)

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
               cluster_rows=F,
               main="397 CN regions (Chi-sq FWER<0.05)",
)

save_pheatmap(fig, './assets/pheatmap-CNregions-integerCN-397regions-PCs.png',width=12,height=8)

# PCA version
save_pheatmap(fig, './assets/pheatmap-CNregions-PC10.png',width=12,height=8)

# plot loadings of top PCs (PC 1 and 2)
loadings = pr$rotation[,1:nPCs]
pheatmap(loadings,cluster_rows = F,cluster_cols=F,show_rownames = F)


# add chromosome names to heatmap
chroms = unlist(lapply(strsplit(features,'\\.'),function(x)x[[1]][[1]]))
chroms.dedup = rep("", length(chroms))
chroms.dedup[1] = chroms[1]
for (i in 2:length(chroms)){
  if (chroms[i] != chroms[i-1]){
    chroms.dedup[i] = chroms[i]
  }
}
# assign plotting names
colnames(heatmap_data) = chroms.dedup
# now set show_rownames=T and re-plot