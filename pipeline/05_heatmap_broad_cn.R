library(pheatmap)
library(RColorBrewer)
source('./utils/save_pheatmap.R')


data = read.delim('./matrices/broad_cn_matrix.tsv',row.names=1)
colnames(data) = substr(colnames(data),11,15)
sortednames = c("1p22","1q21","2p23","2q22","3p22","3q21","4p15","4q31","5p15","5q31","6p22","6q25","7p14","7q22","8p22","8q24","9p13","9q33","10p14","10q23","11p15","11q23","12p13","12q21","13q14","13q34","14q23","15q15","15q26","16p12","16q22","17p13","17q23","18p11","18q21","19p13","19q13","20p12","20q13","21q21","21q22","22q13")
data = data[,sortednames]

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
common_fhr$risk = factor(common_fhr$risk,levels=c(0,1,2),labels=c("SR","GHR","FHR"))

# pick colors for risk labels
annotation_colours = list(risk = c(SR="#F9F9F9",GHR="#22FF22",FHR="#FF2222"))
# annotation_colours = list(risk = c(GHR="#F1A340",SR="#F7F7F7",FHR="#998EC3"))

# pick colors for heatmap
# heatmap_colors = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(100)
heatmap_colors = rev(brewer.pal(n=11, name="RdBu"))

# samples appear as columns, genes appear as rows
fig = pheatmap(t(common_data),
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
               main="WGS-based chromosome-arm level copy number",
)

save_pheatmap(fig, './assets/pheatmap-broad-cn.png')

