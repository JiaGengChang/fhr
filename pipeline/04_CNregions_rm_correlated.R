source('./utils/save_pheatmap.R')
source('./utils/common_fhr_data.R')

cn = read.delim('./matrices/gene_cn_matrix_CNregions.tsv',sep='\t')
public.id = substr(rownames(cn),1,9)
chroms_cn = as.integer(sub('chr([0-9]+)\\.[0-9]+\\.[0-9]+','\\1',colnames(cn)))
cn = cn[!duplicated(public.id), colnames(cn)[order(chroms_cn)]]
rownames(cn) = public.id[!duplicated(public.id)]

# locations of values < -9
missing.ind = arrayInd(which(cn < -9),dim(cn))
missing.row = unique(missing.ind[,1])

# Option 1: drop individuals with missing CN data
# cn = cn[-missing.row,]
# Option 2: clip the -30 values to -8.79
# cn[cn < -9] = min(cn[cn > -9])

# Optional - Convert segment mean to integer CN
# data = cn
source('./utils/copynumber.R')
data = apply(apply(cn, 2, seg_mean_to_cn), 2, as.integer)
rownames(data) = rownames(cn)


# Chromosome map 

# check for correlated features
cormat = cor(data)
chroms = sub('chr([0-9]+).*','\\1',colnames(cormat))
chroms[duplicated(chroms)] = ""

# trick to rotate tick labels by 90
{
  library(grid)
  draw_colnames_90 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 1, hjust = 0.5, rot = 90, gp = gpar(...))
    return(res)}
  
  ## 'Overwrite' default draw_colnames with your own version 
  assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                    ns=asNamespace("pheatmap"))
}
pheatmap(cormat,
         cluster_cols = F, 
         cluster_rows=F, 
         treeheight_row = 0, 
         treeheight_col = 0, 
         show_rownames = T, 
         show_colnames = T,
         labels_row = chroms,
         labels_col=chroms,
         scale = 'none',
         main = "Correlation in genomic copy number segments across 924 samples")


# filter correlated features
{
  cormask = cormat
  cormask[upper.tri(cormask)]=0
  diag(cormask) = 0
}
# number of collinear and non-collinear features
is.collinear = apply(cormask, 2, function(x) sum(abs(x) >= 0.9) > 1); table(is.collinear)

# number of non-collinear features by chromosome
table(chroms_cn[!is.collinear])

# location of non-collinear features
for (idx in which(!is.collinear)){
  cormask[max(0,idx-30):idx, idx] = 1
}
pheatmap(cormask,cluster_cols=F,cluster_rows=F,show_rownames=T,show_colnames=T,labels_row=chroms,labels_col=chroms)

subset = cormat[!is.collinear,!is.collinear]
labels = chroms_cn[!is.collinear]; labels[duplicated(labels)]=""
pheatmap(subset,cluster_cols=F,cluster_rows=F,show_rownames=T,show_colnames=T,labels_row=labels,labels_col=labels)

cn.out = cn[!is.collinear]
write.table(cn.out,'./matrices/gene_cn_matrix_CNregions_uncorrelated.tsv',sep='\t',quote=F)
