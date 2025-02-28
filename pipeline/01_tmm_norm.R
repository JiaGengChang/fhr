library(edgeR)
library(limma)

counts = read.delim('./data/Expression Estimates - Gene Based_MMRF_CoMMpass_IA22_salmon_geneUnstranded_counts.tsv.gz',row.names = 1)
public_ids = substr(colnames(counts),1,9)
colnames(counts) = public_ids
counts = counts[,!duplicated(colnames(counts))]

genedata = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)

# subset to the protein coding genes
counts = counts[rownames(counts) %in% rownames(genedata),]

# convert to DGE object, rows are genes.
genes = genedata[rownames(counts),'Gene.name']
genes = ifelse(genes=="",rownames(counts),genes)
d0 = DGEList(counts, genes=genes)
d0 = calcNormFactors(d0)

# extract log normalized expression 
logcpm = cpm(d0, log=TRUE)
write.table(t(logcpm),'./matrices/gene_exp_logcpm.tsv',sep='\t')

# there is no filtering of lowly expressed genes here. this is done in 04_filter_genes_dge.R