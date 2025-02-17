library(edgeR)
library(limma)

counts = read.delim('./data/Expression Estimates - Gene Based_MMRF_CoMMpass_IA22_salmon_geneUnstranded_counts.tsv.gz',row.names = 1)
public_ids = substr(colnames(counts),1,9)
colnames(counts) = public_ids
counts = counts[,!duplicated(colnames(counts))]

genedata = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)

# subset to the protein coding genes
counts = counts[rownames(counts) %in% rownames(genedata),]

clindata = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
clindata = clindata[clindata$risk != -1,]
clindata$risk = factor(clindata$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))

common_id = intersect(colnames(counts), rownames(clindata))

counts_common = counts[,common_id]
clin_common = clindata[common_id,]

# convert to DGE object, rows are genes.
d0 = DGEList(counts_common, genes=genedata[rownames(counts_common),'Gene.name'])
d0 = calcNormFactors(d0)

# Filtering genes

# group by risk
mm = model.matrix(~0 + risk, data=clin_common)
keep = filterByExpr(d0, mm)
sum(keep) # number of genes retained

# filter out lowly expressed genes
d = d0[keep,]

# voom transformation and calculation of variance weights
y = voom(d, mm, plot=T)
# y = voom(d0, mm, plot=T) # before removing low expr genes

# fit a linear model using weighted least squares for each gene
fit = lmFit(y, mm)
head(coef(fit))

# comparison between groups 
# introduce (FHR - SR) * (FHR - GHR)
contr = makeContrasts(riskFHR - riskSR - riskGHR, levels = colnames(coef(fit)))

# estimate contrast for each gene
tmp = contrasts.fit(fit, contr)

# apply Empirical Bayes smoothing of standard errors of log FC
tmp.smooth = eBayes(tmp)

# multiple testing adjustment
# n=Inf says to produce the topTable for all genes
top.table <- topTable(tmp.smooth, adjust.method = "BH", sort.by = "P", n = Inf)
head(top.table)

# how many DE genes are there (false discovery rate corrected)?
length(which(top.table$adj.P.Val < 0.05)) # 3873 for FHR-GHR, 32 for FHR-SR, 7807 for GHR-SR

write.table(top.table, file="./scores/DGE_FHR_v_GHR.tsv",sep='\t')
write.table(top.table, file="./scores/DGE_FHR_v_SR.tsv",sep='\t')
write.table(top.table, file="./scores/DGE_GHR_v_SR.tsv",sep='\t')

write.table(top.table, file="./scores/DGE_FHR_v_rest.tsv",sep='\t')
