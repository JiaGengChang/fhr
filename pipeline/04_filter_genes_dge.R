library(edgeR)

counts = read.delim('./data/Expression Estimates - Gene Based_MMRF_CoMMpass_IA22_salmon_geneUnstranded_counts.tsv.gz',row.names = 1)
public_ids = substr(colnames(counts),1,9)
colnames(counts) = public_ids
counts = counts[,!duplicated(colnames(counts))]

annot = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)

# subset to the protein coding genes
counts = counts[rownames(counts) %in% rownames(annot),]

clindata = read.delim('./annotations/fhr-annotations.tsv',row.names = 1)
clindata = clindata[clindata$risk != -1,]
clindata$risk = factor(clindata$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))

common_id = intersect(colnames(counts), rownames(clindata))

counts_common = counts[,common_id]
clin_common = clindata[common_id,]

# convert to DGE object, rows are genes.
d0 = DGEList(counts_common, genes=annot[rownames(counts_common),'Gene.name'])
d0 = calcNormFactors(d0)

# Filtering genes

# use risk-related features as variables
# mm.full = model.matrix(~0 + ., data=clin_common)
# use all except risk as variables
# mm.semi = model.matrix(risk~0 + ., data=clin_common)
# use risk itself as variables
mm.risk = model.matrix(~0 + risk, data=clin_common)

keep = filterByExpr(d0,)
sum(keep) # number of genes retained: 15159

# filter out lowly expressed genes
d = d0[keep,]

# voom transformation and calculation of variance weights
y = voom(d, plot=T)

# fit a linear model using weighted least squares for each gene
fit = lmFit(y, mm.risk)
head(coef(fit))

# comparison between groups 
contr = makeContrasts(FHR=riskFHR-(riskGHR+riskSR)/2,
                      GHR=riskGHR - (riskFHR+riskSR)/2,
                      SR=riskSR - (riskGHR+riskFHR)/2,
                      levels = colnames(coef(fit)))

# estimate contrast for each gene
fit2 = contrasts.fit(fit, contr)

# apply Empirical Bayes smoothing of standard errors of log FC
fit2 = eBayes(fit2)

# multiple testing adjustment
# n=Inf says to produce the topTable for all genes
top.table <- topTable(fit2, adjust.method = "BH", sort.by = "F", n = Inf)
head(top.table)

# how many DE genes are there (false discovery rate corrected)?
length(which(top.table$adj.P.Val < 0.05)) # 3873 for FHR-GHR, 32 for FHR-SR, 7807 for GHR-SR

results = decideTests(fit2, method = 'nestedF',lfc=1)
vennDiagram(results)
vennCounts(results)

# 67, or 20, or 274
genes.signif = rownames(results@.Data)[results@.Data[,1]!=0 | results@.Data[,2]!=0 | results@.Data[,3]!=0]

# prepare output dataset
d1 = DGEList(counts, genes=annot[rownames(counts),'Gene.name'])
d1 = calcNormFactors(d1)

# write the gene expression table
cpm.d1 = t(cpm(d1[genes.signif,],log=TRUE)) # 806 by n_genes
write.table(cpm.d1,"./matrices/gene_exp_logcpm_proteincoding_signif274.tsv",sep='\t',quote=F)

# write the test scores
write.table(top.table,file='./scores/DGE_proteincoding_all_contrasts.tsv',sep='\t',quote=F)

write.table(top.table, file="./scores/DGE_FHR_v_GHR.tsv",sep='\t')
write.table(top.table, file="./scores/DGE_FHR_v_SR.tsv",sep='\t')
write.table(top.table, file="./scores/DGE_GHR_v_SR.tsv",sep='\t')
write.table(top.table, file="./scores/DGE_FHR_v_rest.tsv",sep='\t')
