library(edgeR)
library(limma)

counts = read.delim('./data/Expression Estimates - Gene Based_MMRF_CoMMpass_IA22_salmon_geneUnstranded_counts.tsv.gz',row.names = 1)
public_ids = substr(colnames(counts),1,9)
colnames(counts) = public_ids
counts = counts[,!duplicated(colnames(counts))]

annot.mrna = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)
annot.mirna = read.delim('./annotations/ensgid-autosomal-mirna.txt',row.names = 1)
annot.lncrna = read.delim('./annotations/ensgid-autosomal-lncrna.txt',row.names = 1)

# subset to the protein coding genes
clindata = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
clindata = clindata[clindata$risk != -1,]
clindata$risk = factor(clindata$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))

common_id = intersect(colnames(counts), rownames(clindata))
counts_common = counts[,common_id]
clin_common = clindata[common_id,]

d0 = DGEList(counts_common)
d0 = calcNormFactors(d0)

# filter low expression genes
mm = model.matrix(~0 + risk, data=clin_common)
keep = filterByExpr(d0, mm)
sum(keep) # 24694

# y0 = voom(d0,plot=T)
d = d0[keep,]
y = voom(d,plot=T)

logcpm = cpm(d, log=TRUE)

sum(rownames(counts[keep,]) %in% rownames(annot.mirna)) # 9
sum(rownames(counts[keep,]) %in% rownames(annot.lncrna)) # 5187
sum(rownames(counts[keep,]) %in% rownames(annot.mrna)) # 15151

fit = lmFit(y, mm)

contr = makeContrasts(riskFHR-riskSR, riskFHR-riskGHR, levels = colnames(coef(fit)))

fit2 = contrasts.fit(fit, contr)
fit2 = eBayes(fit2)

top.table <- topTable(fit2, adjust.method = "BH", sort.by = "F", n = Inf) # sort by P or F

head(top.table)

# how many DE genes are there (false discovery rate corrected)?
length(which(top.table$adj.P.Val < 0.05)) # 3873 for FHR-GHR, 32 for FHR-SR, 7807 for GHR-SR

#  Large values of eb$t indicate differential expression
results = decideTests(fit2, method="nestedF",lfc=1.2)
vennCounts(results)

genes.fhr = rownames(results@.Data)[(results@.Data[,1] & !results@.Data[,2] & !results@.Data[,3])]
genes.sr = rownames(results@.Data)[(!results@.Data[,1] & !results@.Data[,2] & results@.Data[,3])]
genes.ghr = rownames(results@.Data)[(!results@.Data[,1] & results@.Data[,2] & !results@.Data[,3])]
genes.nonfhr = rownames(results@.Data)[(!results@.Data[,1] & results@.Data[,2] & results@.Data[,3])]

genes = rownames(top.table)[top.table$adj.P.Val < 0.05]
genes = rownames(top.table)[1:100]

logcpm.signif = logcpm[genes,]

dim(logcpm.signif)

write.table(top.table,"./scores/DGE_all_RNA_all_contrasts.tsv")

# prepare logcpm matrix containing all samples, not just those with label
d1 = DGEList(counts)
d1 = calcNormFactors(d1)
y1 = voom(d1[keep,],plot=T)

logcpm1 = cpm(d1[keep,], log=TRUE)
tlogcpm1 = t(logcpm1)
write.table(tlogcpm1,'./matrices/gene_exp_logcpm_totalRNA.tsv')
