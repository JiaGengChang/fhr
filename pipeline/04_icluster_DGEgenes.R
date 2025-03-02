library(GenomicRanges)
library(iClusterPlus)
# heatmap
library(lattice)
library(gplots)

# unsupervised clustering with iClusterBayes
# features are selected via supervision
# using 77 DGE genes instead of 1000 variable genes 
# using 42 broad CN features instead of 1464 CN regions
# mutations is the same, 18 genes binary status

clin = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
clin = clin[clin$risk != -1,]
clin$risk = factor(clin$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))

cn = read.delim('./matrices/broad_cn_matrix.tsv',row.names=1)
colnames(cn) = substr(colnames(cn),11,15)
sortednames = c("1p22","1q21","2p23","2q22","3p22","3q21","4p15","4q31","5p15","5q31","6p22","6q25","7p14","7q22","8p22","8q24","9p13","9q33","10p14","10q23","11p15","11q23","12p13","12q21","13q14","13q34","14q23","15q15","15q26","16p12","16q22","17p13","17q23","18p11","18q21","19p13","19q13","20p12","20q13","21q21","21q22","22q13")
cn = cn[,sortednames]

expr = read.delim('./matrices/gene_exp_logcpm.tsv',row.names=1)
scores1 = read.delim('./scores/DGE_FHR_v_SR.tsv',row.names=1)
genes1 = rownames(scores1[scores1$adj.P.Val < 0.05,])
scores2 = read.delim('./scores/DGE_FHR_v_GHR.tsv',row.names=1)
genes2 = rownames(scores2[1:30,])
scores3 = read.delim('./scores/DGE_GHR_v_SR.tsv',row.names=1)
genes3 = rownames(scores3[1:30,])
genes = union(genes1, union(genes2, genes3))
exprv = expr[, genes]

ensg = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)
gene.names = ensg[colnames(exprv),'Gene.name']
colnames(exprv) = ifelse(gene.names=="",colnames(exprv),gene.names)

mut = read.delim('./matrices/gene_mut_matrix_gt1.tsv.gz',row.names=1)
mut = mut[,2:ncol(mut)]
mut[mut > 1] = 1 # convert to binomial
gene.names.mut = ensg[colnames(mut),'Gene.name']
colnames(mut) = ifelse(gene.names.mut=="",colnames(mut),gene.names.mut)
freq = apply(mut, 2, mean)
mutf = mut[, which(freq>0.05)]

common.ids = Reduce(intersect, list(rownames(exprv),rownames(mutf),rownames(cn),rownames(clin)))

data.mut = as.matrix(mutf[common.ids,])
data.cn = as.matrix(cn[common.ids,])
data.expr = as.matrix(exprv[common.ids,])
data.y = clin[common.ids,'risk']

bayfit = tune.iClusterBayes(cpus=10, data.expr, 
                             type="gaussian",
                             K = 1:10,)

saveRDS(bayfit, file="./scores/iClusterBayesK1to10_DGEgenes.rds")

nK=10
allBIC = sapply(1:nK, function(i) bayfit$fit[[i]]$BIC)
devratio = sapply(1:nK, function(i) bayfit$fit[[i]]$dev.ratio)

# determine the number of clusters
png("./assets/iClusterBayes.tuning.supervised.png",width=2000,height=800,res=200)
par(mar=c(4,4,0.5,0.5),mfrow=c(1,2))
mypch = rep(1,nK); mypch[2] = 19
plot(1:nK, allBIC,type="b",xlab="k",ylab="BIC",pch=mypch)
plot(1:nK, devratio,type="b",xlab="k",ylab="Deviance ratio",pch=mypch)
dev.off()

# color scheme for data types
bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bluered(256) # expr
col.scheme[[2]] = bluered(256) # cn
col.scheme[[3]] = bw.col # mut

# plot heatmap for optimal K
k=2
chr=as.integer(sapply(colnames(cn), function(...) strsplit(strsplit(..., "p")[[1]][[1]],"q")[[1]][[1]]))
plotHeatmap(bayfit$fit[[k]], list(data.expr),
            type=c("gaussian"),
            threshold=c(0.25),
            col.scheme = col.scheme,
            chr=chr,
            width=6,
            row.order=c(F),
            plot.chr = c(F),
            sparse = c(F),
            cap=c(F))

# plot posterior probabilities
par(mar=c(2,4,1,1), mfrow=c(3,1))
plot(bayfit$fit[[k]]$beta.pp[[1]],xlab="Genes",ylab="Posterior proba",
     main="Expression")
plot(bayfit$fit[[k]]$beta.pp[[2]],xlab="Genomic region",ylab="Posterior proba",
     main="Copy Number")
plot(bayfit$fit[[k]]$beta.pp[[3]],xlab="Genes",ylab="Posterior proba",
       main="Mutation")

# tabulate number of FHR/GHR/SR per cluster
table(data.y, bayfit$fit[[k]]$clusters)
# data.y   1   2   3
# SR   92 182 109
# GHR  50  59 101
# FHR   7  35  38
