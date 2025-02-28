library(GenomicRanges)
library(iClusterPlus)
# heatmap
library(lattice)
library(gplots)
source('./utils/copynumber.R')

# unsupervised clustering with iClusterBayes

clin = read.delim('./annotations/fhr-annotations-raw.tsv',row.names = 1)
clin = clin[clin$risk != -1,]
clin$risk = factor(clin$risk, levels=c(0,1,2), labels=c("SR","GHR",'FHR'))

# option 1 - cn segments, segment means
cn = read.delim('./matrices/gene_cn_matrix_CNregions.tsv',sep='\t')
names_cn = substr(rownames(cn),1,9)
chroms_cn = sapply(colnames(cn), function(...) as.integer(strsplit(strsplit(...,"\\.")[[1]][[1]],"chr")[[1]][[2]]))
cn = cn[!duplicated(names_cn), colnames(cn)[order(chroms_cn)]]
rownames(cn) = substr(rownames(cn),1,9)
cn[cn < -9] = min(cn[cn > -9]) # clip the -30 values to -8.79

# option 2 - broad CN features, integer level
cn = read.delim('./matrices/broad_cn_matrix.tsv',sep='\t',row.names=1)
colnames(cn) = sub(".*Cp_([0-9]+[pq][0-9]+).*", "\\1", colnames(cn))
sortednames = c("1p22","1q21","2p23","2q22","3p22","3q21","4p15","4q31","5p15","5q31","6p22","6q25","7p14","7q22","8p22","8q24","9p13","9q33","10p14","10q23","11p15","11q23","12p13","12q21","13q14","13q34","14q23","15q15","15q26","16p12","16q22","17p13","17q23","18p11","18q21","19p13","19q13","20p12","20q13","21q21","21q22","22q13")
cn = get_cn2_df(cn)[,sortednames]
chroms_cn = as.integer(sub("([0-9]+).*", "\\1", colnames(cn)))

expr = read.delim('./matrices/gene_exp_logcpm_proteincoding.tsv',row.names=1)
# pick the K most variable genes (K ~ N)
vars = apply(expr, 2, var)
exprv = expr[,order(vars,decreasing=T)[1:1000]]

ensg = read.delim('./annotations/ensgid-autosomal-proteincoding.txt',row.names = 1)
gene.names = ensg[colnames(exprv),'Gene.name']
colnames(exprv) = ifelse(gene.names=="",colnames(exprv),gene.names)

mut = read.delim('./matrices/gene_mut_matrix_gt1.tsv.gz',row.names=1)
mut = mut[,2:ncol(mut)]
mut[mut > 1] = 1 # convert to binomial
gene.names.mut = ensg[colnames(mut),'Gene.name']
colnames(mut) = ifelse(gene.names.mut=="",colnames(mut),gene.names.mut)
freq = apply(mut, 2, mean)
mutf = mut[, which(freq>=0.02)]

common.ids = Reduce(intersect, list(rownames(exprv),rownames(mutf),rownames(cn),rownames(clin)))

data.mut = as.matrix(mutf[common.ids,])
data.cn = as.matrix(cn[common.ids,])
data.expr = as.matrix(exprv[common.ids,])
data.y = clin[common.ids,'risk']

nK=10
bayfit = tune.iClusterBayes(cpus=nK, dt1=data.expr, dt2=data.cn, dt3=data.mut, 
                             type=c("gaussian","poisson","binomial"),
                             K = 1:nK,)

saveRDS(bayfit, file="./scores/iClusterBayesK1to10.rds")

allBIC = sapply(1:nK, function(i) bayfit$fit[[i]]$BIC)
devratio = sapply(1:nK, function(i) bayfit$fit[[i]]$dev.ratio)

# determine the number of clusters
png('./assets/iClusterBayes.tuning.png',width=1280,height=480,res=200)
par(mar=c(4.0,4.0,0.5,0.5),mfrow=c(1,2))
mypch = rep(1,nK); mypch[5] = 19
plot(1:nK, allBIC,type="b",xlab="k",ylab="BIC",pch=mypch)
plot(1:nK, devratio,type="b",xlab="k",ylab="Deviance ratio",pch=mypch)
dev.off()

# color scheme for data types
bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bluered(256) # expr
col.scheme[[2]] = bluered(5) # cn
col.scheme[[3]] = bw.col # mut

# plot heatmap for optimal K
k=5
plotHeatmap(bayfit$fit[[k]], list(data.expr, data.cn, data.mut),
            type=c("gaussian","poisson","binomial"),
            threshold=c(0.1,0.25,0.1),
            col.scheme = col.scheme,
            chr=chroms_cn,
            width=10,
            row.order=c(T,F,T),
            plot.chr = c(F,T,F),
            sparse = c(T,F,T),
            cap=c(F,T,F))


# plot posterior probabilities
par(mar=c(2,4,1,1), mfrow=c(3,1))
plot(bayfit$fit[[k]]$beta.pp[[1]],xlab="Genes",ylab="Posterior proba",
     main="Expression")
plot(bayfit$fit[[k]]$beta.pp[[2]],xlab="Genomic region",ylab="Posterior proba",
     main="Copy Number")
plot(bayfit$fit[[k]]$beta.pp[[3]],xlab="Genes",ylab="Posterior proba",
       main="Mutation")

table(risk=data.y, cluster=bayfit$fit[[k]]$clusters)
