library(iClusterPlus)
library(lattice)
library(gplots)
source('./utils/datasets.R')

load_common_data()

nK=10

bayfit = tune.iClusterBayes(cpus=nK, 
                            dt1=common.data.rna,
                            dt2=common.data.mut,
                            dt3=common.data.cn,
                            dt4=common.data.ig,
                            type=c("gaussian","binomial","poisson","binomial"),
                            K = 1:nK,)

saveRDS(bayfit, file="./scores/iClusterBayesK1to10.multivariate.rds")

allBIC = sapply(1:nK, function(i) bayfit$fit[[i]]$BIC)
devratio = sapply(1:nK, function(i) bayfit$fit[[i]]$dev.ratio)

par(mar=c(4.0,4.0,0.5,0.5),mfrow=c(1,2))
mypch = rep(1,nK); mypch[5] = 19
plot(1:nK, allBIC,type="b",xlab="k",ylab="BIC",pch=mypch)
plot(1:nK, devratio,type="b",xlab="k",ylab="Deviance ratio",pch=mypch)

# color scheme for data types
bw.col = colorpanel(2,low="white",high="black")
col.scheme = alist()
col.scheme[[1]] = bluered(256) # expr
col.scheme[[2]] = bw.col # mut
col.scheme[[3]] = bluered(5) # cn
col.scheme[[4]] = bw.col # IG

# plot heatmap for optimal K
k=5
plotHeatmap(bayfit$fit[[k]], list(common.data.rna, common.data.mut,common.data.cn,common.data.ig),
            type=c("gaussian","binomial","poisson","binomial"),
            threshold=c(1,1,1,1),
            col.scheme = col.scheme,
            chr=as.integer(sub('^X([0-9]+)[pq].*','\\1',colnames(common.data.cn))),
            width=6,
            row.order=c(T,T,F,T),
            plot.chr = c(F,F,T,F),
            sparse = c(F,F,F,F),
            cap=c(F,F,T,F))


# plot posterior probabilities
par(mar=c(2,4,1,1), mfrow=c(4,1))
plot(bayfit$fit[[k]]$beta.pp[[1]],xlab="Genes",ylab="Posterior proba",
     main="Expression")
plot(bayfit$fit[[k]]$beta.pp[[2]],xlab="Genes",ylab="Posterior proba",
     main="Mutation")
plot(bayfit$fit[[k]]$beta.pp[[3]],xlab="Genomic region",ylab="Posterior proba",
     main="Copy Number")
plot(bayfit$fit[[k]]$beta.pp[[4]],xlab="IG",ylab="Posterior proba",
     main="Canonical IG")


table(risk=common.data.label, cluster=bayfit$fit[[k]]$clusters)

#' Need to visualize this with pheatmap
#' given the labels and data, arrange column and rows in same way
