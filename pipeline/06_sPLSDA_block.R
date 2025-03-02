library(mixOmics)
source('./utils/datasets.R')

load_common_data()

X = list(rna=common.data.rna, cna=common.data.cn, mut=common.data.mut, ig=common.data.ig)
Y = common.data.label$risk

design = matrix(1, ncol=length(X), nrow=length(X), dimnames=list(names(X),names(X)))
diag(design) = 0

list.keepX = list(rna=c(772),
                  cna=c(1,3,42),
                  mut=c(1,31),
                  ig=c(1,3))

BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)

varnames = c(names(list.keepX),'Y')
design = matrix(c(0,1,0,1,1,
                  1,0,0,0,1,
                  0,0,0,0,1,
                  1,0,0,0,1,
                  1,1,1,1,0),nrow = 5, ncol=5, dimnames = list(varnames,varnames))

tune.block.splsda.results = tune.block.splsda(X, Y,
                                              ncomp=c(2),
                                              test.keepX=list.keepX,
                                              design = design,
                                              validation="Mfold",
                                              folds=5,
                                              nrepeat=10,
                                              BPPARAM=BPPARAM)

plot(tune.block.splsda.results, col = color.jet(2))

tune.block.splsda.results$choice.ncomp$ncomp
tune.block.splsda.results$choice.keepX

optimal.ncomp = tune.block.splsda.results$choice.ncomp$ncomp
optimal.keepX = lapply(tune.block.splsda.results$choice.keepX,function(...)...[1:optimal.ncomp])

final.block.splsda = block.splsda(X, Y,
                                  ncomp=optimal.ncomp,
                                  keepX=optimal.keepX,
                                  design=design)

# 
plotIndiv(final.block.splsda, comp = c(1,2), # plot samples from final model
          ind.names = FALSE, # colour by class label
          ellipse = T, legend = TRUE, # include 95% confidence ellipse
          title = 'Block sPLS-DA, comp 1 & 2')


train <- sample(1:length(Y), 670*0.8) # randomly select 80% samples in training
test <- setdiff(1:length(Y), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- lapply(X, function(.).[train,])
X.test <- lapply(X, function(.) .[test,])
Y.train <- Y[train]
Y.test <- Y[test]


train.block.splsda = block.splsda(X.train, Y.train,
                                  near.zero.var = TRUE,
                                  ncomp=optimal.ncomp, 
                                  keepX=optimal.keepX)

predict.block.splsda <- predict(train.block.splsda, X.test, dist = "mahalanobis.dist")

# Get the balanced error rate. Usually improves as more components are used.

# evaluate the prediction accuracy for each modality
predict.cna <- predict.block.splsda$class$mahalanobis.dist$cna[,2]
confmat.cna = table(Y.test,Y.pred=factor(predict.cna, levels = levels(Y)))

predict.mut <- predict.block.splsda$class$mahalanobis.dist$mut[,2]
confmat.mut = table(Y.test,Y.pred=factor(predict.mut, levels = levels(Y)))

predict.rna <- predict.block.splsda$class$mahalanobis.dist$rna[,2]
confmat.rna = table(Y.test,Y.pred=factor(predict.rna, levels = levels(Y)))

predict.ig <- predict.block.splsda$class$mahalanobis.dist$ig[,2]
confmat.ig = table(Y.test,Y.pred=factor(predict.ig, levels = levels(Y)))

get.BER(confmat.cna)
get.BER(confmat.mut)
get.BER(confmat.ig)
get.BER(confmat.rna)

# plot AUROC curves
auc.block.splsda = auroc(final.block.splsda,
                   # newdata=X.test,
                   # outcome.test=Y.test, 
                   roc.block = 4,
                   roc.comp = 2,
                   print=F)

# plot heatmap
colors=unlist(list(SR="#FFFFFF",GHR="#00FF00",FHR="#FF0000")[optimal.splsda$Y])
cim(optimal.splsda$variates$X[,1:4],dist.method = c('euclidean','euclidean'),row.sideColors = colors,transpose = TRUE,cut.tree=c(0,0), save='png', name.save = './assets/splsda.rna.772.cim')


# plot individual genes that are important
plotVar(final.block.splsda)
