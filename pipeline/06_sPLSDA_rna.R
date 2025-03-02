library(mixOmics)
source('./utils/datasets.R')

load_common_data()

X = list(rna=common.data.rna, cna=common.data.cn, mut=common.data.mut, ig=common.data.ig)
Y = model.matrix(~0+., data=common.data.label)

design = matrix(1, ncol=length(X), nrow=length(X), dimnames=list(names(X),names(X)))
diag(design) = 0

ncomp=c(2)

list.keepX = list(rna=c(5, 15, 30, 67),
                  cna=c(5, 15, 42),
                  mut=c(5, 15, 32),
                  ig=c(1, 2, 3))

results.block.plsda = block.plsda(X, Y=common.data.label$risk, ncomp=ncomp, design='full')
plotIndiv(results.block.plsda, group=common.data.label$risk, ind.names = T)
plotVar(results.block.plsda)

BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)

tune.block.splsda.results = 
  tune.block.splsda(X, 
                    Y=common.data.label$risk, 
                    ncomp=2,
                    test.keepX=list.keepX,
                    design = 'full',
                    validation="Mfold",
                    folds=5,
                    nrepeat=2,
                    BPPARAM=BPPARAM)

plot(tune.block.splsda.results, col = color.jet(2)) 

X = common.data.rna
Y = common.data.label$risk
tune.splsda.results = tune.splsda(X, Y, ncomp=4,
                                  test.keepX = c(1:10, seq(20, 100, 10), seq(200, 700, 100), 772),
                                  validation="Mfold",
                                  folds=5,
                                  nrepeat=10,
                                  cpus=3)
plot(tune.splsda.results, col = color.jet(4))
tune.splsda.results$choice.ncomp
tune.splsda.results$choice.keepX


optimal.splsda = splsda(X, Y, 
                        ncomp=tune.splsda.results$choice.ncomp$ncomp, 
                        keepX=tune.splsda.results$choice.keepX)

perf.optimal.splsda = perf(optimal.splsda,
                           folds=5,
                           nrepeat=10,
                           auc=TRUE,
                           progressBar=TRUE,
                           cpus=3)

# visualize classification error rates
plot(perf.optimal.splsda)

# splsda.rna.772.indiv
plotIndiv(optimal.splsda, comp = c(1,2), # plot samples from final model
          ind.names = FALSE, # colour by class label
          ellipse = TRUE, legend = TRUE, # include 95% confidence ellipse
          title = 'sPLS-DA on 772 gene expression, comp 1 & 2')


train <- sample(1:nrow(X), 670*0.8) # randomly select 80% samples in training
test <- setdiff(1:nrow(X), train) # rest is part of the test set

# store matrices into training and test set:
X.train <- X[train,]
X.test <- X[test,]
Y.train <- Y[train]
Y.test <- Y[test]


train.splsda = splsda(X.train, Y.train, 
                      ncomp=tune.splsda.results$choice.ncomp$ncomp, 
                      keepX=tune.splsda.results$choice.keepX)

predict.splsda <- predict(train.splsda, X.test, dist = "all")

# Get the balanced error rate. Usually improves as more components are used.

# evaluate the prediction accuracy for the first two components
predict.comp2 <- predict.splsda$class$mahalanobis.dist[,2]
get.BER(table(Y.test,Y.pred=factor(predict.comp2, levels = levels(Y))))

## -------------------------------------------------------------------------------------------------------------------
# evaluate the prediction accuracy for the first three components
predict.comp3 <- predict.splsda$class$mahalanobis.dist[,3]
get.BER(table(Y.test, Y.pred=factor(predict.comp3, levels = levels(Y))))

## -------------------------------------------------------------------------------------------------------------------
# evaluate the prediction accuracy for all four components
predict.comp4 <- predict.splsda$class$mahalanobis.dist[,4]
get.BER(table(Y.test, Y.pred=factor(predict.comp4, levels = levels(Y))))


# plot AUROC curves
auc.splsda = auroc(train.splsda, 
                   newdata=X.test, 
                   outcome.test=Y.test, 
                   roc.comp = 4, 
                   print=T)

# plot heatmap
colors=unlist(list(SR="#FFFFFF",GHR="#00FF00",FHR="#FF0000")[optimal.splsda$Y])
cim(optimal.splsda$variates$X[,1:4],dist.method = c('euclidean','euclidean'),row.sideColors = colors,transpose = TRUE,cut.tree=c(0,0), save='png', name.save = './assets/splsda.rna.772.cim')


# plot individual genes that are important

