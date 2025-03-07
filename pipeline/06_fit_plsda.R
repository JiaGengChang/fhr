setwd('/home/users/nus/e1083772/fhr/')
source('./utils/datasets.R')
source('./utils/cv_utils.R')
suppressPackageStartupMessages(library(mixOmics))
suppressPackageStartupMessages(library(smotefamily))

args = commandArgs(trailingOnly = TRUE)
shuffle = as.integer(args[1])
fold = as.integer(args[2])


ids = load_split_data(shuffle,fold)

# delete `all` tables
suppressMessages(rm(all.rna, all.label, all.cn.arm, all.cn.seg, all.ig, all.mut))
suppressMessages(gc())

# FEATURE SELECTION
# subset to variable genes
train.rna.vars = apply(train.rna, 2, var)
var.keep = train.rna.vars > 2
print(table(var.keep))

do_smote=FALSE

# triple the number of FHR classes
if (do_smote){
  genData = SMOTE(data.frame(train.rna[,var.keep]),train.label$risk, dup_size = 5, K=5)
  ncD = ncol(genData$data)
  X = genData$data[,-ncD]
  Y = genData$data[,ncD]
  } else {
    X = data.frame(train.rna[,var.keep])
    Y = train.label$risk
  }

# visualize number of classes
print(table(Y))

# args for the model training
args = list(exptname="RNA_HVG_noSMOTE",
            X=X,
            Y=Y,
            Xvalid=data.frame(valid.rna[,var.keep]),
            Yvalid=valid.label$risk,
            ncomp=3,
            keepX=c(seq(1,10,1),seq(200,1000,100)),
            measure='BER',
            validation='Mfold',
            folds=5,
            nrepeat=5)

# tune the sparsity of inputs
print(args$keepX)

# one function per model
fn.fit = function(args){
  require(mixOmics)
  # cross validation to tune n.components and n.features
  tune.result = tune.splsda(args$X, 
                            args$Y, 
                            ncomp=args$ncomp,
                            test.keepX = args$keepX, 
                            measure = args$measure,
                            validation = args$validation,
                            folds = args$folds,
                            nrepeat = args$nrepeat)
  # train full model
  choice.keepX = tune.result$choice.keepX
  choice.ncomp = tune.result$choice.ncomp$ncomp
  print(choice.keepX)
  sprintf('Chosen number of components: %d',choice.ncomp)
  # ncomp may be NULL if looCV or nrepeat too low
  choice.ncomp = ifelse(is.null(choice.ncomp),2,choice.ncomp)
  fitted.model = splsda(args$X, args$Y, ncomp=choice.ncomp, keepX=choice.keepX)
  return(fitted.model)
}

fn.transform = function(fitted.model, args){
  prediction = predict(fitted.model, newdata=args$Xvalid, dist='mahalanobis.dist')
  print(dim(args$Xvalid))
  Y.pred = prediction$class$mahalanobis.dist
  # print(Y.pred)
  Y.pred = Y.pred[,ncol(Y.pred)] # use the last component
  return(Y.pred)
}

# confusion matrix of validation examples
cm = cv(shuffle, fold, fn.fit, fn.transform, args)
print(cm)

# save confusion matrix
outfile=sprintf('./scores/models/sPLSDA/%s/valid-%d-%d.csv',args$exptname,shuffle,fold)
if (!dir.exists(dirname(outfile))) {
  dir.create(dirname(outfile), recursive = TRUE)
}
write.csv(cm,outfile,quote=F)
