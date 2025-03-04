
# inner loop of cross validation
cv.helper = function(ncomp, keepX, shuffle,fold){
  train=read.delim(sprintf('./data/splits/%d/train_%d.txt',shuffle,fold))$PUBLIC_ID
  test=read.delim(sprintf('./data/splits/%d/valid_%d.txt',shuffle,fold))$PUBLIC_ID
  Y.train=common.data.label[train,'risk']
  Y.test=common.data.label[test,'risk']

  if (length(names(keepX))>1){
    # perform BLOCK splsda
    X.train = lapply(X, function(x) x[train,])
    X.test = lapply(X, function(x) x[test,])
    train.splsda = block.plsda(X.train, Y.train, ncomp=ncomp, near.zero.var = T)
    predict.splsda = predict(train.splsda, X.test, dist='mahalanobis.dist')
    Y.pred=predict.splsda$WeightedVote$mahalanobis.dist[,ncomp]
  } else {
    X.train=X[train,]
    X.test=X[test,]
    train.splsda = splsda(X.train, Y.train, ncomp=ncomp, keepX = keepX)
    predict.splsda <- predict(train.splsda, X.test, dist = "mahalanobis.dist")
    Y.pred=predict.splsda$class$mahalanobis.dist[,ncomp]
  }
  cm=get.confusion_matrix(Y.test,predicted = Y.pred)
  class(cm) = c("confusion", class(cm))
  return(cm)
}

# wrapper for cross validation
cv = function(ncomp, keepX, n_shuffles=10,n_folds=5){
  cv.results = lapply(1:n_shuffles, function(shuffle){
    result.shuffle=lapply(1:n_folds, function(fold){
      result.fold = cv.helper(ncomp, keepX, shuffle,fold)
      class(result.fold) = c("fold",class(result.fold))
      return(result.fold)
    })
    class(result.shuffle) = c("shuffle",class(result.shuffle))
    return(result.shuffle)
  })
  class(cv.results) = c("mixo_splsda_experiment",class(cv.results))
  return(cv.results)
}


# Function to calculate precision, accuracy, and F1 score from a confusion matrix
calculate_metrics = function(cm){
  precision = diag(cm)/colSums(cm)
  recall = diag(cm)/rowSums(cm)
  f1_score = ifelse(precision+recall==0, 0, 2 * (precision * recall) / (precision + recall))
  f1_score[is.na(f1_score)] = 0
  return(list(precision=precision,
              recall=recall,
              f1_score=f1_score))
}

# calculate summary statistics - precision, recall, f1 score
summary.mixo_splsda_experiment = function(experiment) {
    # metrics for each fold
    results = lapply(unlist(experiment,recursive = F), calculate_metrics)
    
    # mean and standard deviation of each metric
    precision.raw = do.call(rbind, lapply(results,function(x) x$precision))
    precision.mean = apply(precision.raw, 2, mean)
    precision.std = apply(precision.raw, 2, function(x)sd(x))
    
    recall.raw = do.call(rbind, lapply(results,function(x) x$recall))
    recall.mean = apply(recall.raw, 2, mean)
    recall.std = apply(recall.raw, 2, function(x)sd(x))
    
    f1_score.raw = do.call(rbind, lapply(results,function(x) x$f1_score))
    f1_score.mean = apply(f1_score.raw, 2, mean)
    f1_score.std = apply(f1_score.raw, 2, function(x)sd(x))
    
    result = list(precision.mean=precision.mean,
                precision.std=precision.std,
                recall.mean=recall.mean,
                recall.std=recall.std,
                f1_score.mean=f1_score.mean,
                f1_score.std=f1_score.std)
    class(result) = c("mixo_splsda_experiment_summary", class(result))
    return(result)
}

print.mixo_splsda_experiment_summary = function(experiment.summary){
  result = cbind(experiment.summary$precision.mean,
                 experiment.summary$precision.std,
                 experiment.summary$recall.mean,
                 experiment.summary$recall.std,
                 experiment.summary$f1_score.mean,
                 experiment.summary$f1_score.std)
  colnames(result) = c("precision.mu","precision.sd",
                       "recall.mu","recall.sd",
                       "f1.mu","f1.sd")
  rownames(result) = c("SR","GHR","FHR")
  print(result)
  return(invisible(result))
}

# Example usage
# my.experiment = cv()
# my.summary = summary(my.experiment)
# Example output
# > print(my.summary)
# precision.mu precision.sd recall.mu  recall.sd     f1.mu      f1.sd
# SR     0.7925033   0.03356358 0.7488889 0.03828501 0.7695301 0.02958217
# GHR    0.8034909   0.05468321 0.6691923 0.07026339 0.7280691 0.05195984
# FHR    0.1880977   0.06071247 0.3441758 0.12540096 0.2415945 0.07930788
