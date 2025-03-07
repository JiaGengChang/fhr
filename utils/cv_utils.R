source('./utils/datasets.R')

# train and validate a single model
# fit on train data, and transform on validation data
cv = function(shuffle,fold,fn.fit,fn.transform,args){
  # fn.fit: input `args`, returns a trained_model
  # fn.fit - fits on train set and hyperparameters in args 
  # fn.transform: input trained_model and `args`, returns a prediction vector
  # fn.transform - transform on validation set in args
  # value: a 3x3 confusion matrix, truth are rows, preds are cols
  load_split_data(shuffle,fold)
  fitted.model = fn.fit(args)
  Y.pred = fn.transform(fitted.model, args)
  cm=get.confusion_matrix(args$Yvalid, predicted = Y.pred)
  class(cm) = c("confusion", class(cm))
  return(cm)
}

# wrapper for cross validation
# returns a nested list of CV results of n.shuffle shuffles and n.fold folds
# each list element is a confusion matrix
repeated.cv = function(n.shuffle=10,n.fold=5){
  cv.results = lapply(1:n.shuffle, function(shuffle){
    result.shuffle=lapply(1:n.fold, function(fold){
      result.fold = cv(shuffle,fold)
      class(result.fold) = c("fold",class(result.fold))
      return(result.fold)
    })
    class(result.shuffle) = c("shuffle",class(result.shuffle))
    return(result.shuffle)
  })
  class(cv.results) = c("experiment",class(cv.results))
  return(cv.results)
}
