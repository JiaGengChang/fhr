# functions to summarise score across a repeated cross validation experiment

score_repeated_cv = function(scoredir,n_shuffle=10,n_fold=5) {
  # scoredir - a directory containing csv files of confusion matrices
  # with filenames in the format valid-[1-10]-[1-5].csv
  
  # value: precision and recall statistics

  m.precision = numeric()
  m.recall = numeric()
  
  for (shuffle in 1:n_shuffle){
    for (fold in 1:n_fold){
      cm = as.matrix(read.csv(sprintf('%s/valid-%d-%d.csv',scoredir,shuffle,fold),row.names=1))
      precision = diag(cm)/colSums(cm)
      recall = diag(cm)/rowSums(cm)
      m.precision = rbind(m.precision,precision)
      m.recall = rbind(m.recall, recall)
    }
  }
  rownames(m.precision) = 1:nrow(m.precision)
  rownames(m.recall) = 1:nrow(m.recall)
  
  precision.mean = apply(m.precision, 2, mean)
  precision.sd = apply(m.precision, 2, sd)
  
  recall.mean = apply(m.recall, 2, mean)
  recall.sd = apply(m.recall, 2, sd)
  
  result = list(precision.mean=precision.mean,
                precision.sd=precision.sd,
                recall.mean=recall.mean,
                recall.sd=recall.sd,
                m.precision=m.precision, # raw values
                m.recall=m.recall, # raw values,
                scoredir=scoredir
                )
  class(result) = c("experiment", class(result))
  return (result)
}

# Print mean and stdev of precision and recall
summary.experiment = function(experiment){
  print(sprintf('Experiment results for %s:',experiment$scoredir))
  print('Precision')
  print(sprintf('%s: %.2f (%.3f)', names(experiment$precision.mean), experiment$precision.mean, experiment$precision.sd))
  print('Recall')
  print(sprintf('%s: %.2f (%.3f)', names(experiment$recall.mean), experiment$recall.mean, experiment$recall.sd))
  return(invisible(experiment))
}

## Example usage
# scoredir="./scores/models/sPLSDA/RNA"
# result = score_repeated_cv(scoredir,n_shuffle=10,n_fold=5)
# summary(result)
# [1] "Experiment results for ./scores/models/sPLSDA/RNA:"
# [1] "Precision"
# [1] "predicted.as.SR: 0.77 (0.029)"  "predicted.as.GHR: 0.79 (0.073)" "predicted.as.FHR: 0.21 (0.116)"
# [1] "Recall"
# [1] "SR: 0.80 (0.049)"  "GHR: 0.64 (0.085)" "FHR: 0.24 (0.123)"
