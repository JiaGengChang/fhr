# get samples with both features and labels
subset_to_labelled = function(data, fhr){
  # identify the names of labeled patients
  common_ids = intersect(rownames(data), rownames(fhr))
  
  # subset to common ids
  common_data = data[rownames(data) %in% common_ids,]
  common_fhr = fhr[rownames(fhr) %in% common_ids,]
  
  # match the order of rows
  common_data = common_data[match(common_ids, rownames(common_data)),]
  common_fhr = common_fhr[match(common_ids, rownames(common_fhr)),]
  
  # convert integer-coded risk to factors
  common_fhr$risk = factor(common_fhr$risk,levels=c(0,1,2),labels=c("SR","GHR","FHR"))
  
  assign("common_ids",common_ids,envir=.GlobalEnv)
  assign("common_data",common_data,envir=.GlobalEnv)
  assign("common_fhr",common_fhr,envir=.GlobalEnv)
  
  return(0)
}
