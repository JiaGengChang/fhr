# convert segment mean to integer copy number
seg_mean_to_cn = function(x){
  ifelse(x > 0.66, 2,
         ifelse(x > 0.38, 1,
                ifelse(x > -0.5, 0,
                       ifelse(x > -1.5, -1, -2))))
}

# convert a segment mean matrix to a integer copy number matrix
# 0 means normal, -1, -2 are deletions, +1, +2 are gains
get_cn_df = function(...) {
  out = apply(apply(...,2, seg_mean_to_cn),2,as.integer)
  rownames(out) = rownames(...)
  return(out)
}

# 2 means normal, 1, 0 are deletions, 3, 4 are gains 
get_cn2_df = function(...) {
  out = apply(apply(...,2, seg_mean_to_cn) + 2,2,as.integer)
  rownames(out) = rownames(...)
  return(out)
}
