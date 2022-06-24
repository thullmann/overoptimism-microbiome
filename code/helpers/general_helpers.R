# The code for the functions norm_counts, trans_to_diss and trans_to_sim is taken from the NetCoMi package,
# available at https://github.com/stefpeschel/NetCoMi

norm_counts = function(countMat, normMethod, normParam = NULL) {
  
  if(normMethod %in% c("fractions", "TSS")){
    
    countMat_norm = t(apply(countMat, 1, function(x) x/sum(x)))
    
  } else if(normMethod == "VST"){
    
    countMat.tmp = countMat
    countMat = apply(countMat.tmp, 2, as.integer)
    rownames(countMat) = rownames(countMat.tmp)
    
    normParam$object = countMat
    countMat_norm = do.call("varianceStabilizingTransformation", normParam)
    
  } else if(normMethod == "clr"){
    
    normParam$x.f = countMat
    normParam$mar = 1
    if(is.null(normParam$base)) normParam$base = exp(1)
    
    countMat_norm = t(do.call(SpiecEasi::clr, normParam))
    
  } else if(normMethod == "mclr"){
    
    normParam$dat = countMat
    fun_mclr = get("mclr", asNamespace("SPRING"))
    
    countMat_norm = do.call(fun_mclr, normParam)
  }
  
  return(countMat_norm)
}

trans_to_diss = function(x, dissFunc) {
  xvec = x[lower.tri(x)]
  
  if(dissFunc == "signed"){
    dissvec = sqrt(0.5 * (1-xvec))
  } else if(dissFunc == "unsigned"){
    dissvec = sqrt(1-xvec^2)
  }
  dissMat = x
  dissMat[lower.tri(dissMat)] = dissvec
  dissMat[upper.tri(dissMat)] = t(dissMat)[upper.tri(t(dissMat))]
  diag(dissMat) = 0
  
  return(dissMat)
}

trans_to_sim = function(x) {
  
  x.tmp = x[upper.tri(x)]
  x.tmp = x[!is.infinite(x)]
  
  simMat = 1-x
  simMat[x == Inf] = 0
  
  return(simMat)
}

calc_jaccard = function(group1, group2) {
  uni = length(union(group1, group2))
  inters = length(intersect(group1, group2))
  jacc = ifelse(uni>0, inters/uni, 0)
}