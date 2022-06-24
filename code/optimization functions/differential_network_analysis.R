differential_network_analysis = function(df_discov_no_antibiotics, df_discov_antibiotics, 
                             df_valid_no_antibiotics, df_valid_antibiotics, nrsplit, n) {
  
  ## applying different method combinations for network construction ###########################################
  
  res = data.frame(measure = character(), zeroMethod = character(), normMethod = character(), 
                   sparsMethod = character(),
                   gcd = numeric())

  for (i in 1:nrow(assoc_networkbased)) {
    if (assoc_networkbased$normMethod[i] == "VST") {
      normPar = list(fitType = "local")
    } else {
      normPar = NULL
    }
    
    nc = netConstruct(data = df_discov_no_antibiotics, 
                      data2 = df_discov_antibiotics,
                      measure = assoc_networkbased$measure[i],
                      zeroMethod = assoc_networkbased$zeroMethod[i],
                      normMethod = assoc_networkbased$normMethod[i], normPar = normPar,
                      sparsMethod = assoc_networkbased$sparsMethod[i], 
                      thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                      dissFunc = "signed", weighted = TRUE,
                      seed = 123456)
    
    adjaMat1 = nc$adjaMat1
    adjaMat2 = nc$adjaMat2
    
    net1 = igraph::graph_from_adjacency_matrix(adjaMat1, weighted=T,
                                               mode="undirected", diag=F)
    net2 = igraph::graph_from_adjacency_matrix(adjaMat2, weighted=T,
                                                        mode="undirected", diag=F)
    
    # calculate GCD
    
    edgelist1 = igraph::get.edgelist(net1, names = FALSE)
    edgelist2 = igraph::get.edgelist(net2, names = FALSE)
    edgelist1 = apply(edgelist1, 2, as.integer)
    edgelist2 = apply(edgelist2, 2, as.integer)
    orbitCounts1 = orca::count4(edgelist1)[,c(1:3,5:12)] #  only non-redundant orbits, i.e. the fourth orbit (which is orbit "O3", orbit enumeration starts at 0) is removed 
    orbitCounts2 = orca::count4(edgelist2)[,c(1:3,5:12)]
    
    buffer1 = matrix(0, nrow=nrow(adjaMat1)-nrow(orbitCounts1), ncol=ncol(orbitCounts1))
    orbitCounts1 = rbind(orbitCounts1, buffer1)
    buffer2 = matrix(0, nrow=nrow(adjaMat2)-nrow(orbitCounts2), ncol=ncol(orbitCounts2))
    orbitCounts2 = rbind(orbitCounts2, buffer2)
    
    cor1 = cor(orbitCounts1, method = "spearman")
    cor2 = cor(orbitCounts2, method = "spearman")
    cor1 = cor1[upper.tri(cor1)] 
    cor2 = cor2[upper.tri(cor2)] 
    gcd = sqrt(sum((cor1 - cor2)^2))
      
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i], 
                                normMethod = assoc_networkbased$normMethod[i], 
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                gcd = gcd))
    
  }
  
  ## optimization ####################################################################################################
  
  best_result = which(res$gcd == max(res$gcd))[1]
  # best_result = which(res$ari == min(res$ari))[1]
  # best_result = which(res$connectDiff == max(res$connectDiff))[1]
  measure = res[best_result,"measure"]
  zeroMethod = res[best_result, "zeroMethod"]
  normMethod = res[best_result, "normMethod"]
  sparsMethod = res[best_result, "sparsMethod"]
  connectDiffDiscov = res[best_result, "connectDiff"]
  ariDiscov = res[best_result, "ari"]
  clustDiscov1 = res$clust_no_antibiotics[[best_result]]
  clustDiscov2 = res$clust_antibiotics[[best_result]]
  gcdDiscov = res[best_result, "gcd"]
  
  ## validation #####################################################################################################
  
  if (normMethod == "VST") {
    normPar = list(fitType = "local")
  } else {
    normPar = NULL
  }
  
  nc = netConstruct(data = df_valid_no_antibiotics, 
                    data2 = df_valid_antibiotics, 
                    measure = measure,
                    zeroMethod = zeroMethod,
                    normMethod = normMethod, normPar = normPar,
                    sparsMethod = sparsMethod, 
                    thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                    dissFunc = "signed", 
                    seed = 123456)

  adjaMat1 = nc$adjaMat1
  adjaMat2 = nc$adjaMat2

  net1 = igraph::graph_from_adjacency_matrix(adjaMat1, weighted=T,
                                             mode="undirected", diag=F)
  net2 = igraph::graph_from_adjacency_matrix(adjaMat2, weighted=T,
                                             mode="undirected", diag=F)
  
  # calculate GCD
  
  edgelist1 = igraph::get.edgelist(net1, names = FALSE)
  edgelist2 = igraph::get.edgelist(net2, names = FALSE)
  edgelist1 = apply(edgelist1, 2, as.integer)
  edgelist2 = apply(edgelist2, 2, as.integer)
  orbitCounts1 = orca::count4(edgelist1)[,c(1:3,5:12)] # only non-redundant orbits
  orbitCounts2 = orca::count4(edgelist2)[,c(1:3,5:12)]
  cor1 = cor(orbitCounts1, method = "spearman")
  cor2 = cor(orbitCounts2, method = "spearman")
  cor1 = cor1[upper.tri(cor1)] 
  cor2 = cor2[upper.tri(cor2)] 
  gcdValid = sqrt(sum((cor1 - cor2)^2))
  
  val_result = data.frame(measure = measure, zeroMethod = zeroMethod, normMethod = normMethod, 
                          sparsMethod = sparsMethod,
                          gcdDiscov = gcdDiscov, gcdValid = gcdValid)
  return(list(res = res, val_result = val_result))
  
}
