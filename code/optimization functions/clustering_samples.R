clustering_samples = function(df_discov, df_valid, nrsplit, n, k_fix = NULL) {
  
  ## applying different method combinations for generating clusterings ###########################################
  
  res = data.frame(measure = character(), zeroMethod = character(), normMethod = character(), 
                   sparsMethod = character(),
                   clustAlgo = character(), k = numeric(), linkage = character(),
                   sw = numeric())
  
  # DMM clustering
  # requires count matrix with samples in the rows, bacteria in the columns 
  
  if (df_discov@otu_table@taxa_are_rows) {
      countMat = t(df_discov@otu_table@.Data)
    } else {
      countMat = df_discov@otu_table@.Data
  }
  
  clustdmm = lapply(2:10, dmn, count = countMat, verbose=TRUE)
  lplc = sapply(clustdmm, laplace)
  best = clustdmm[[which.min(lplc)]]
  clust = mixture(best, assign = TRUE)
  
  dissMat = netConstruct(df_discov, measure = "bray",
                    zeroMethod = "none",
                    normMethod = "mclr", normPar = NULL,
                    sparsMethod = "none", 
                    weighted = TRUE,
                    seed = 123456)$dissMat1
  
  res = rbind(res, data.frame(measure = "countMat", zeroMethod = "none",
                              normMethod = "none",
                              clustAlgo = "dmm", k = length(unique(clust)), linkage = "none",
                              sparsMethod = "none",
                              sw = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE),
                              clust = I(list(clust))))
  
  for (i in 1:nrow(assoc_distbased)) {
    
    if (assoc_distbased$normMethod[i] == "VST") {
      normPar = list(fitType = "local")
    } else {
      normPar = NULL
    }
    
    # first, construct dissimilarity matrix
    
    nc = netConstruct(df_discov, measure = assoc_distbased$measure[i],
                             zeroMethod = assoc_distbased$zeroMethod[i],
                             normMethod = assoc_distbased$normMethod[i], normPar = normPar,
                             sparsMethod = "none", 
                             weighted = TRUE,
                             seed = 123456)
    dissMat = nc$dissMat1
    adjaMat = nc$adjaMat1
    
    # second, cluster the data
    
    # PAM clustering
    if (is.null(k_fix)) {
      avg_silwidth = numeric()
      for (k in 2:10) {
        clust = clustering(dissMat = dissMat, algo = "pam", clustPar = list(k = k))$cluster
        avg_silwidth[k-1] = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE)
      }
      best_k = which(avg_silwidth == max(avg_silwidth))[1] + 1
      clust = clustering(dissMat = dissMat, algo = "pam", clustPar = list(k = best_k))$cluster
    } else {
      clust = clustering(dissMat = dissMat, algo = "pam", clustPar = list(k = k_fix))$cluster 
    }
    
    #sanity check
    if(!all(names(clust) == colnames(df_discov@otu_table@.Data))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_distbased$measure[i], zeroMethod = assoc_distbased$zeroMethod[i], 
                                normMethod = assoc_distbased$normMethod[i], 
                                clustAlgo = "pam", k = length(unique(clust)), linkage = "none", 
                                sparsMethod = "none",
                                sw = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE),
                                clust = I(list(clust))))
    
    
    # spectral clustering
    if (is.null(k_fix)) {
      avg_silwidth = numeric()
      for (k in 2:10) {
        clust = clustering(adjaMat = adjaMat, algo = "spectral", clustPar = list(k = k))$cluster
        avg_silwidth[k-1] = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE)
      }
      best_k = which(avg_silwidth == max(avg_silwidth))[1] + 1
      clust = clustering(adjaMat = adjaMat, algo = "spectral", clustPar = list(k = best_k))$cluster
    } else {
      clust = clustering(adjaMat = adjaMat, algo = "spectral", clustPar = list(k = k_fix))$cluster 
    }
    
    #sanity check
    if(!all(names(clust) == colnames(df_discov@otu_table@.Data))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_distbased$measure[i], zeroMethod = assoc_distbased$zeroMethod[i], 
                                normMethod = assoc_distbased$normMethod[i], 
                                clustAlgo = "spectral", k = length(unique(clust)), linkage = "none",
                                sparsMethod = "none",
                                sw = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE),
                                clust = I(list(clust))))
  }
  
  for (i in 1:nrow(assoc_networkbased)) {
    
    if (assoc_networkbased$normMethod[i] == "VST") {
      normPar = list(fitType = "local")
    } else {
      normPar = NULL
    }
    
    # first, construct association/adjacency matrix
    
    nc = netConstruct(df_discov, measure = assoc_networkbased$measure[i],
                      zeroMethod = assoc_networkbased$zeroMethod[i],
                      normMethod = assoc_networkbased$normMethod[i], normPar = normPar,
                      sparsMethod = assoc_networkbased$sparsMethod[i], 
                      thresh = 0.85, kNeighbor = 3, knnMutual = FALSE,
                      weighted = TRUE,
                      seed = 123456)
    dissMat = nc$dissMat1
    # NetCoMi assigns "inf" to sparsified distances. Set these values to 1 ("maximum" distance). 
    dissMat[is.infinite(dissMat)] = 1
    adjaMat = nc$adjaMat1
    
    # second, cluster the data
    
    # fast greedy modularity clustering
    clust = clustering(adjaMat = adjaMat, algo = "cluster_fast_greedy", clustPar = list(k = NULL))$cluster
    #sanity check
    if(!all(names(clust) == colnames(df_discov@otu_table@.Data))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i], 
                                normMethod = assoc_networkbased$normMethod[i], 
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                clustAlgo = "cluster_fast_greedy", k = length(unique(clust)), linkage = "none",
                                sw = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE),
                                clust = I(list(clust))))
    
    # Louvain method 
    clust = clustering(adjaMat = adjaMat, algo = "cluster_louvain", clustPar = list(k = NULL))$cluster
    #sanity check
    if(!all(names(clust) == colnames(df_discov@otu_table@.Data))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i],
                                normMethod = assoc_networkbased$normMethod[i],
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                clustAlgo = "cluster_louvain", k = length(unique(clust)), linkage = "none",
                                sw = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE),
                                clust = I(list(clust))))
    
    
  }
  
  ## optimization ####################################################################################################
  
  best_result = which(res$sw == max(res$sw))[1]
  measure = res[best_result, "measure"]
  zeroMethod = res[best_result, "zeroMethod"]
  normMethod = res[best_result, "normMethod"]
  clustAlgo = res[best_result, "clustAlgo"]
  k = res[best_result, "k"]
  linkage = res[best_result, "linkage"]
  sparsMethod = res[best_result, "sparsMethod"]
  sw = res[best_result, "sw"]
  clustDiscov = res$clust[[best_result]]
  
  ## validation #####################################################################################################
  
  if (normMethod == "VST") {
    normPar = list(fitType = "local")
  } else {
    normPar = NULL
  }
  
  if (clustAlgo == "dmm") { 
    if (df_valid@otu_table@taxa_are_rows) {
      countMat = t(df_valid@otu_table@.Data)
    } else {
      countMat = df_valid@otu_table@.Data
    }
    
    clustdmm = lapply(2:10, dmn, count = countMat, verbose=TRUE)
    lplc = sapply(clustdmm, laplace)
    best = clustdmm[[which.min(lplc)]]
    clust = mixture(best, assign = TRUE)
    
    dissMat = netConstruct(df_valid, measure = "bray",
                           zeroMethod = "none",
                           normMethod = "mclr", normPar = NULL,
                           sparsMethod = "none", 
                           weighted = TRUE,
                           seed = 123456)$dissMat1
    
  } else { # for everything apart from DMM, use NetCoMi
    nc = netConstruct(df_valid, measure = measure,
                      zeroMethod = zeroMethod,
                      normMethod = normMethod, normPar = normPar,
                      sparsMethod = sparsMethod, 
                      thresh = 0.85, kNeighbor = 3, knnMutual = FALSE,
                      seed = 123456)
    dissMat = nc$dissMat1
    # NetCoMi assigns "inf" to sparsified distances. Set these values to 1 ("maximum" distance). 
    dissMat[is.infinite(dissMat)] = 1
    adjaMat = nc$adjaMat1
    
    if (clustAlgo == "hierarchical") {
      clustPar = list(k = k, linkage = "average")
    } else if (clustAlgo == "pam") {
      clustPar = list(k = k)
    } else if (clustAlgo == "spectral") {
      clustPar = list(k = k) 
    } else if (clustAlgo == "cluster_louvain") {
      clustPar = list(k = NULL) 
    } else if (clustAlgo == "cluster_fast_greedy") {
      clustPar = list(k = NULL)
    }
    clust = clustering(dissMat = dissMat, adjaMat = adjaMat, algo = clustAlgo, clustPar = clustPar)$cluster
  }
  
  swValid = mean(cluster::silhouette(clust, dissMat)[,3], na.rm = TRUE)
  
  val_result = data.frame(measure = measure, zeroMethod = zeroMethod, normMethod = normMethod, 
                          sparsMethod = sparsMethod,
                          clustAlgo = clustAlgo, kDiscov = k, kValid = length(unique(clust)), linkage = linkage,
                          swDiscov = sw, swValid = swValid)
  return(list(res = res, val_result = val_result))
  
}