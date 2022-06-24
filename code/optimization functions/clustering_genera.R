clustering_genera = function(df_discov, df_valid, nrsplit, n, k_fix = NULL) {
  
  ## applying different method combinations for generating clusterings ###########################################
  
  res = data.frame(measure = character(), zeroMethod = character(), normMethod = character(), 
                   sparsMethod = character(),
                   clustAlgo = character(), k = numeric(), linkage = character(),
                   ari = numeric())
  
  for (i in 1:nrow(assoc_distbased)) {
    
    if (assoc_distbased$normMethod[i] == "VST") {
      normPar = list(fitType = "local")
    } else {
      normPar = NULL
    }
    
    # first, construct dissimilarity matrix
    
    if (assoc_distbased$measure[i] == "latentcor") { # for latentcor (semi-parametric rank based), use mixedCCA
      
      # mixedCCA needs count matrix with samples in the rows, bacteria in the columns 
      
      if(df_discov@otu_table@taxa_are_rows) {
        countMat = t(df_discov@otu_table@.Data) 
      } else {
        countMat = df_discov@otu_table@.Data
      }
      # no zero treatment for mclr
      df_norm = norm_counts(countMat = countMat, normMethod = assoc_distbased$normMethod[i],
                            normParam = normPar)
      
      assoMat = mixedCCA::estimateR(df_norm, type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
      dissMat = trans_to_diss(assoMat, dissFunc = "signed")
      rm(countMat, df_norm, assoMat)
      
    } else { # for everything apart from latentcor, use NetCoMi
      dissMat = netConstruct(df_discov, measure = assoc_distbased$measure[i],
                             zeroMethod = assoc_distbased$zeroMethod[i],
                             normMethod = assoc_distbased$normMethod[i], normPar = normPar,
                             sparsMethod = "none", 
                             dissFunc = "signed", weighted = TRUE,
                             seed = 123456)$dissMat1
      dissMat[is.infinite(dissMat)] = sqrt(0.5) # NetCoMi assigns diss=infinity to OTUs with association=0
    }
    
    # second, cluster the data
    
    # hierarchical clustering
    clust = clustering(dissMat = dissMat, algo = "hierarchical", clustPar = list(k = k_fix, linkage = "average"))$cluster 
    #sanity check
    if(!all(names(clust) == names(df_discov@tax_table@.Data[,5]))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_distbased$measure[i], zeroMethod = assoc_distbased$zeroMethod[i], 
                                normMethod = assoc_distbased$normMethod[i], 
                                clustAlgo = "hierarchical", k = length(unique(clust)), linkage = "average", 
                                sparsMethod = "none",
                                ari = adjustedRandIndex(clust, as.numeric(as.factor(df_discov@tax_table@.Data[,5]))),
                                clust = I(list(clust))))
    
    # spectral clustering
    clust = clustering(adjaMat = trans_to_sim(dissMat), algo = "spectral", clustPar = list(k = k_fix))$cluster
    #sanity check
    if(!all(names(clust) == names(df_discov@tax_table@.Data[,5]))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_distbased$measure[i], zeroMethod = assoc_distbased$zeroMethod[i], 
                                normMethod = assoc_distbased$normMethod[i], 
                                clustAlgo = "spectral", k = length(unique(clust)), linkage = "none",
                                sparsMethod = "none",
                                ari = adjustedRandIndex(clust, as.numeric(as.factor(df_discov@tax_table@.Data[,5]))),
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
                      thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                      dissFunc = "signed", weighted = TRUE,
                      seed = 123456)
    assoMat = nc$assoMat1 
    adjaMat = nc$adjaMat1
    
    # second, cluster the data
    
    # fast greedy modularity clustering
    clust = clustering(adjaMat = adjaMat, algo = "cluster_fast_greedy", clustPar = list(k = NULL))$cluster
    #sanity check
    if(!all(names(clust) == names(df_discov@tax_table@.Data[,5]))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i], 
                                normMethod = assoc_networkbased$normMethod[i], 
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                clustAlgo = "cluster_fast_greedy", k = length(unique(clust)), linkage = "none",
                                ari = adjustedRandIndex(clust, as.numeric(as.factor(df_discov@tax_table@.Data[,5]))),
                                clust = I(list(clust))))
    
    # Louvain method 
    clust = clustering(adjaMat = adjaMat, algo = "cluster_louvain", clustPar = list(k = NULL))$cluster
    #sanity check
    if(!all(names(clust) == names(df_discov@tax_table@.Data[,5]))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i],
                                normMethod = assoc_networkbased$normMethod[i],
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                clustAlgo = "cluster_louvain", k = length(unique(clust)), linkage = "none",
                                ari = adjustedRandIndex(clust, as.numeric(as.factor(df_discov@tax_table@.Data[,5]))),
                                clust = I(list(clust))))
    
    
    # clustering with manta
    clust = clustering(assoMat = assoMat, algo = "manta", clustPar = list(filepath_s = nrsplit, filepath_i = i, filepath_n = n, k = NULL))$cluster
    #sanity check
    if(!all(names(clust) == names(df_discov@tax_table@.Data[,5]))) {
      print("names not equal")
    }
    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i],
                                normMethod = assoc_networkbased$normMethod[i],
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                clustAlgo = "manta", k = length(unique(clust)), linkage = "none",
                                ari = adjustedRandIndex(clust, as.numeric(as.factor(df_discov@tax_table@.Data[,5]))),
                                clust = I(list(clust))))
  }
  
  ## optimization ####################################################################################################

  best_result = which(res$ari == max(res$ari))[1]
  measure = res[best_result, "measure"]
  zeroMethod = res[best_result, "zeroMethod"]
  normMethod = res[best_result, "normMethod"]
  clustAlgo = res[best_result, "clustAlgo"]
  k = res[best_result, "k"]
  linkage = res[best_result, "linkage"]
  sparsMethod = res[best_result, "sparsMethod"]
  ari = res[best_result, "ari"]
  clustDiscov = res$clust[[best_result]]
  
  ## validation #####################################################################################################
  
  if (normMethod == "VST") {
    normPar = list(fitType = "local")
  } else {
    normPar = NULL
  }
  
  if (measure == "latentcor") { # for latentcor (semi-parametric rank based), use mixedCCA
    
    # mixedCCA needs count matrix with samples in the rows, bacteria in the columns 
    
    if(df_valid@otu_table@taxa_are_rows) {
      countMat = t(df_valid@otu_table@.Data) 
    } else {
      countMat = df_valid@otu_table@.Data
    }
    
    normMethod = "mclr"
    # no zero treatment for mclr
    df_norm = norm_counts(countMat = countMat, normMethod = normMethod,
                          normParam = normPar)
    
    assoMat = mixedCCA::estimateR(df_norm, type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
    dissMat = trans_to_diss(assoMat, dissFunc = "signed")
    
  } else { # for everything apart from latentcor, use NetCoMi
    nc = netConstruct(df_valid, measure = measure,
                      zeroMethod = zeroMethod,
                      normMethod = normMethod, normPar = normPar,
                      sparsMethod = sparsMethod, 
                      thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                      dissFunc = "signed", 
                      seed = 123456)
    dissMat = nc$dissMat1
    assoMat = nc$assoMat1
    adjaMat = nc$adjaMat1
  }
  
  if (clustAlgo == "hierarchical") {
    dissMat[is.infinite(dissMat)] = sqrt(0.5) #NetCoMi assigns diss=infinity to OTUs with association=0
    clustPar = list(k = k_fix, linkage = "average")
  } else if (clustAlgo == "spectral") {
    dissMat[is.infinite(dissMat)] = sqrt(0.5) #NetCoMi assigns diss=infinity to OTUs with association=0
    adjaMat = trans_to_sim(dissMat) # "adjacency matrix" with a_ij = 1-sqrt(0.5) if association between i,j is 0. In contrast, nc$adjaMat1 sets a_ij to 0. 
    clustPar = list(k = k_fix) 
  } else if (clustAlgo == "manta") {
    clustPar = list(filepath_s = nrsplit, filepath_i = "valid", k = NULL)
  } else if (clustAlgo == "cluster_louvain") {
    clustPar = list(k = NULL) 
  } else if (clustAlgo == "cluster_fast_greedy") {
    clustPar = list(k = NULL)
  }
  clust = clustering(dissMat = dissMat, assoMat = assoMat, adjaMat = adjaMat, algo = clustAlgo, clustPar = clustPar)$cluster
  
  ariValid = adjustedRandIndex(clust, as.numeric(as.factor(df_valid@tax_table@.Data[,5])))
  
  ariStab = adjustedRandIndex(clust, clustDiscov)
  
  val_result = data.frame(measure = measure, zeroMethod = zeroMethod, normMethod = normMethod, 
                          sparsMethod = sparsMethod,
                          clustAlgo = clustAlgo, kDiscov = k, kValid = length(unique(clust)), linkage = linkage,
                          ariDiscov = ari, ariValid = ariValid, ariStab = ariStab)
  return(list(res = res, val_result = val_result))

}

