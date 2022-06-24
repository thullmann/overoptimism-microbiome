# Hierarchical clustering ############################################################################

hier.cluster = function(dissMat, k, linkage = "average") {
  h_clust = hclust(as.dist(dissMat), method = linkage)
  cluster = cutree(h_clust, k = k)
  names(cluster) = dimnames(dissMat)[[1]]
  return(list("cluster" = cluster))
}


# Spectral clustering #################################################################################

# The following code for spectral clustering was written by Badri, Kurtz et al. (2020) and can be found 
# at https://github.com/MichelleBadri/NormCorr-manuscript

## Create affinity matrix from nearest neighbor similarity matrix
# kNN graph 
make.affinity = function(S, k=3, sym="or") {
  S = S-diag(diag(S))
  if (k >= ncol(S)) {
    return(S)
  } else {
    thr = apply(S, 1, function(s) s[order(s, decreasing=TRUE)[k]])
    A = S>=thr
    if (sym == "or") A = sign(A + t(A))
    else if (sym == "and") A = sign(A * t(A))
    return(S*A)
  }
}

## Normalized spectral clustering, Ng, Jordan & Weiss (2002)
spectral.cluster = function(simMat, k1 = 2, kmax = 16, k = NULL) {
  # R => Correlation/proportionality matrix
  # k1 => nearest neighbors
  
  A = make.affinity(simMat , k1, 'or')
  ## degree matrix
  d = apply(A, 1, sum)
  D = diag(d)
  D_isqrt = diag(1/sqrt(d))
  ## normalized, symmetric Laplacian
  L = D_isqrt %*% (D - A) %*% D_isqrt
  evL = eigen(L, symmetric=TRUE)
  spectrum = rev(evL$values)[1:(kmax-1)]
  zeroEvals = which(abs(spectrum) < 1e-9)
  # tmp = rle(diff(zeroEvals)) 
  # k2 = tmp$lengths[1]+1
  k2 = length(zeroEvals)
  if (!is.null(k)) k2 = k 
  if (k2 != 1) {
    T  = evL$vectors[,1:k2] # eigenvectors of the k2 largest eigenvalues
    # T  = evL$vectors[,1:k]
    T  = T/sqrt(rowSums(T^2))
    # define 0/0 := 0
    T[!is.finite(T)] = 0
    set.seed(1000)
    km = kmeans(T, centers=k2, nstart=3000)
    # km = kmeans(T, centers=k, nstart=3000)
    
    ## Return eigendecomposition for spectrum plots
    # km$Leig = evL
    # km$spectrum = spectrum
    cluster = km$cluster
    
  } else {
    cluster = rep(1, times = dim(simMat)[1])
  }
  names(cluster) = dimnames(simMat)[[1]]
  
  return(list("cluster" = cluster))
}


# fast greedy modularity optimization and the Louvain method for community detection #################

igraph.cluster = function(adjaMat, clustMethod, k = NULL) {
  net = igraph::graph_from_adjacency_matrix(adjaMat, weighted=TRUE,
                                            mode="undirected", diag=FALSE)
  clust = do.call(getExportedValue("igraph", clustMethod), 
                  c(list(net)))
  cluster = clust$membership
  if (!is.null(k)) if(length(unique(cluster)) < k) cluster = cut_at(clust, no = k)
  names(cluster) = clust$names
  return(list("cluster" = cluster))
}


# manta ##############################################################################################

manta = function(assoMat, filepath_s, filepath_i, filepath_n, k = NULL) {
  
  # filepath_s indicates number of current split in discov/valid data
  # filepath_i indicates number of current iteration over 1:nrow(assoc_networkbased)
  # filepath_n indicates number of samples drawn
  
  # use igraph::write_graph() and igraph::read_graph() for converting to/from graphml (input and output of manta)
  net = igraph::graph_from_adjacency_matrix(assoMat, weighted=TRUE,
                                            mode="undirected", diag=FALSE)
  
  file_input = paste0("results/mantaGraphml/", n, "/split", filepath_s, "i", filepath_i, "net.graphml")
  file_output = paste0("results/mantaGraphml/", n, "/split", filepath_s, "i", filepath_i, "manta_output")
  write_graph(net, file = file_input, format = "graphml")
  
  if (!is.null(k)) {
    system(paste("manta -i", file_input, "-o", file_output, "-f graphml", "-min", k, "-seed 123456", sep = " "))
  } else {
    system(paste("manta -i", file_input, "-o", file_output, "-f graphml -seed 123456", sep = " "))
  }
  net = read_graph(paste0(file_output, ".graphml"), format = "graphml")
  cluster = V(net)$cluster
  names(cluster) = V(net)$name
  return(list("cluster" = cluster))
  
}


# function for accessing all cluster algorithms ####################################################

clustering = function(dissMat = NULL, assoMat = NULL, adjaMat = NULL, algo = "hierarchical", clustPar = NULL) {
  if (algo == "hierarchical") {
    clustPar$dissMat = dissMat
    cluster = do.call("hier.cluster", clustPar)
    # clustPar must contain k = k, linkage = "average"
  } else if (algo == "spectral") {
    clustPar$simMat = adjaMat
    cluster = do.call("spectral.cluster", clustPar)
  } else if (algo == "cluster_fast_greedy" | algo == "cluster_louvain") {
    clustPar$adjaMat = adjaMat
    clustPar$clustMethod = algo
    cluster = do.call("igraph.cluster", clustPar)
  } else if (algo == "manta") {
    clustPar$assoMat = assoMat
    cluster = do.call("manta", clustPar)
  }
}

