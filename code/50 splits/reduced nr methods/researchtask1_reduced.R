########################################################################################################
### Clustering of bacterial genera #####################################################################
### 50 repetitions of optimization on the discovery data and validation on the validation data #########
########################################################################################################

library(NetCoMi)
library(phyloseq)
library(mixedCCA)
library(igraph)
library(mclust)
library(readxl)
library(parallel)
library(dplyr)

library(reticulate)
# If you use reticulate for the first time, check whether Python is installed: 
# py_available(initialize=TRUE)
# py_discover_config()

# install the manta module  
#py_install("git")
#py_install("git+https://github.com/ramellose/manta.git", pip = TRUE)

system("echo $PATH") # path might be without python interpreter
try(reticulate::py_run_file("hello.py")) # throws error but necessary to add interpreter path to PATH variable
system("echo $PATH") # should now have additional python path
py_run_string("import os")
py_run_string("os.environ['OPENBLAS_NUM_THREADS'] = '1'")

# two lists of association methods: 
# 1.) input for distance-based cluster algorithms (without sparsification),
# 2.) input for network-based cluster algorithms (with sparsification)
assoc_distbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "task1_distbased")
assoc_networkbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "task1_networkbased")

source("code/helpers/general_helpers.R")
source("code/helpers/cluster_algorithms.R")
source("code/optimization functions/clustering_genera.R")

ag.genus = readRDS("data/ag.genus.rds")

for (n in c(100, 250, 500, 1000, 4000)) {
  
  seed = c(123, 456, 789, 101112, 131415)
  clust_results = vector(mode = "list", length = 5)
  
  for (i in 1:5) {
    # In order to optimize speed of computation, draw the fifty discovery/validation samples in five batches of ten,
    # by setting five seeds.
    # Evaluating all fifty splits in a single parLapply conflicts with using ten cores,
    # as too many OpenBLAS threads are generated. Hence, the number of cores would have to be reduced. 
    # This might be different on a machine with more OpenBLAS thread capacity.
    system.time({
      cl = makeCluster(10, type = 'FORK')
      clusterSetRNGStream(cl, seed[i]) 
      clust_results1 = parLapply(cl, 1:10, function(nrsplit) {
        
        # split into discovery and validation data
        discov_ids = sample(1:9631, n)
        valid_ids = sample(setdiff(1:9631, discov_ids), n)
        df_discov = phyloseq(otu_table(ag.genus)[,discov_ids], tax_table(ag.genus))
        df_valid = phyloseq(otu_table(ag.genus)[,valid_ids], tax_table(ag.genus)) 
        
        clustResult = clustering_genera(df_discov, df_valid, nrsplit = (i-1)*10 + nrsplit, n = n, k_fix = 10)
        
        return(list(res = clustResult$res, val_result = clustResult$val_result))
      })
      
      stopCluster(cl)
      
      clust_results[[i]] = clust_results1
    })
  }
  saveRDS(clust_results, paste0("results/reduced nr methods/clustering genera/clust_results_", n, ".rds"))
  
}