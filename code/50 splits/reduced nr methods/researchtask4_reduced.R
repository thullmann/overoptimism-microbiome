########################################################################################################
### Clustering of samples #####################################################################
### 50 repetitions of optimization on the discovery data and validation on the validation data #########
########################################################################################################

library(NetCoMi)
library(phyloseq)
library(igraph)
library(cluster)
library(DirichletMultinomial)
library(mclust)
library(readxl)
library(parallel)
library(dplyr)

# two lists of association methods: 
# 1.) input for distance-based cluster algorithms (without sparsification),
# 2.) input for network-based cluster algorithms (with sparsification)
assoc_distbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "task4_distbased")
assoc_networkbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "task4_networkbased")

source("code/helpers/general_helpers.R")
source("code/helpers/cluster_algorithms.R")
source("code/optimization functions/clustering_samples.R")

ag.genus = readRDS("data/ag.genus.rds")
sample_data(ag.genus)$age = as.numeric(sample_data(ag.genus)$AGE_CORRECTED)
# ag.age = subset_samples(ag.genus, AGE_CAT %in% c("teen", "20s", "30s", "40s", "50s", "60s"))
ag.age = subset_samples(ag.genus, age >= 20 & age <= 65)

for (n in c(100, 250, 500, 1000, 3500)) {
  
  seed = c(123, 456, 789, 101112, 131415)
  clust_samples_results = vector(mode = "list", length = 5)
  
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
        discov_ids = sample(1:7145, n)
        valid_ids = sample(setdiff(1:7145, discov_ids), n)
        df_discov = phyloseq(otu_table(ag.age)[,discov_ids], tax_table(ag.age))
        df_valid = phyloseq(otu_table(ag.age)[,valid_ids], tax_table(ag.age)) 
        
        clustResult = clustering_samples(df_discov, df_valid, nrsplit = (i-1)*10 + nrsplit, n = n, k_fix = NULL)
        
        return(list(res = clustResult$res, val_result = clustResult$val_result))
      })
      
      stopCluster(cl)
      
      clust_samples_results[[i]] = clust_results1
    })
  }
  saveRDS(clust_samples_results, paste0("results/reduced nr methods/clustering samples/clust_samples_results_", n, ".rds"))
  
}
