########################################################################################################
### Hub detection ######################################################################################
### 50 repetitions of optimization on the discovery data and validation on the validation data #########
########################################################################################################

library(NetCoMi)
library(phyloseq)
library(mixedCCA)
library(igraph)
library(readxl)
library(parallel)
library(dplyr)

# lists of association methods
assoc_networkbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "tasks2_3")

source("code/helpers/general_helpers.R")
source("code/optimization functions/hub_detection.R")

ag.genus = readRDS("data/ag.genus.rds")

for (n in c(100, 250, 500, 1000, 4000)) {
  
  seed = c(123, 456, 789, 101112, 131415)
  hub_results = vector(mode = "list", length = 5)
  
  for (i in 1:5) {
    # Draw the fifty discovery/validation samples in five batches of ten, by setting five seeds.
    system.time({
      cl = makeCluster(10, type = 'FORK')
      clusterSetRNGStream(cl, seed[i]) 
      hub_results1 = parLapply(cl, 1:10, function(nrsplit) {
        
        # split into discovery and validation data
        discov_ids = sample(1:9631, n)
        valid_ids = sample(setdiff(1:9631, discov_ids), n)
        df_discov = phyloseq(otu_table(ag.genus)[,discov_ids], tax_table(ag.genus))
        df_valid = phyloseq(otu_table(ag.genus)[,valid_ids], tax_table(ag.genus))
        
        hubResult = hub_detection(df_discov, df_valid, nrsplit = (i-1)*10 + nrsplit, n = n)
        
        return(list(res = hubResult$res, val_result = hubResult$val_result))
      })
      
      stopCluster(cl)
      
      hub_results[[i]] = hub_results1
    })
  }
  saveRDS(hub_results, paste0("results/reduced nr methods/hub detection/hub_results_", n, ".rds"))
}