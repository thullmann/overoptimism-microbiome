########################################################################################################
### Differential network analysis ######################################################################
### 50 repetitions of optimization on the discovery data and validation on the validation data #########
########################################################################################################

library(NetCoMi)
library(phyloseq)
library(mixedCCA)
library(igraph)
library(readxl)
library(parallel)
library(dplyr)
library(orca)

# lists of association methods
assoc_networkbased = read_excel("code/optimization functions/reduced nr methods/assoc_methods_reduced.xlsx", sheet = "tasks2_3")

source("code/helpers/general_helpers.R")
source("code/optimization functions/differential_network_analysis.R")

ag.genus = readRDS("data/ag.genus.rds")

# filter dataset according to antibiotics status
ag.antibiotics = subset_samples(ag.genus, ANTIBIOTIC_HISTORY %in% c("I have not taken antibiotics in the past year.",
                                                                    "Month", "Week"))

for (n in c(100, 250, 500)) {
  
  seed = c(123, 456, 789, 101112, 131415)
  antibiotics_results = vector(mode = "list", length = 5)
  
  for (i in 1:5) {
    # Draw the fifty discovery/validation samples in five batches of ten, by setting five seeds.
    system.time({
      cl = makeCluster(10, type = 'FORK')
      clusterSetRNGStream(cl, seed[i]) 
      antibiotics_results1 = parLapply(cl, 1:10, function(nrsplit) {
        
        # split into discovery and validation data
        discov_ids_no_antibiotics = sample(which(sample_data(ag.antibiotics)$ANTIBIOTIC_HISTORY == "I have not taken antibiotics in the past year."), n/2)
        discov_ids_antibiotics = sample(which(sample_data(ag.antibiotics)$ANTIBIOTIC_HISTORY %in% c("Month", "Week")), n/2)
        valid_ids_no_antibiotics = sample(setdiff(which(sample_data(ag.antibiotics)$ANTIBIOTIC_HISTORY == "I have not taken antibiotics in the past year."), discov_ids_no_antibiotics), n/2)
        valid_ids_antibiotics = sample(setdiff(which(sample_data(ag.antibiotics)$ANTIBIOTIC_HISTORY %in% c("Month", "Week")), discov_ids_antibiotics), n/2)
        
        df_discov_no_antibiotics = phyloseq(otu_table(ag.genus)[,discov_ids_no_antibiotics], tax_table(ag.genus))
        df_discov_antibiotics = phyloseq(otu_table(ag.genus)[,discov_ids_antibiotics], tax_table(ag.genus))
        df_valid_no_antibiotics = phyloseq(otu_table(ag.genus)[,valid_ids_no_antibiotics], tax_table(ag.genus)) 
        df_valid_antibiotics = phyloseq(otu_table(ag.genus)[,valid_ids_antibiotics], tax_table(ag.genus)) 
        
        antibiotics_result = differential_network_analysis(df_discov_no_antibiotics, df_discov_antibiotics, 
                                                           df_valid_no_antibiotics, df_valid_antibiotics, nrsplit = (i-1)*10 + nrsplit, n = n)
        
        return(list(res = antibiotics_result$res, val_result = antibiotics_result$val_result))
      })
      
      stopCluster(cl)
      
      antibiotics_results[[i]] = antibiotics_results1
    })
  }
  saveRDS(antibiotics_results, paste0("results/reduced nr methods/differential network analysis/antibiotics_results_", n, ".rds"))
  
}