########################################################################################################
### Clustering genera ##################################################################################
########################################################################################################

# read and restructure results ##

clust_results_100 = readRDS("results/reduced nr methods/clustering genera/clust_results_100.rds")
clust_results_250 = readRDS("results/reduced nr methods/clustering genera/clust_results_250.rds")
clust_results_500 = readRDS("results/reduced nr methods/clustering genera/clust_results_500.rds")
clust_results_1000 = readRDS("results/reduced nr methods/clustering genera/clust_results_1000.rds")
clust_results_4000 = readRDS("results/reduced nr methods/clustering genera/clust_results_4000.rds")

for (n in c(100, 250, 500, 1000, 4000)) {
  clust_results = vector(mode = "list", length = 50)
  for (i in 1:5) {
    for (s in 1:10) {
      clust_results[[(i-1)*10+s]] = get(paste0("clust_results_", n))[[i]][[s]] 
    }
  }
  
  val_result = clust_results[[1]]$val_result
  for (i in 2:50) {
    val_result = rbind(val_result, clust_results[[i]]$val_result)
  }
  val_result$ariDiff = val_result$ariValid - val_result$ariDiscov
  
  assign(paste0("clust_results_", n), clust_results)
  assign(paste0("val_result_", n), val_result)
  
}


# descriptive statistics for ARI difference ##

ari_diff = data.frame(n = c(100, 250, 500, 1000, 4000), 
                      mean_diff = c(mean(val_result_100$ariDiff), mean(val_result_250$ariDiff),
                                    mean(val_result_500$ariDiff), mean(val_result_1000$ariDiff),
                                    mean(val_result_4000$ariDiff)), 
                      median_diff = c(median(val_result_100$ariDiff), median(val_result_250$ariDiff),
                                      median(val_result_500$ariDiff), median(val_result_1000$ariDiff),
                                      median(val_result_4000$ariDiff)), 
                      sd_diff = c(sd(val_result_100$ariDiff), sd(val_result_250$ariDiff),
                                  sd(val_result_500$ariDiff), sd(val_result_1000$ariDiff),
                                  sd(val_result_4000$ariDiff)),
                      mean_scaled_diff = c(mean(val_result_100$ariDiff/val_result_100$ariDiscov), mean(val_result_250$ariDiff/val_result_250$ariDiscov),
                                           mean(val_result_500$ariDiff/val_result_500$ariDiscov), mean(val_result_1000$ariDiff/val_result_1000$ariDiscov),
                                           mean(val_result_4000$ariDiff/val_result_4000$ariDiscov)),
                      median_scaled_diff = c(median(val_result_100$ariDiff/val_result_100$ariDiscov), median(val_result_250$ariDiff/val_result_250$ariDiscov),
                                             median(val_result_500$ariDiff/val_result_500$ariDiscov), median(val_result_1000$ariDiff/val_result_1000$ariDiscov),
                                             median(val_result_4000$ariDiff/val_result_4000$ariDiscov)),
                      sd_scaled_diff = c(sd(val_result_100$ariDiff/val_result_100$ariDiscov), sd(val_result_250$ariDiff/val_result_250$ariDiscov),
                                         sd(val_result_500$ariDiff/val_result_500$ariDiscov), sd(val_result_1000$ariDiff/val_result_1000$ariDiscov),
                                         sd(val_result_4000$ariDiff/val_result_4000$ariDiscov)))

ari_diff$effect_diff = ari_diff$mean_diff/ari_diff$sd_diff
ari_diff$effect_scaled_diff = ari_diff$mean_scaled_diff/ari_diff$sd_scaled_diff

saveRDS(ari_diff, file = "descriptive statistics/reduced nr methods/clustering_ari_diff.rds")

# additional stability analysis

ari_stability = data.frame(n = c(100, 250, 500, 1000, 4000), 
                           mean = c(mean(val_result_100$ariStab), mean(val_result_250$ariStab),
                                    mean(val_result_500$ariStab), mean(val_result_1000$ariStab),
                                    mean(val_result_4000$ariStab)), 
                           median = c(median(val_result_100$ariStab), median(val_result_250$ariStab),
                                      median(val_result_500$ariStab), median(val_result_1000$ariStab),
                                      median(val_result_4000$ariStab)), 
                           sd = c(sd(val_result_100$ariStab), sd(val_result_250$ariStab),
                                  sd(val_result_500$ariStab), sd(val_result_1000$ariStab),
                                  sd(val_result_4000$ariStab)))

saveRDS(ari_stability, file = "descriptive statistics/reduced nr methods/clustering_stability.rds")


########################################################################################################
### hub detection #####################################################################################
########################################################################################################

# read and restructure results ##

hub_results_100 = readRDS("results/reduced nr methods/hub detection/hub_results_100.rds")
hub_results_250 = readRDS("results/reduced nr methods/hub detection/hub_results_250.rds")
hub_results_500 = readRDS("results/reduced nr methods/hub detection/hub_results_500.rds")
hub_results_1000 = readRDS("results/reduced nr methods/hub detection/hub_results_1000.rds")
hub_results_4000 = readRDS("results/reduced nr methods/hub detection/hub_results_4000.rds")

for (n in c(100, 250, 500, 1000, 4000)) {
  hub_results = vector(mode = "list", length = 10)
  for (i in 1:5) {
    for (s in 1:10) {
      hub_results[[(i-1)*10+s]] = get(paste0("hub_results_", n))[[i]][[s]] 
    }
  }
  
  val_result = hub_results[[1]]$val_result
  for (i in 2:50) {
    val_result = rbind(val_result, hub_results[[i]]$val_result)
  }
  val_result$nrHubsDiff = val_result$nrHubsValid - val_result$nrHubsDiscov
  val_result$degreeCentDiff = val_result$degreeCentValid - val_result$degreeCentDiscov 
  val_result$betweenCentDiff = val_result$betweenCentValid - val_result$betweenCentDiscov 
  val_result$closeCentDiff = val_result$closeCentValid - val_result$closeCentDiscov 
  
  assign(paste0("hub_results_", n), hub_results)
  assign(paste0("val_result_", n), val_result)
  
}

# descriptive statistics for nrHubs difference ##

nrhubs_diff = data.frame(n = c(100, 250, 500, 1000, 4000), 
                         mean_diff = c(mean(val_result_100$nrHubsDiff), mean(val_result_250$nrHubsDiff),
                                       mean(val_result_500$nrHubsDiff), mean(val_result_1000$nrHubsDiff),
                                       mean(val_result_4000$nrHubsDiff)), 
                         median_diff = c(median(val_result_100$nrHubsDiff), median(val_result_250$nrHubsDiff),
                                         median(val_result_500$nrHubsDiff), median(val_result_1000$nrHubsDiff),
                                         median(val_result_4000$nrHubsDiff)), 
                         sd_diff = c(sd(val_result_100$nrHubsDiff), sd(val_result_250$nrHubsDiff),
                                     sd(val_result_500$nrHubsDiff), sd(val_result_1000$nrHubsDiff),
                                     sd(val_result_4000$nrHubsDiff)), 
                         mean_scaled_diff = c(mean(val_result_100$nrHubsDiff/val_result_100$nrHubsDiscov), mean(val_result_250$nrHubsDiff/val_result_250$nrHubsDiscov),
                                              mean(val_result_500$nrHubsDiff/val_result_500$nrHubsDiscov), mean(val_result_1000$nrHubsDiff/val_result_1000$nrHubsDiscov),
                                              mean(val_result_4000$nrHubsDiff/val_result_4000$nrHubsDiscov)),
                         median_scaled_diff = c(median(val_result_100$nrHubsDiff/val_result_100$nrHubsDiscov), median(val_result_250$nrHubsDiff/val_result_250$nrHubsDiscov),
                                                median(val_result_500$nrHubsDiff/val_result_500$nrHubsDiscov), median(val_result_1000$nrHubsDiff/val_result_1000$nrHubsDiscov),
                                                median(val_result_4000$nrHubsDiff/val_result_4000$nrHubsDiscov)),
                         sd_scaled_diff = c(sd(val_result_100$nrHubsDiff/val_result_100$nrHubsDiscov), sd(val_result_250$nrHubsDiff/val_result_250$nrHubsDiscov),
                                            sd(val_result_500$nrHubsDiff/val_result_500$nrHubsDiscov), sd(val_result_1000$nrHubsDiff/val_result_1000$nrHubsDiscov),
                                            sd(val_result_4000$nrHubsDiff/val_result_4000$nrHubsDiscov)))

nrhubs_diff$effect_diff = nrhubs_diff$mean_diff/nrhubs_diff$sd_diff
nrhubs_diff$effect_scaled_diff = nrhubs_diff$mean_scaled_diff/nrhubs_diff$sd_scaled_diff

saveRDS(nrhubs_diff, file = "descriptive statistics/reduced nr methods/hub_detection_nr_hubs_diff.rds")

# additional stability analysis

hubs_stability = data.frame(n = c(100, 250, 500, 1000, 4000), 
                            mean_jaccard = c(mean(val_result_100$jaccard), mean(val_result_250$jaccard),
                                             mean(val_result_500$jaccard), mean(val_result_1000$jaccard),
                                             mean(val_result_4000$jaccard)), 
                            median_jaccard = c(median(val_result_100$jaccard), median(val_result_250$jaccard),
                                               median(val_result_500$jaccard), median(val_result_1000$jaccard),
                                               median(val_result_4000$jaccard)), 
                            sd_jaccard = c(sd(val_result_100$jaccard), sd(val_result_250$jaccard),
                                           sd(val_result_500$jaccard), sd(val_result_1000$jaccard),
                                           sd(val_result_4000$jaccard)),
                            mean_cosineFamilies = c(mean(val_result_100$cosineFamilies), mean(val_result_250$cosineFamilies),
                                                    mean(val_result_500$cosineFamilies), mean(val_result_1000$cosineFamilies),
                                                    mean(val_result_4000$cosineFamilies)), 
                            median_cosineFamilies = c(median(val_result_100$cosineFamilies), median(val_result_250$cosineFamilies),
                                                      median(val_result_500$cosineFamilies), median(val_result_1000$cosineFamilies),
                                                      median(val_result_4000$cosineFamilies)), 
                            sd_cosineFamilies = c(sd(val_result_100$cosineFamilies), sd(val_result_250$cosineFamilies),
                                                  sd(val_result_500$cosineFamilies), sd(val_result_1000$cosineFamilies),
                                                  sd(val_result_4000$cosineFamilies)))

saveRDS(hubs_stability, file = "descriptive statistics/reduced nr methods/hub_detection_stability.rds")

########################################################################################################
### differential network analysis ######################################################################
########################################################################################################

# read and restructure results ##

antibiotics_results_100 = readRDS("results/reduced nr methods/differential network analysis/antibiotics_results_100.rds")
antibiotics_results_250 = readRDS("results/reduced nr methods/differential network analysis/antibiotics_results_250.rds")
antibiotics_results_500 = readRDS("results/reduced nr methods/differential network analysis/antibiotics_results_500.rds")

for (n in c(100, 250, 500)) {
  antibiotics_results = vector(mode = "list", length = 10)
  for (i in 1:5) {
    for (s in 1:10) {
      antibiotics_results[[(i-1)*10+s]] = get(paste0("antibiotics_results_", n))[[i]][[s]] 
    }
  }
  
  val_result = antibiotics_results[[1]]$val_result
  for (i in 2:50) {
    val_result = rbind(val_result, antibiotics_results[[i]]$val_result)
  }
  val_result$gcdDiff = val_result$gcdValid - val_result$gcdDiscov
  assign(paste0("antibiotics_results_", n), antibiotics_results)
  assign(paste0("val_result_", n), val_result)
  
}

# descriptive statistics for GCD difference ##

gcd_diff = data.frame(n = c(100, 250, 500), 
                      mean_diff = c(mean(val_result_100$gcdDiff), mean(val_result_250$gcdDiff),
                                    mean(val_result_500$gcdDiff)), 
                      median_diff = c(median(val_result_100$gcdDiff), median(val_result_250$gcdDiff),
                                      median(val_result_500$gcdDiff)), 
                      sd_diff = c(sd(val_result_100$gcdDiff), sd(val_result_250$gcdDiff),
                                  sd(val_result_500$gcdDiff)), 
                      mean_scaled_diff = c(mean(val_result_100$gcdDiff/val_result_100$gcdDiscov), mean(val_result_250$gcdDiff/val_result_250$gcdDiscov),
                                           mean(val_result_500$gcdDiff/val_result_500$gcdDiscov)),
                      median_scaled_diff = c(median(val_result_100$gcdDiff/val_result_100$gcdDiscov), median(val_result_250$gcdDiff/val_result_250$gcdDiscov),
                                             median(val_result_500$gcdDiff/val_result_500$gcdDiscov)),
                      sd_scaled_diff = c(sd(val_result_100$gcdDiff/val_result_100$gcdDiscov), sd(val_result_250$gcdDiff/val_result_250$gcdDiscov),
                                         sd(val_result_500$gcdDiff/val_result_500$gcdDiscov)))

gcd_diff$effect_diff = gcd_diff$mean_diff/gcd_diff$sd_diff
gcd_diff$effect_scaled_diff = gcd_diff$mean_scaled_diff/gcd_diff$sd_scaled_diff

saveRDS(gcd_diff, file = "descriptive statistics/reduced nr methods/diff_network_gcd_diff.rds")

########################################################################################################
### Clustering samples ##################################################################################
########################################################################################################

# read and restructure results ##

clust_samples_results_100 = readRDS("results/reduced nr methods/clustering samples/clust_samples_results_100.rds")
clust_samples_results_250 = readRDS("results/reduced nr methods/clustering samples/clust_samples_results_250.rds")
clust_samples_results_500 = readRDS("results/reduced nr methods/clustering samples/clust_samples_results_500.rds")
clust_samples_results_1000 = readRDS("results/reduced nr methods/clustering samples/clust_samples_results_1000.rds")
clust_samples_results_3500 = readRDS("results/reduced nr methods/clustering samples/clust_samples_results_3500.rds")

for (n in c(100, 250, 500, 1000, 3500)) {
  clust_samples_results = vector(mode = "list", length = 50)
  for (i in 1:5) {
    for (s in 1:10) {
      clust_samples_results[[(i-1)*10+s]] = get(paste0("clust_samples_results_", n))[[i]][[s]] 
    }
  }
  
  val_result = clust_samples_results[[1]]$val_result
  for (i in 2:50) {
    val_result = rbind(val_result, clust_samples_results[[i]]$val_result)
  }
  val_result$swDiff = val_result$swValid - val_result$swDiscov
  
  assign(paste0("clust_samples_results_", n), clust_samples_results)
  assign(paste0("val_result_", n), val_result)
  
}

# descriptive statistics for SW difference ##

sw_diff = data.frame(n = c(100, 250, 500, 1000, 3500), 
                     mean_diff = c(mean(val_result_100$swDiff), mean(val_result_250$swDiff),
                                   mean(val_result_500$swDiff), mean(val_result_1000$swDiff),
                                   mean(val_result_3500$swDiff)), 
                     median_diff = c(median(val_result_100$swDiff), median(val_result_250$swDiff),
                                     median(val_result_500$swDiff), median(val_result_1000$swDiff),
                                     median(val_result_3500$swDiff)), 
                     sd_diff = c(sd(val_result_100$swDiff), sd(val_result_250$swDiff),
                                 sd(val_result_500$swDiff), sd(val_result_1000$swDiff),
                                 sd(val_result_3500$swDiff)), 
                     mean_scaled_diff = c(mean(val_result_100$swDiff/val_result_100$swDiscov), mean(val_result_250$swDiff/val_result_250$swDiscov),
                                          mean(val_result_500$swDiff/val_result_500$swDiscov), mean(val_result_1000$swDiff/val_result_1000$swDiscov),
                                          mean(val_result_3500$swDiff/val_result_3500$swDiscov)),
                     median_scaled_diff = c(median(val_result_100$swDiff/val_result_100$swDiscov), median(val_result_250$swDiff/val_result_250$swDiscov),
                                            median(val_result_500$swDiff/val_result_500$swDiscov), median(val_result_1000$swDiff/val_result_1000$swDiscov),
                                            median(val_result_3500$swDiff/val_result_3500$swDiscov)),
                     sd_scaled_diff = c(sd(val_result_100$swDiff/val_result_100$swDiscov), sd(val_result_250$swDiff/val_result_250$swDiscov),
                                        sd(val_result_500$swDiff/val_result_500$swDiscov), sd(val_result_1000$swDiff/val_result_1000$swDiscov),
                                        sd(val_result_3500$swDiff/val_result_3500$swDiscov)))

sw_diff$effect_diff = sw_diff$mean_diff/sw_diff$sd_diff
sw_diff$effect_scaled_diff = sw_diff$mean_scaled_diff/sw_diff$sd_scaled_diff

saveRDS(sw_diff, file = "descriptive statistics/reduced nr methods/clustering_sw_diff.rds")

