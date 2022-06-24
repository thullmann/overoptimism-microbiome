library(dplyr)
library(ggplot2)
library(gridExtra)
# devtools::install_github("eclarke/ggbeeswarm")
library(ggbeeswarm)
library(ggnewscale)

########################################################################################################
### Clustering genera ##################################################################################
########################################################################################################

# read and restructure results ##

clust_results_100 = readRDS("results/clustering genera/clust_results_100.rds")
clust_results_250 = readRDS("results/clustering genera/clust_results_250.rds")
clust_results_500 = readRDS("results/clustering genera/clust_results_500.rds")
clust_results_1000 = readRDS("results/clustering genera/clust_results_1000.rds")
clust_results_4000 = readRDS("results/clustering genera/clust_results_4000.rds")

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

saveRDS(ari_diff, file = "descriptive statistics/clustering_ari_diff.rds")

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

saveRDS(ari_stability, file = "descriptive statistics/clustering_stability.rds")


# plots ##

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
colors = gg_color_hue(11)
names(colors) = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10")

clust_results = list(clust_results_100, clust_results_250, clust_results_500, clust_results_1000, clust_results_4000)
val_results = list(val_result_100, val_result_250, val_result_500, val_result_1000, val_result_4000)
sample_size = c(100, 250, 500, 1000, 4000)
p1a = list() # results for clustering on the discovery data, not separated by either normalization or sparsification methods
p1b = list() # results for clustering on the discovery data, separated by normalization, but not sparsification methods
p1c = list() # results for network-based clustering on the discovery data, separated by normalization and sparsification methods
p2 = list() # best ARIs for the clustering of bacterial genera on the discovery data, compared with the results on validation data; methods separated by
# normalization, but not sparsification methods

for (n in 1:5) {
  
  val_results[[n]]$clustAlgo = recode(val_results[[n]]$clustAlgo, "cluster_fast_greedy" = "fg.modular", "cluster_louvain" = "louvain", 
                                      "hierarchical" = "hierarch")
  val_results[[n]]$method = paste(val_results[[n]]$clustAlgo, val_results[[n]]$measure, sep = ", ")
  val_results[[n]]$methodNorm = paste(val_results[[n]]$clustAlgo, val_results[[n]]$measure, val_results[[n]]$normMethod, sep = ", ")
  val_results[[n]]$methodNorm = gsub(pattern = "spring, none", replacement = "spring, mclr", x = val_results[[n]]$methodNorm)
  val_results[[n]]$methodSpars = paste(val_results[[n]]$clustAlgo, val_results[[n]]$measure, val_results[[n]]$normMethod, val_results[[n]]$sparsMethod, sep = ", ")
  val_results[[n]]$ari = val_results[[n]]$ariDiscov
  
  clust_results[[n]][[1]]$res$split = 1
  clust_results[[n]][[1]]$res$bestAri = val_results[[n]][1,"ariDiscov"]
  clust_results[[n]][[1]]$res$bestMethod = val_results[[n]][1,"method"]
  clust_results[[n]][[1]]$res$bestMethodNorm = val_results[[n]][1,"methodNorm"]
  clust_results[[n]][[1]]$res$bestMethodSpars = val_results[[n]][1,"methodSpars"]
  res = clust_results[[n]][[1]]$res
  for (i in 2:50) {
    clust_results[[n]][[i]]$res$split = i
    clust_results[[n]][[i]]$res$bestAri = val_results[[n]][i,"ariDiscov"]
    clust_results[[n]][[i]]$res$bestMethod = val_results[[n]][i,"method"]
    clust_results[[n]][[i]]$res$bestMethodNorm = val_results[[n]][i,"methodNorm"]
    clust_results[[n]][[i]]$res$bestMethodSpars = val_results[[n]][i,"methodSpars"]
    res = rbind(res, clust_results[[n]][[i]]$res)
  }
  
  res$clustAlgo = recode(res$clustAlgo, "cluster_fast_greedy" = "fg.modular", "cluster_louvain" = "louvain", 
                         "hierarchical" = "hierarch")
  res$method = paste(res$clustAlgo, res$measure, sep = ", ")
  res$method = factor(res$method, ordered = TRUE, levels = c("hierarch, pearson", "hierarch, spearman", "hierarch, latentcor", "hierarch, propr",
                                                             "spectral, pearson", "spectral, spearman", "spectral, latentcor", "spectral, propr",
                                                             "fg.modular, pearson", "fg.modular, spearman", "fg.modular, spring", "fg.modular, propr",
                                                             "louvain, pearson", "louvain, spearman", "louvain, spring", "louvain, propr",
                                                             "manta, pearson", "manta, spearman", "manta, spring", "manta, propr"))
  
  res$k_mod = ifelse(res$k <= 10, res$k, ">10")
  res$k_mod = factor(res$k_mod, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
  
  p1a[[n]] = ggplot(data = res, aes(x = method, y = ari)) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.6, outlier.shape = 4) + 
    geom_beeswarm(aes(colour = k_mod), cex = 0.3, size = 0.05, corral = "random", corral.width = 0.9) +
    scale_color_manual(values = colors, name = "k") +
    guides(color = guide_legend(override.aes = list(size = 0.8))) +
    geom_vline(xintercept = c(4.5, 8.5, 12.5, 16.5)) + 
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("ARI") +
    labs(title = paste0("n = ", sample_size[n])) +
    annotate("text", label = "hierarchical", x = 1.8, y = 0.26, size = 2) +
    annotate("text", label = "spectral", x = 5.6, y = 0.26, size = 2) +
    annotate("text", label = "fast greedy modularity", x = 10.4, y = 0.26, size = 2) +
    annotate("text", label = "Louvain", x = 13.6, y = 0.26, size = 2) +
    annotate("text", label = "manta", x = 17.4, y = 0.26, size = 2) +
    ylim(-0.1, 0.26) + 
    geom_beeswarm(aes(x = bestMethod, y = bestAri), cex = 0.3, shape = 22, color = "red", size = 0.8, stat = "unique") 

  
  levels_norm = c("hierarch, pearson, clr", "hierarch, pearson, VST", "hierarch, pearson, mclr",
                  "hierarch, spearman, clr", "hierarch, spearman, VST", "hierarch, spearman, mclr",
                  "hierarch, latentcor, mclr", 
                  "hierarch, propr, none",
                  "spectral, pearson, clr", "spectral, pearson, VST", "spectral, pearson, mclr", 
                  "spectral, spearman, clr", "spectral, spearman, VST","spectral, spearman, mclr",
                  "spectral, latentcor, mclr", 
                  "spectral, propr, none",
                  "fg.modular, pearson, clr", "fg.modular, pearson, VST", "fg.modular, pearson, mclr",
                  "fg.modular, spearman, clr", "fg.modular, spearman, VST", "fg.modular, spearman, mclr", 
                  "fg.modular, spring, mclr", 
                  "fg.modular, propr, none",
                  "louvain, pearson, clr", "louvain, pearson, VST", "louvain, pearson, mclr", 
                  "louvain, spearman, clr", "louvain, spearman, VST", "louvain, spearman, mclr",  
                  "louvain, spring, mclr", 
                  "louvain, propr, none",
                  "manta, pearson, clr", "manta, pearson, VST", "manta, pearson, mclr", 
                  "manta, spearman, clr", "manta, spearman, VST", "manta, spearman, mclr", 
                  "manta, spring, mclr", 
                  "manta, propr, none")
  
  res$methodNorm = paste(res$clustAlgo, res$measure, res$normMethod, sep = ", ")
  res$methodNorm = gsub(pattern = "spring, none", replacement = "spring, mclr", x = res$methodNorm)
  res$methodNorm = factor(res$methodNorm, ordered = TRUE, levels = levels_norm)
  
  res$bestMethodNorm = gsub(pattern = "spring, none", replacement = "spring, mclr", x = res$bestMethodNorm)
  res$bestMethodNorm = factor(res$bestMethodNorm, ordered = TRUE, levels = levels_norm)
  
  
  p1b[[n]] = ggplot(data = res, aes(x = methodNorm, y = ari)) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.6, outlier.shape = 4) + 
    geom_beeswarm(aes(colour = k_mod), cex = 0.2, size = 0.1, corral = "random", corral.width = 0.9) +
    scale_color_manual(values= colors, name = "k") +
    guides(color = guide_legend(override.aes = list(size = 0.8))) +
    geom_vline(xintercept = c(8.5, 16.5, 24.5, 32.5)) + 
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("ARI") +
    labs(title = paste0("n = ", sample_size[n])) +
    annotate("text", label = "hierarchical", x = 2.3, y = 0.26, size = 2) +
    annotate("text", label = "spectral", x = 10.1, y = 0.26, size = 2) +
    annotate("text", label = "fast greedy modularity", x = 20.1, y = 0.26, size = 2) +
    annotate("text", label = "Louvain", x = 26.1, y = 0.26, size = 2) +
    annotate("text", label = "manta", x = 33.9, y = 0.26, size = 2) +
    ylim(-0.1, 0.26) + 
    geom_beeswarm(aes(x = bestMethodNorm, y = bestAri), cex = 0.3, shape = 22, color = "red", size = 0.8, stat = "unique") 

  val_results[[n]]$k_mod_discov = ifelse(val_results[[n]]$kDiscov <= 10, val_results[[n]]$kDiscov, ">10")
  val_results[[n]]$k_mod_discov = factor(val_results[[n]]$k_mod_discov, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
  val_results[[n]]$k_mod_valid = ifelse(val_results[[n]]$kValid <= 10, val_results[[n]]$kValid, ">10")
  val_results[[n]]$k_mod_valid = factor(val_results[[n]]$k_mod_valid, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))
  
  val_results[[n]]$split = 1:50
  val_result1 = val_results[[n]]
  val_result1$DiscovValid = "discov"
  val_result2 = val_results[[n]]
  val_result2$DiscovValid = "valid"
  
  val_result1$ari = val_result1$ariDiscov
  val_result2$ari = val_result2$ariValid
  
  val_result1$k_mod = val_result1$k_mod_discov
  val_result2$k_mod = val_result2$k_mod_valid
  
  val_result = rbind(val_result1, val_result2)
  
  val_result$methodNorm = factor(val_result$methodNorm, ordered = TRUE, levels = levels_norm)
  
  
  p2[[n]] = ggplot(val_result) +
    geom_point(aes(x = DiscovValid, y = ari, group = split, color = k_mod)) +
    scale_color_manual(values = colors, name = "k") +
    guides(color = guide_legend(override.aes = list(size = 0.8))) +
    geom_line(aes(x = DiscovValid, y = ari, group = split), alpha = 0.2)+
    facet_grid(.~methodNorm, scales="fixed", switch = "x", labeller = label_wrap_gen(width=10)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + labs(x = "", y = "ARI") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(0.02, 0.25) +
    theme(panel.spacing.x = unit(0,"line")) +
    theme(strip.placement = 'outside',
          strip.background.x = element_blank())
  
  
  res = res[res$clustAlgo %in% c("louvain", "fg.modular", "manta"),]
  res$isSparse = "sparse"
  res$isSparse[!grepl(pattern = 'fg.modular|louvain|manta', x = res$bestMethodSpars)] = "not sparse"
  res$isSparse = factor(res$isSparse, ordered = TRUE, levels = c("sparse", "not sparse"))
  
  levels_sparse = c("fg.modular, pearson, clr, t-test", "fg.modular, pearson, clr, threshold",
                    "fg.modular, pearson, VST, t-test",  "fg.modular, pearson, VST, threshold", 
                    "fg.modular, pearson, mclr, t-test", "fg.modular, pearson, mclr, threshold",
                    "fg.modular, spearman, clr, t-test", "fg.modular, spearman, clr, threshold", 
                    "fg.modular, spearman, VST, t-test", "fg.modular, spearman, VST, threshold", 
                    "fg.modular, spearman, mclr, t-test", "fg.modular, spearman, mclr, threshold", 
                    "fg.modular, spring, mclr, NB selection",
                    "fg.modular, propr, none, threshold",
                    "louvain, pearson, clr, t-test", "louvain, pearson, clr, threshold", 
                    "louvain, pearson, VST, t-test", "louvain, pearson, VST, threshold", 
                    "louvain, pearson, mclr, t-test", "louvain, pearson, mclr, threshold", 
                    "louvain, spearman, clr, t-test", "louvain, spearman, clr, threshold", 
                    "louvain, spearman, VST, t-test", "louvain, spearman, VST, threshold", 
                    "louvain, spearman, mclr, t-test", "louvain, spearman, mclr, threshold", 
                    "louvain, spring, mclr, NB selection", 
                    "louvain, propr, none, threshold",
                    "manta, pearson, clr, t-test", "manta, pearson, clr, threshold", 
                    "manta, pearson, VST, t-test", "manta, pearson, VST, threshold", 
                    "manta, pearson, mclr, t-test", "manta, pearson, mclr, threshold", 
                    "manta, spearman, clr, t-test", "manta, spearman, clr, threshold", 
                    "manta, spearman, VST, t-test", "manta, spearman, VST, threshold", 
                    "manta, spearman, mclr, t-test", "manta, spearman, mclr, threshold", 
                    "manta, spring, mclr, NB selection", 
                    "manta, propr, none, threshold")
  
  res$methodSpars = paste(res$clustAlgo, res$measure, res$normMethod, res$sparsMethod, sep = ", ")
  res$methodSpars = gsub(pattern = "spring, none, none", replacement = "spring, mclr, NB selection", x = res$methodSpars)
  res$methodSpars = factor(res$methodSpars, ordered = TRUE, levels = levels_sparse)
  
  res$bestMethodSpars = gsub(pattern = "spring, none, none", replacement = "spring, mclr, NB selection", x = res$bestMethodSpars)
  # for the splits in which hierarchical or spectral clustering was chosen as the best method,
  # replace bestMethodSpars with an auxiliary "fake value" (here: "manta, pearson, clr, t-test")
  # These fake values will not be seen in the final plot,
  # and merely serve to keep the correct order of the method levels.
  res$bestMethodSpars[!grepl(pattern = 'fg.modular|louvain|manta', x = res$bestMethodSpars)] = "manta, pearson, clr, t-test"
  res$bestMethodSpars = factor(res$bestMethodSpars, ordered = TRUE, levels = levels_sparse)
  
  p1c[[n]] = ggplot(data = res, aes(x = methodSpars, y = ari)) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.6, outlier.shape = 4) + 
    geom_beeswarm(aes(colour = k_mod), cex = 0.2, size = 0.1, corral = "random", corral.width = 0.9) +
    scale_color_manual(values= colors, name = "k") +
    geom_vline(xintercept = c(14.5, 28.5)) + 
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("ARI") +
    labs(title = paste0("n = ", sample_size[n])) +
    annotate("text", label = "fast greedy modularity", x = 4.2, y = 0.26, size = 2) +
    annotate("text", label = "Louvain", x = 16.1, y = 0.26, size = 2) +
    annotate("text", label = "manta", x = 30.0, y = 0.26, size = 2) +
    ylim(-0.1, 0.26) + 
    new_scale_color() +
    geom_beeswarm(aes(x = bestMethodSpars, y = bestAri, colour = isSparse), cex = 0.3, shape = 22, size = 0.8, stat = "unique") +
    scale_color_manual(values = c("red", NA), guide = FALSE) # the auxiliary "fake" values are made transparent 
    
  ggsave(paste0("plots/clustering genera/p1a_", sample_size[n],".png"), plot = p1a[[n]], dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(paste0("plots/clustering genera/p1b_", sample_size[n],".png"), plot = p1b[[n]], dpi = 600, width = 18, height = 12, units = "cm")
  ggsave(paste0("plots/clustering genera/p1c_", sample_size[n],".png"), plot = p1c[[n]], dpi = 600, width = 18, height = 12, units = "cm")
  if (n == 1) {
    ggsave(paste0("plots/clustering genera/p2_", sample_size[n],".png"), plot = p2[[n]], dpi = 600, width = 25, height = 9, units = "cm")
  } else {
    ggsave(paste0("plots/clustering genera/p2_", sample_size[n],".png"), plot = p2[[n]], dpi = 600, width = 17, height = 9, units = "cm")
  }
  
}

########################################################################################################
### hub detection #####################################################################################
########################################################################################################

# read and restructure results ##

hub_results_100 = readRDS("results/hub detection/hub_results_100.rds")
hub_results_250 = readRDS("results/hub detection/hub_results_250.rds")
hub_results_500 = readRDS("results/hub detection/hub_results_500.rds")
hub_results_1000 = readRDS("results/hub detection/hub_results_1000.rds")
hub_results_4000 = readRDS("results/hub detection/hub_results_4000.rds")

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

saveRDS(nrhubs_diff, file = "descriptive statistics/hub_detection_nr_hubs_diff.rds")

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

saveRDS(hubs_stability, file = "descriptive statistics/hub_detection_stability.rds")


# plots ##

hub_results = list(hub_results_100, hub_results_250, hub_results_500, hub_results_1000, hub_results_4000)
val_results = list(val_result_100, val_result_250, val_result_500, val_result_1000, val_result_4000)
sample_size = c(100, 250, 500, 1000, 4000)
p1a = list() #results for hub detection on the discovery data, separated by normalization, but not sparsification methods
p1b = list() # results for hub detection on the discovery data, separated by both normalization and sparsification methods
p2 = list() # highest numbers of hubs on the discovery data, compared with the results on validation data; methods separated by 
# both normalization and sparsification methods

for (n in 1:5) {
  
  val_results[[n]]$method = paste(val_results[[n]]$measure, val_results[[n]]$normMethod, sep = ", ")
  val_results[[n]]$methodSpars = paste(val_results[[n]]$measure, val_results[[n]]$normMethod, val_results[[n]]$sparsMethod, sep = ", ")
  val_results[[n]]$nrHubs = val_results[[n]]$nrHubsDiscov
  val_results[[n]]$method[val_results[[n]]$method == "spring, none"] = "spring, mclr"
  val_results[[n]]$methodSpars[val_results[[n]]$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  
  hub_results[[n]][[1]]$res$split = 1
  hub_results[[n]][[1]]$res$highestNrHubs = val_results[[n]][1,"nrHubsDiscov"]
  hub_results[[n]][[1]]$res$bestMethod = val_results[[n]][1,"method"]
  hub_results[[n]][[1]]$res$bestMethodSpars = val_results[[n]][1,"methodSpars"]
  res = hub_results[[n]][[1]]$res
  for (i in 2:50) {
    hub_results[[n]][[i]]$res$split = i
    hub_results[[n]][[i]]$res$highestNrHubs = val_results[[n]][i,"nrHubsDiscov"]
    hub_results[[n]][[i]]$res$bestMethod = val_results[[n]][i,"method"]
    hub_results[[n]][[i]]$res$bestMethodSpars = val_results[[n]][i,"methodSpars"]
    res = rbind(res, hub_results[[n]][[i]]$res)
  }
  
  res$method = paste(res$measure, res$normMethod, sep = ", ")
  res$method[res$method == "spring, none"] = "spring, mclr"
  res$method = factor(res$method, ordered = TRUE, levels = c("pearson, clr", "pearson, VST", "pearson, mclr", 
                                                             "spearman, clr", "spearman, VST", "spearman, mclr",
                                                             "spring, mclr", "propr, none"))
  res$methodSpars = paste(res$measure, res$normMethod, res$sparsMethod, sep = ", ")
  res$methodSpars[res$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  res$methodSpars = factor(res$methodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                       "pearson, VST, t-test", "pearson, VST, threshold",
                                                                       "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                       "spearman, clr, t-test", "spearman, clr, threshold",
                                                                       "spearman, VST, t-test", "spearman, VST, threshold",
                                                                       "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                       "spring, mclr, NB selection", "propr, none, threshold"))
  res$bestMethodSpars[res$bestMethodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  res$bestMethodSpars = factor(res$bestMethodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                               "pearson, VST, t-test", "pearson, VST, threshold",
                                                                               "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                               "spearman, clr, t-test", "spearman, clr, threshold",
                                                                               "spearman, VST, t-test", "spearman, VST, threshold",
                                                                               "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                               "spring, mclr, NB selection", "propr, none, threshold"))
  
  
  p1a[[n]] = ggplot(data = res, aes(x = method, y = nrHubs)) +
    # geom_beeswarm(cex = 0.3, size = 0.2) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.5, outlier.shape = 4) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("nr of hubs") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(-0.1, 15) + 
    geom_beeswarm(data = res[seq(1, 700, by = 14),], aes(x = bestMethod, y = highestNrHubs), groupOnX = TRUE, shape = 15, color = "red", size = 0.6) 
  #geom_text(data = val_results[[n]], aes(label=split), stat = "unique", hjust=-0.5, vjust=0.5, size = 1, check_overlap = TRUE)
  
  p1b[[n]] = ggplot(data = res, aes(x = methodSpars, y = nrHubs)) +
    # geom_beeswarm(groupOnX = TRUE, cex = 0.3, size = 0.2) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.5, outlier.shape = 4) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("nr of hubs") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(-0.1, 15) + 
    geom_beeswarm(data = res[seq(1, 700, by = 14),], aes(x = bestMethodSpars, y = highestNrHubs), groupOnX = TRUE, shape = 15, color = "red", size = 0.6) 
  #geom_text(data = val_results[[n]], aes(label=split), stat = "unique", hjust=-0.5, vjust=0.5, size = 1, check_overlap = TRUE)
  
  
  val_results[[n]]$split = 1:50
  val_result1 = val_results[[n]]
  val_result1$DiscovValid = "discov"
  val_result2 = val_results[[n]]
  val_result2$DiscovValid = "valid"
  
  val_result1$nrHubs = val_result1$nrHubsDiscov 
  val_result2$nrHubs = val_result2$nrHubsValid    
  
  val_result = rbind(val_result1, val_result2)
  
  val_result$method = factor(val_result$method, ordered = TRUE, levels = levels(res$method))
  val_result$methodSpars[val_result$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  val_result$methodSpars = factor(val_result$methodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                                     "pearson, VST, t-test", "pearson, VST, threshold",
                                                                                     "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                                     "spearman, clr, t-test", "spearman, clr, threshold",
                                                                                     "spearman, VST, t-test", "spearman, VST, threshold",
                                                                                     "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                                     "spring, mclr, NB selection", "propr, none, threshold"))
  levels(val_result$methodSpars) = c("pears., clr, t-test", "pears., clr, thresh.",
                                     "pears., VST, t-test", "pears., VST, thresh.",
                                     "pears., mclr, t-test", "pears., mclr, thresh.",
                                     "spear., clr, t-test", "spear., clr, thresh.",
                                     "spear., VST, t-test", "spear., VST, thresh.",
                                     "spear., mclr, t-test", "spear., mclr, thresh.",
                                     "spring, mclr, NB select.", "propr, none, thresh.")
  
  
  p2[[n]] = ggplot(val_result) +
    geom_point(aes(x = DiscovValid, y = nrHubs, group = split, color = DiscovValid, shape = DiscovValid), size = 0.8) +
    scale_color_manual(values=c("red", "black")) +
    scale_shape_manual(values = c(15, 20)) + 
    theme(legend.position = "none") +
    geom_line(aes(x = DiscovValid, y = nrHubs, group = split), alpha = 0.2)+
    facet_grid(.~methodSpars, scales="fixed", switch = "x", labeller = label_wrap_gen(width=10)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + #, strip.text.x = element_text(angle = 70, hjust = -0.5)) + 
    labs(x = "", y = "nr of hubs") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(0.00, 15) +
    #theme_bw() +
    #theme(panel.border = element_rect(fill = NA, color = "black")) +
    theme(panel.spacing.x = unit(0,"line")) +
    theme(strip.placement = 'outside',
          strip.background.x = element_blank())
  
  ggsave(paste0("plots/hub detection/p1a_", sample_size[n],".png"), plot = p1a[[n]], dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(paste0("plots/hub detection/p1b_", sample_size[n],".png"), plot = p1b[[n]], dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(paste0("plots/hub detection/p2_", sample_size[n],".png"), plot = p2[[n]], dpi = 600, width = 15, height = 9, units = "cm")
  
  
}

########################################################################################################
### differential network analysis ######################################################################
########################################################################################################

# read and restructure results ##

antibiotics_results_100 = readRDS("results/differential network analysis/antibiotics_results_100.rds")
antibiotics_results_250 = readRDS("results/differential network analysis/antibiotics_results_250.rds")
antibiotics_results_500 = readRDS("results/differential network analysis/antibiotics_results_500.rds")

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

saveRDS(gcd_diff, file = "descriptive statistics/diff_network_gcd_diff.rds")

# plots ##

antibiotics_results = list(antibiotics_results_100, antibiotics_results_250, antibiotics_results_500)
val_results = list(val_result_100, val_result_250, val_result_500)
sample_size = c(100, 250, 500)
p1a = list() #results for GCD on the discovery data, separated by normalization, but not sparsification methods
p1b = list() # results for GCD on the discovery data, separated by both normalization and sparsification methods
p2 = list() # largest GCDs on the discovery data, compared with the results on validation data; methods separated by 
# both normalization and sparsification methods

for (n in 1:3) {
  
  val_results[[n]]$method = paste(val_results[[n]]$measure, val_results[[n]]$normMethod, sep = ", ")
  val_results[[n]]$methodSpars = paste(val_results[[n]]$measure, val_results[[n]]$normMethod, val_results[[n]]$sparsMethod, sep = ", ")
  val_results[[n]]$gcd = val_results[[n]]$gcdDiscov
  val_results[[n]]$method[val_results[[n]]$method == "spring, none"] = "spring, mclr"
  val_results[[n]]$methodSpars[val_results[[n]]$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  
  antibiotics_results[[n]][[1]]$res$split = 1
  antibiotics_results[[n]][[1]]$res$largestGCD = val_results[[n]][1,"gcdDiscov"]
  antibiotics_results[[n]][[1]]$res$bestMethod = val_results[[n]][1,"method"]
  antibiotics_results[[n]][[1]]$res$bestMethodSpars = val_results[[n]][1,"methodSpars"]
  res = antibiotics_results[[n]][[1]]$res
  for (i in 2:50) {
    antibiotics_results[[n]][[i]]$res$split = i
    antibiotics_results[[n]][[i]]$res$largestGCD = val_results[[n]][i,"gcdDiscov"]
    antibiotics_results[[n]][[i]]$res$bestMethod = val_results[[n]][i,"method"]
    antibiotics_results[[n]][[i]]$res$bestMethodSpars = val_results[[n]][i,"methodSpars"]
    res = rbind(res, antibiotics_results[[n]][[i]]$res)
  }
  
  res$method = paste(res$measure, res$normMethod, sep = ", ")
  res$method[res$method == "spring, none"] = "spring, mclr"
  res$method = factor(res$method, ordered = TRUE, levels = c("pearson, clr", "pearson, VST", "pearson, mclr", 
                                                             "spearman, clr", "spearman, VST", "spearman, mclr",
                                                             "spring, mclr", "propr, none"))
  res$methodSpars = paste(res$measure, res$normMethod, res$sparsMethod, sep = ", ")
  res$methodSpars[res$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  res$methodSpars = factor(res$methodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                       "pearson, VST, t-test", "pearson, VST, threshold",
                                                                       "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                       "spearman, clr, t-test", "spearman, clr, threshold",
                                                                       "spearman, VST, t-test", "spearman, VST, threshold",
                                                                       "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                       "spring, mclr, NB selection", "propr, none, threshold"))
  res$bestMethodSpars[res$bestMethodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  res$bestMethodSpars = factor(res$bestMethodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                               "pearson, VST, t-test", "pearson, VST, threshold",
                                                                               "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                               "spearman, clr, t-test", "spearman, clr, threshold",
                                                                               "spearman, VST, t-test", "spearman, VST, threshold",
                                                                               "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                               "spring, mclr, NB selection", "propr, none, threshold"))
  
  
  p1a[[n]] = ggplot(data = res, aes(x = method, y = gcd)) +
    #geom_point(size = 0.2) +
    geom_boxplot(lwd = 0.3, outlier.size = 0.5, outlier.shape = 4) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("GCD") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(-0.1, 3.5) +
    geom_point(data = res[seq(1, 700, by = 14),], aes(x = bestMethod, y = largestGCD), shape = 15, size = 0.5, color = "red", stat = "unique") 
  
  p1b[[n]] = ggplot(data = res, aes(x = methodSpars, y = gcd)) +
    #geom_point(size = 0.2) +    
    geom_boxplot(lwd = 0.3, outlier.size = 0.5, outlier.shape = 4) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + ylab("GCD") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(-0.1, 3.5) + 
    geom_point(data = res[seq(1, 700, by = 14),], aes(x = bestMethodSpars, y = largestGCD), shape = 15, size = 0.5, color = "red", stat = "unique") 
  #geom_text(data = val_results[[n]], aes(label=split), stat = "unique", hjust=-0.5, vjust=0.5, size = 1, check_overlap = TRUE)
  
  
  val_results[[n]]$split = 1:50
  val_result1 = val_results[[n]]
  val_result1$DiscovValid = "discov"
  val_result2 = val_results[[n]]
  val_result2$DiscovValid = "valid"
  
  val_result1$gcd = val_result1$gcdDiscov 
  val_result2$gcd = val_result2$gcdValid
  
  val_result = rbind(val_result1, val_result2)
  
  val_result$method = factor(val_result$method, ordered = TRUE, levels = levels(res$method))
  val_result$methodSpars[val_result$methodSpars == "spring, none, none"] = "spring, mclr, NB selection"
  val_result$methodSpars = factor(val_result$methodSpars, ordered = TRUE, levels = c("pearson, clr, t-test", "pearson, clr, threshold",
                                                                                     "pearson, VST, t-test", "pearson, VST, threshold",
                                                                                     "pearson, mclr, t-test", "pearson, mclr, threshold",
                                                                                     "spearman, clr, t-test", "spearman, clr, threshold",
                                                                                     "spearman, VST, t-test", "spearman, VST, threshold",
                                                                                     "spearman, mclr, t-test", "spearman, mclr, threshold",
                                                                                     "spring, mclr, NB selection", "propr, none, threshold"))
  levels(val_result$methodSpars) = c("pears., clr, t-test", "pears., clr, thresh.",
                                     "pears., VST, t-test", "pears., VST, thresh.",
                                     "pears., mclr, t-test", "pears., mclr, thresh.",
                                     "spear., clr, t-test", "spear., clr, thresh.",
                                     "spear., VST, t-test", "spear., VST, thresh.",
                                     "spear., mclr, t-test", "spear., mclr, thresh.",
                                     "spring, mclr, NB select.", "propr, none, thresh.")
  
  
  p2[[n]] = ggplot(val_result) +
    geom_point(aes(x = DiscovValid, y = gcd, group = split, color = DiscovValid, shape = DiscovValid), size = 0.6) +
    scale_color_manual(values=c("red", "black")) +
    scale_shape_manual(values = c(15, 20)) + 
    theme(legend.position = "none") +
    geom_line(aes(x = DiscovValid, y = gcd, group = split), alpha = 0.2)+
    facet_grid(.~methodSpars, scales="fixed", switch = "x", labeller = label_wrap_gen(width=8)) +
    theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 70, hjust = 1)) + labs(x = "", y = "GCD") +
    labs(title = paste0("n = ", sample_size[n])) +
    ylim(-0.1, 3.5) +
    #theme_bw() +
    #theme(panel.border = element_rect(fill = NA, color = "black")) +
    theme(panel.spacing.x = unit(0,"line")) +
    theme(strip.placement = 'outside',
          strip.background.x = element_blank())
  
  ggsave(paste0("plots/differential network analysis/p1a_", sample_size[n],".png"), plot = p1a[[n]], dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(paste0("plots/differential network analysis/p1b_", sample_size[n],".png"), plot = p1b[[n]], dpi = 600, width = 15, height = 12, units = "cm")
  ggsave(paste0("plots/differential network analysis/p2_", sample_size[n],".png"), plot = p2[[n]], dpi = 600, width = 15, height = 9, units = "cm")
  
}
