hub_detection = function(df_discov, df_valid, nrsplit, n) {
  
  ## applying different method combinations for network construction ###########################################
  
  res = data.frame(measure = character(), zeroMethod = character(), normMethod = character(), 
                   sparsMethod = character(),
                   nrHubs = numeric(), degreeCent = numeric(), betweenCent = numeric(), closeCent = numeric(),
                   hubs = I(list()))
  
  for (i in 1:nrow(assoc_networkbased)) {

    if (assoc_networkbased$normMethod[i] == "VST") {
      normPar = list(fitType = "local")
    } else {
      normPar = NULL
    }
    
    nc = netConstruct(df_discov, measure = assoc_networkbased$measure[i],
                      zeroMethod = assoc_networkbased$zeroMethod[i],
                      normMethod = assoc_networkbased$normMethod[i], normPar = normPar,
                      sparsMethod = assoc_networkbased$sparsMethod[i], 
                      thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                      dissFunc = "signed",
                      seed = 123456)
    
    # determine nodes with centrality values above the 0.95 quantile for 
    # degree, betweenness and closeness centrality simultaneously 
    analysis = netAnalyze(net = nc, connectivity  = FALSE, hubPar = c("degree", "betweenness", "closeness"),
               hubQuant = 0.95, lnormFit = FALSE)
    nrHubs = length(analysis$hubs$hubs1)
    degreeCent = analysis$centralities$degree1[names(analysis$centralities$degree1) %in% analysis$hubs$hubs1]
    betweenCent = analysis$centralities$between1[names(analysis$centralities$between1) %in% analysis$hubs$hubs1]
    closeCent = analysis$centralities$close1[names(analysis$centralities$close1) %in% analysis$hubs$hubs1]

    res = rbind(res, data.frame(measure = assoc_networkbased$measure[i], zeroMethod = assoc_networkbased$zeroMethod[i], 
                                normMethod = assoc_networkbased$normMethod[i], 
                                sparsMethod = assoc_networkbased$sparsMethod[i],
                                nrHubs = nrHubs, degreeCent = mean(degreeCent), betweenCent = mean(betweenCent),
                                closeCent = mean(closeCent), hubs = I(list(analysis$hubs$hubs1))))
    
    
  }
  
  
  ## optimization ####################################################################################################
  
  best_results = which(res$nrHubs == max(res$nrHubs))
  if(length(best_results) > 1) {
    rankMatrix = data.frame(rank1 = rank(res[best_results, "degreeCent"]), rank2 = rank(res[best_results, "betweenCent"]),
                            rank3 = rank(res[best_results, "closeCent"]))
    meanRanks = apply(rankMatrix, MARGIN = 1, FUN = mean)
    # return method with best mean rank of the mean centralities 
    best_rank = which(meanRanks == max(meanRanks))[1]
    best_result = best_results[best_rank]
  } else {
    best_result = best_results
  }
  measure = res[best_result,"measure"]
  zeroMethod = res[best_result, "zeroMethod"]
  normMethod = res[best_result, "normMethod"]
  sparsMethod = res[best_result, "sparsMethod"]
  nrHubsDiscov = res[best_result, "nrHubs"]
  degreeCentDiscov = res[best_result, "degreeCent"]
  betweenCentDiscov = res[best_result, "betweenCent"]
  closeCentDiscov = res[best_result, "closeCent"]
  hubsDiscov = res$hubs[[best_result]]
  
  
  ## validation #####################################################################################################
  
  if (normMethod == "VST") {
    normPar = list(fitType = "local")
  } else {
    normPar = NULL
  }
  
  nc = netConstruct(df_valid, measure = measure,
                      zeroMethod = zeroMethod,
                      normMethod = normMethod, normPar = normPar,
                      sparsMethod = sparsMethod, 
                      thresh = 0.15, alpha = 0.05, adjust = "lfdr",
                      dissFunc = "signed", weighted = TRUE,
                      seed = 123456)
  analysis = netAnalyze(net = nc, connectivity = FALSE, hubPar = c("degree", "betweenness", "closeness"),
                        hubQuant = 0.95, lnormFit = FALSE)
  
  hubsValid = analysis$hubs$hubs1
  nrHubsValid = length(hubsValid)
  jaccard = as.numeric(calc_jaccard(hubsDiscov, hubsValid))
  
  allFamilies = unique(ag.genus@tax_table@.Data[,5][c(hubsDiscov, hubsValid)])
  tableDiscov = as.numeric(table(factor(ag.genus@tax_table@.Data[,5][hubsDiscov], levels = allFamilies)))
  tableValid = as.numeric(table(factor(ag.genus@tax_table@.Data[,5][hubsValid], levels = allFamilies)))
  cosineFamilies = sum(tableValid*tableDiscov)/sqrt(sum(tableValid^2)*sum(tableDiscov^2))

  degreeCentValid = analysis$centralities$degree1[names(analysis$centralities$degree1) %in% hubsValid]
  betweenCentValid = analysis$centralities$between1[names(analysis$centralities$between1) %in% hubsValid]
  closeCentValid = analysis$centralities$close1[names(analysis$centralities$close1) %in% hubsValid]
  
  val_result = data.frame(measure = measure, zeroMethod = zeroMethod, normMethod = normMethod, 
                          sparsMethod = sparsMethod, 
                          hubsDiscov = I(list(hubsDiscov)), hubsValid = I(list(hubsValid)),
                          nrHubsDiscov = nrHubsDiscov, nrHubsValid = nrHubsValid, jaccard = jaccard,
                          cosineFamilies = cosineFamilies,
                          degreeCentDiscov = degreeCentDiscov, degreeCentValid = mean(degreeCentValid),
                          betweenCentDiscov = betweenCentDiscov, betweenCentValid = mean(betweenCentValid),
                          closeCentDiscov = closeCentDiscov, closeCentValid = mean(closeCentValid))
  return(list(res = res, val_result = val_result))
  
}

