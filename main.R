# Header ----
# gc()

# devtools::install_github("PetoLau/TSrepr") # https://github.com/PetoLau/TSrepr
pkgs <- c("parallel", "doParallel", "data.table", "smooth", "forecast",
          "clusterCrit", "rpart", "party", "ggplot2", "dbscan", "TSrepr")
lapply(pkgs, library, character.only = TRUE) # load them

# Read time series data ----
# rows of a matrix must represent time series
data_ts <- fread("file.csv")

seas <- 48 # seasonality

# Dynamic clustering ----
source("optimalClustering.R")

data.clust <- data.frame(N.slid.win = 0, N.K = 0, Load = 0)
clust.res <- data.frame(N.slid.win = 0, Class = 0, ID = 0)

SetTrain_sum <- colSums(data_ts)
win <- 21
days_for <- length(SetTrain_sum)/seas - win

for(k in 0:(days_for-1)){

  oomDT.select <- data_ts[, ((k*seas)+1):((k+win)*seas), with = F]
  oomDT.sel.scale <- t(apply(oomDT.select, 1, norm_z))

  repr_z <- repr_matrix(oomDT.sel.scale, func = repr_lm, args = list(method = "lm", freq = c(seas, seas*7)))

  km_res <- clusterOptim(repr_z, 10, 18)
  
  km_sums <- t(sapply(1:length(unique(km_res)), function(x) colSums(oomDT.select[km_res == x,])))
  
  for(l in 1:length(unique(km_res))){
    data.clust <- rbind(data.clust,
                        data.frame(N.slid.win = k, N.K = l, Load = km_sums[l,]))
  }
  
  clust.res <- rbind(clust.res,
                     data.frame(N.slid.win = k, Class = km_res, ID = 1:nrow(oomDT.select)))
  
  print(k)
}

data.clust <- data.clust[-1,]
clust.res <- clust.res[-1,]
summary(data.clust)
summary(clust.res)

write.table(data.clust, "LM_data_clust.csv", row.names = F, col.names = T, quote = F)
write.table(clust.res, "LM_data_clust_info.csv", row.names = F, col.names = T, quote = F)

# data.clust_1 <- fread("LM_second_Irske.csv")

max_K <- sapply(0:(days_for-1), function(x) max(data.clust[data.clust[,1] %in% x, 2]))
table(max_K)

# Testing ----
# aggregated time series ensemble learning forecasting ----
source("TestingPredictions.R")

# Aggregate time series
data_sum <- colSums(data_ts)

res_ens_1 <- ForecastAggregatedEnsemble(data_sum, FUN = predEnsembles)
err_ens_1 <- computeMape(data_sum, res_ens_1)
gc()

res_ens_my_1 <- ForecastAggregatedEnsemble(data_sum, FUN = predEnsemblesMy)
err_ens_my_1 <- computeMape(data_sum, res_ens_my_1)
gc()

res_ens_km_1 <- ForecastAggregatedEnsemble(data_sum, FUN = predEnsemblesKM)
err_ens_km_1 <- computeMape(data_sum, res_ens_km_1)
gc()

# Compare
c(err_ens_1$Whole)
c(err_ens_my_1$Whole)
c(err_ens_km_1$Whole)

write.table(res_ens_1, "res_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_my_1, "res_ens_my.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_km_1, "res_ens_km.csv", row.names = F, col.names = F, quote = F)

# Simple forecasting of aggregated time series
res_sim <- ForecastAggregatedSimple(data_sum)
err_sim <- computeMape(data_sum, res_sim)
c(err_sim$Whole)

write.table(res_sim, "res_sim.csv", row.names = F, col.names = F, quote = F)

# Forecasting on clusters - ensemble learning ----
# read clustered data
data.clust <- fread("LM_data_clust.csv")
data_sum <- colSums(data_ts)

res_ens_clus <- ForecastClustersEnsemble(data.clust, FUN = predEnsembles)
gc()
res_ens_clus_my <- ForecastClustersEnsemble(data.clust, FUN = predEnsemblesMy)
gc()
res_ens_clus_km <- ForecastClustersEnsemble(data.clust, FUN = predEnsemblesKM)
gc()

write.table(res_ens_clus, "res_lm_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_clus_my, "res_lm_my_ens.csv", row.names = F, col.names = F, quote = F)
write.table(res_ens_clus_km, "res_lm_km_ens.csv", row.names = F, col.names = F, quote = F)

# simple forecasting on clusters ----
res_sim_clust <- ForecastClustersSimple(data.clust)
err_sim_clust <- computeMape(data_sum, res_sim_clust)
gc()

c(err_sim_clust$Whole)

write.table(res_sim_clust, "res_lm_sim.csv", row.names = F, col.names = F, quote = F)
