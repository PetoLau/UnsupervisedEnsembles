# Automatic forecasting with sliding window approach ----
source("ensembles.R")

pickwindow <- function(data, N.win){
  
  data_sub <- data[N.slid.win %in% N.win, lapply(Load, as.vector), by = N.K]
  data_sub[, N.K := NULL]
  data_sub <- data.matrix(data_sub)
  
  data_sub
}

ForecastClustersEnsemble <- function(dataset, FUN = predEnsembles){
  
  n_day <- length(unique(dataset$N.slid.win))
  
  data_list <- lapply(0:(n_day-1), function(x) pickwindow(dataset, x))
  
  cl <- makeForkCluster(detectCores()-1, outfile = "debug.txt")
  registerDoParallel(cl)
  pred_clusters <- parLapply(cl, 1:(length(data_list)),
                             function(i) lapply(1:dim(data_list[[i]])[1], function(j)
                               FUN(data_list[[i]][j,], 100)))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  res <- sapply(seq_len(length(pred_clusters[[1]][[1]])),
                function(k) as.vector(sapply(seq_len(length(pred_clusters)),
                                             function(i) rowSums(sapply(seq_len(length(pred_clusters[[i]])),
                                                                        function(j) as.vector(pred_clusters[[i]][[j]][[k]]))))))
  
  return(res)
}

ForecastClustersSimple <- function(dataset){
  
  n_day <- length(unique(dataset$N.slid.win))
  
  data_list <- lapply(0:(n_day-1), function(x) pickwindow(dataset, x))
  
  pred_clusters <- lapply(1:(length(data_list)),
                          function(i) lapply(1:dim(data_list[[i]])[1], 
                              function(j) predSimAll(data_list[[i]][j,])))

  res <- sapply(seq_len(length(pred_clusters[[1]][[1]])),
                function(k) as.vector(sapply(seq_len(length(pred_clusters)),
                                             function(i) rowSums(sapply(seq_len(length(pred_clusters[[i]])),
                                                                        function(j) as.vector(pred_clusters[[i]][[j]][[k]]))))))
  
  return(res)
}

ForecastAggregatedEnsemble <- function(dataset, FUN = predEnsembles){
  
  win <- 21
  days_for <- length(dataset)/seas - win
  
  cl <- makeForkCluster(detectCores()-1, outfile = "debug.txt")
  registerDoParallel(cl)
  pred_sums <- parLapply(cl, 0:(days_for-1), function(i) FUN(dataset[((i*seas)+1):((seas*i)+(win*seas))]))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  predictions <- sapply(seq_len(length(pred_sums[[1]])),
                        function(k) sapply(seq_len(days_for),
                                           function(j) as.vector(pred_sums[[j]][[k]])))
  
  return(predictions)
}

ForecastAggregatedSimple <- function(dataset){
  
  win <- 21
  days_for <- length(dataset)/seas - win
  
  pred_sums <- lapply(0:(days_for-1),
                      function(i) predSimAll(dataset[((i*seas)+1):((seas*i)+(win*seas))]))
  
  predictions <- sapply(seq_len(length(pred_sums[[1]])),
                        function(k) sapply(seq_len(days_for),
                                           function(j) as.vector(pred_sums[[j]][[k]])))
  
  return(predictions)
}

# Compute MAPE ----
computeMape <- function(real, predictions){
  
  win <- 21
  data_test <- real[-(1:(win*seas))]
  n_day <- length(data_test)/seas
  
  err_byday <- sapply(seq_len(dim(predictions)[2]),
                      function(i) sapply(0:(n_day-1),
                        function(j) mape(as.vector(data_test[((j*seas)+1):(seas*(j+1))]), predictions[((j*seas)+1):(seas*(j+1)), i])))
  
  
  err_whole <- sapply(seq_len(dim(predictions)[2]),
                      function(i) mape(as.vector(data_test), predictions[, i]))
  
  return(list(ByDay = err_byday, Whole = matrix(err_whole)))
}
