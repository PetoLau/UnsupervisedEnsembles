source("myBootstrap.R")
# Methods are implemented just for double seasonal time series with daily and weekly seasonalities

# Simple Methods ----
simRpart <- function(Y, K = 2, freq = 48){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  N <- length(Y)
  daily <- rep(1:seas, N/seas)
  week <- N/(seas*7)
  weekly <- rep(rep(c(1:7), week), each = seas)
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             Daily = daily,
                             Weekly = weekly)
  
  tree_1 <- rpart(Load ~ ., data = matrix_train,
                  control = rpart.control(minsplit = 2,
                                          maxdepth = 30,
                                          cp = 0.000001))
  
  # new data and prediction
  pred_tree <- predict(tree_1, matrix_train[1:freq]) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simCtreeLag <- function(Y, freq = 48, h = 48, K = 2){
  
  N <- length(Y)
  window <- (N / freq) - 1
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                freq.Ts = train_s)
  
  test_s <- tail(decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1], h)
  
  matrix_test <- data.table(fuur_test,
                            freq.Ts = test_s)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_trigonom,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.925,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, matrix_test) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simCtreeFur <- function(Y, freq = 48, h = 48, K = 4){
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.table(fourier(data_ts, K = c(K, K*2), h = freq))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  tree_2 <- party::ctree(Load ~ ., data = matrix_train,
                         controls = party::ctree_control(teststat = c("quad"),
                                                         testtype = c("Teststatistic"),
                                                         mincriterion = 0.925,
                                                         minsplit = 1,
                                                         minbucket = 1,
                                                         mtry = 0, maxdepth = 0))
  
  pred_tree <- predict(tree_2, fuur_test) + mean(trend_for)
  
  return(as.vector(pred_tree))
}

simEs <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  pred <- es(data_ts, model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h,
             cfType = "MAE", intervals = "none", silent = "all")$forecast
  
  return(as.vector(pred))
}

simSTLArima <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  stl.decom <- stl(data_ts, s.window = "periodic", robust = T)
  
  pred <- forecast(stl.decom, method = "arima", h = h)$mean
  
  return(as.vector(pred))
}

simSTLExp <- function(Y, freq = 48, h = 48) {
  
  data_ts <- ts(Y, freq = freq * 7)
  
  stl.decom <- stl(data_ts, s.window = "periodic", robust = T)
  
  pred <- forecast(stl.decom, method = "ets", h = h)$mean
  
  return(as.vector(pred))
}

# All simple methods in one FUN ----

predSimAll <- function(Y){
  
  pred_rpart <- simRpart(Y)
  pred_ctree_lag <- simCtreeLag(Y)
  pred_ctree_fur <- simCtreeFur(Y)
  pred_arima <- simSTLArima(Y)
  pred_ets <- simSTLExp(Y)
  pred_es <- simEs(Y)

  return(list(RPART = pred_rpart,
              CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              ARIMA = pred_arima,
              ETS = pred_ets,
              ES = pred_es))
}

# Bootstrapped agg. base methods ----
baggRpart <- function(Y, ntrees = 100, freq = 48, h = 48){
  
  # creation of training dataset
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  N <- length(Y)
  daily <- rep(1:seas, N/seas)
  week <- N/(seas*7)
  weekly <- rep(rep(c(1:7), week), each = seas)
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             Daily = daily,
                             Weekly = weekly)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = freq)
  for(i in 1:ntrees){
    matrixSam <- matrix_train[sample(1:N, floor(N * sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- rpart(Load ~ ., data = matrixSam,
                      control = rpart.control(minsplit = sample(2:3, 1),
                                              maxdepth = sample(27:30, 1),
                                              cp = sample(seq(0.0000009, 0.00001, by = 0.0000001), 1)))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, matrix_train[1:h]) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggCtreeLag <- function(Y, ntrees = 100, freq = 48, h = 48, K = 2){
  
  # creation of training dataset
  N <- length(Y)
  window <- (N / freq) - 1
  
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K))
  fuur_test <- as.data.frame(fourier(data_ts, K = c(K, K), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, h)$mean)
  
  train_s <- decom_1$time.series[1:(freq*window), 1]
  
  matrix_trigonom <- data.table(Load = tail(new_load, window*freq),
                                fuur[(freq+1):N,],
                                freq.Ts = train_s)
  
  test_s <- tail(decom_1$time.series[((freq*window)+1):(freq*(window+1)), 1], h)
  
  matrix_test <- data.table(fuur_test,
                            freq.Ts = test_s)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = freq)
  for(i in 1:ntrees){
    matrixSam <- matrix_trigonom[sample(1:nrow(matrix_trigonom), floor(nrow(matrix_trigonom)*sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- party::ctree(Load ~ ., data = matrixSam,
                             controls = party::ctree_control(teststat = c("quad"),
                                                             testtype = c("Teststatistic"),
                                                             mincriterion = sample(seq(0.88, 0.97, by = 0.005), 1),
                                                             minsplit = 1,
                                                             minbucket = 1,
                                                             mtry = 0, maxdepth = 0))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, matrix_test) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggCtreeFur <- function(Y, ntrees = 100, freq = 48, h = 48, K = 4){
  
  # creation of training dataset
  data_ts <- msts(Y, seasonal.periods = c(freq, freq*7))
  
  fuur <- fourier(data_ts, K = c(K, K*2))
  fuur_test <- as.data.table(fourier(data_ts, K = c(K, K*2), h = h))
  
  data_ts <- ts(Y, freq = freq*7)
  decom_1 <- stl(data_ts, s.window = "periodic", robust = T)
  new_load <- rowSums(decom_1$time.series[, c(1,3)])
  trend_part <- ts(decom_1$time.series[,2])
  
  trend_fit <- auto.arima(trend_part)
  trend_for <- as.vector(forecast(trend_fit, freq)$mean)
  
  matrix_train <- data.table(Load = new_load,
                             fuur)
  
  # model fit
  pred_mat <- matrix(0, nrow = ntrees, ncol = freq)
  for(i in 1:ntrees){
    matrixSam <- matrix_train[sample(1:nrow(matrix_train), floor(nrow(matrix_train) * sample(seq(0.7, 0.9, by = 0.01), 1)), replace = TRUE)]
    tree_bag <- party::ctree(Load ~ ., data = matrixSam,
                             controls = party::ctree_control(teststat = c("quad"),
                                                             testtype = c("Teststatistic"),
                                                             mincriterion = sample(seq(0.88, 0.97, by = 0.005), 1),
                                                             minsplit = 1,
                                                             minbucket = 1,
                                                             mtry = 0, maxdepth = 0))
    
    # new data and prediction
    pred_mat[i,] <- predict(tree_bag, fuur_test) + mean(trend_for)
  }
  
  return(pred_mat)
}

baggEs <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- bld.mbb.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  mat_pred <- sapply(ind.good, function(i)
      es(data.b[[i]], model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h, cfType = "MAE", intervals = "none", silent = "all")$forecast)
  
  return(t(mat_pred))
}

mybaggEs <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- smo.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  mat_pred <- sapply(ind.good, function(i)
    es(data.b[[i]], model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h, cfType = "MAE", intervals = "none", silent = "all")$forecast)
  
  return(t(mat_pred))
}

KMbaggEs <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- KMboot(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  mat_pred <- sapply(ind.good, function(i)
    es(data.b[[i]], model = c("AAA", "ANA", "AAdA"), ic = c("AIC"), h = h, cfType = "MAE", intervals = "none", silent = "all")$forecast)
  
  return(t(mat_pred))
}

baggSTLArima <- function(Y, ntimes = 100, freq = 48, h = 48){

  data_ts <- ts(Y, freq = freq*7)
  data.b <- bld.mbb.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "arima", h = h)$mean)
  
  return(t(mat_pred))
}

baggSTLExp <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- bld.mbb.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "ets", h = h)$mean)
  
  return(t(mat_pred))
}

mybaggSTLArima <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- smo.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "arima", h = h)$mean)
  
  return(t(mat_pred))
}

mybaggSTLExp <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- smo.bootstrap(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "ets", h = h)$mean)
  
  return(t(mat_pred))
}

KMbaggSTLArima <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- KMboot(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "arima", h = h)$mean)
  
  return(t(mat_pred))
}

KMbaggSTLExp <- function(Y, ntimes = 100, freq = 48, h = 48){
  
  data_ts <- ts(Y, freq = freq*7)
  data.b <- KMboot(data_ts, ntimes)
  
  res.inf <- sapply(1:ntimes, function(i) sum(as.integer(is.infinite(data.b[[i]]))))
  
  ind.good <- which(res.inf == 0)
  
  stl.decom <- lapply(ind.good, function(i) stl(data.b[[i]], s.window = "periodic", robust = T))
  
  mat_pred <- sapply(seq_along(stl.decom), function(i) forecast(stl.decom[[i]], method = "ets", h = h)$mean)
  
  return(t(mat_pred))
}

# All bagging functions in one FUN ----

predAllBagg <- function(Y, n_bagg = 100){
  
  pred_rpart <- baggRpart(Y, n_bagg)
  pred_ctree_lag <- baggCtreeLag(Y, n_bagg)
  pred_ctree_fur <- baggCtreeFur(Y, n_bagg)
  pred_arima <- baggSTLArima(Y, n_bagg)
  pred_ets <- baggSTLExp(Y, n_bagg)
  pred_es <- baggEs(Y, n_bagg)

  return(list(RPART = pred_rpart, CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              ARIMA = pred_arima, ETS = pred_ets, ES = pred_es))
}

# All bagging functions (with my smothed boostrapping) in one FUN ----
predAllBaggMy <- function(Y, n_bagg = 100){
  
  pred_rpart <- baggRpart(Y, n_bagg)
  pred_ctree_lag <- baggCtreeLag(Y, n_bagg)
  pred_ctree_fur <- baggCtreeFur(Y, n_bagg)
  pred_arima <- mybaggSTLArima(Y, n_bagg)
  pred_ets <- mybaggSTLExp(Y, n_bagg)
  pred_es <- mybaggEs(Y, n_bagg)

  return(list(RPART = pred_rpart, CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              ARIMA = pred_arima, ETS = pred_ets, ES = pred_es))
}

# All bagging functions (with my K-means boostrapping) in one FUN ----
predAllBaggKM <- function(Y, n_bagg = 100){
  
  pred_rpart <- baggRpart(Y, n_bagg)
  pred_ctree_lag <- baggCtreeLag(Y, n_bagg)
  pred_ctree_fur <- baggCtreeFur(Y, n_bagg)
  pred_arima <- KMbaggSTLArima(Y, n_bagg)
  pred_ets <- KMbaggSTLExp(Y, n_bagg)
  pred_es <- KMbaggEs(Y, n_bagg)

  return(list(RPART = pred_rpart, CTREE_LAG = pred_ctree_lag,
              CTREE_FUR = pred_ctree_fur,
              ARIMA = pred_arima, ETS = pred_ets, ES = pred_es))
}

# Helper function
check.matrix <- function(mat){
  if(is.matrix(mat) == TRUE){
    return(mat)
  } else
    return(t(as.matrix(mat)))
}

# Unsupervised Ensemble Learning Methods ----

predEnsembles <- function(Y, n_bagg = 100){
  
  # compute all predictions from above functions
  predictions <- predAllBagg(Y, n_bagg)
  
  # transformation to data.matrix
  pred_matrix <- data.matrix(rbindlist(lapply(1:length(predictions), function(x) as.data.table(predictions[[x]]))))
  
  # Compute first 3 PC
  if (sd(Y) == 0){
    pred_pc <- prcomp(pred_matrix, center = T, scale. = F)$x[, 1:3]
  } else {
    pred_pc <- prcomp(pred_matrix, center = T, scale. = T)$x[, 1:3]
  }
  
  if (nrow(unique(pred_pc)) == 1){
    pred_clus_res <- rep(1, nrow(pred_matrix))
  } else {
    # Optimal clustering by DB index and K-means++
    pred_clus_res <- clusterOptimKmeans(pred_pc, 3, 8)
  }
  
  # DBSCAN
  opt_res <- dbscan::optics(pred_pc, eps = 2.5, minPts = 8)
  
  db_res <- extractDBSCAN(opt_res, eps_cl = 1.45)$cluster

  xi_res <- extractXi(opt_res, xi = 0.045, minimum = F, correctPredecessors = T)

  # Final predictions
  ens_stack_pred_km <- rowMeans(sapply(unique(pred_clus_res), function(x) colMeans(check.matrix(pred_matrix[pred_clus_res == x,]))))
  
  if(max(unique(db_res)) == 0) {
    ens_stack_pred_db <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db <- rowMeans(sapply(1:max(unique(db_res)), function(x) apply(check.matrix(pred_matrix[db_res == x,]), 2, median)))
  }
  
  if(is.null(xi_res$cluster)) {
    ens_stack_pred_xi_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_med <- apply(sapply(1:max(unique(xi_res$cluster)), function(x) apply(check.matrix(pred_matrix[xi_res$cluster == x,]), 2, median)), 1, median)
  }
  
  ens_medians <- rowMeans(sapply(1:length(predictions), function(i) apply(predictions[[i]], 2, median)))
  
  return(list(RPART.median = apply(predictions$RPART, 2, median),
              CTREE_LAG.median = apply(predictions$CTREE_LAG, 2, median),
              CTREE_FUR.median = apply(predictions$CTREE_FUR, 2, median),
              ARIMA.median = apply(predictions$ARIMA, 2, median),
              ETS.median = apply(predictions$ETS, 2, median),
              ES.median = apply(predictions$ES, 2, median),
              Simple.mean = apply(pred_matrix, 2, mean),
              Simple.median = apply(pred_matrix, 2, median),
              Ens.medians = ens_medians,
              Ensemble_km = ens_stack_pred_km,
              Ensemble_db = ens_stack_pred_db,
              Ensemble_xi_med = ens_stack_pred_xi_med))
}

predEnsemblesMy <- function(Y, n_bagg = 100){
  
  # compute all predictions from above functions
  predictions <- predAllBaggMy(Y, n_bagg)
  
  # transformation to data.matrix
  pred_matrix <- data.matrix(rbindlist(lapply(1:length(predictions), function(x) as.data.table(predictions[[x]]))))
  
  # Compute first 3 PC
  if (sd(Y) == 0){
    pred_pc <- prcomp(pred_matrix, center = T, scale. = F)$x[, 1:3]
  } else {
    pred_pc <- prcomp(pred_matrix, center = T, scale. = T)$x[, 1:3]
  }
  
  if (nrow(unique(pred_pc)) == 1){
    pred_clus_res <- rep(1, nrow(pred_matrix))
  } else {
    # Optimal clustering by DB index and K-means++
    pred_clus_res <- clusterOptimKmeans(pred_pc, 3, 8)
  }
  
  # DBSCAN
  opt_res <- dbscan::optics(pred_pc, eps = 2.5, minPts = 8)

  db_res <- extractDBSCAN(opt_res, eps_cl = 1.45)$cluster

  xi_res <- extractXi(opt_res, xi = 0.045, minimum = F, correctPredecessors = T)

  # Final predictions
  ens_stack_pred_km <- rowMeans(sapply(unique(pred_clus_res), function(x) colMeans(check.matrix(pred_matrix[pred_clus_res == x,]))))
  
  if(max(unique(db_res)) == 0) {
    ens_stack_pred_db <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db <- rowMeans(sapply(1:max(unique(db_res)), function(x) apply(check.matrix(pred_matrix[db_res == x,]), 2, median)))
  }
  
  if(is.null(xi_res$cluster)) {
    ens_stack_pred_xi_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_med <- apply(sapply(1:max(unique(xi_res$cluster)), function(x) apply(check.matrix(pred_matrix[xi_res$cluster == x,]), 2, median)), 1, median)
  }
  
  ens_medians <- rowMeans(sapply(1:length(predictions), function(i) apply(predictions[[i]], 2, median)))
  
  return(list(RPART.median = apply(predictions$RPART, 2, median),
              CTREE_LAG.median = apply(predictions$CTREE_LAG, 2, median),
              CTREE_FUR.median = apply(predictions$CTREE_FUR, 2, median),
              ARIMA.median = apply(predictions$ARIMA, 2, median),
              ETS.median = apply(predictions$ETS, 2, median),
              ES.median = apply(predictions$ES, 2, median),
              Simple.mean = apply(pred_matrix, 2, mean),
              Simple.median = apply(pred_matrix, 2, median),
              Ens.medians = ens_medians,
              Ensemble_km = ens_stack_pred_km,
              Ensemble_db = ens_stack_pred_db,
              Ensemble_xi_med = ens_stack_pred_xi_med))
}

predEnsemblesKM <- function(Y, n_bagg = 100){
  
  # compute all predictions from above functions
  predictions <- predAllBaggKM(Y, n_bagg)
  
  # transformation to data.matrix
  pred_matrix <- data.matrix(rbindlist(lapply(1:length(predictions), function(x) as.data.table(predictions[[x]]))))
  
  # Compute first 3 PC
  if (sd(Y) == 0){
    pred_pc <- prcomp(pred_matrix, center = T, scale. = F)$x[, 1:3]
  } else {
    pred_pc <- prcomp(pred_matrix, center = T, scale. = T)$x[, 1:3]
  }
  
  if (nrow(unique(pred_pc)) == 1){
    pred_clus_res <- rep(1, nrow(pred_matrix))
  } else {
    # Optimal clustering by DB index and K-means++
    pred_clus_res <- clusterOptimKmeans(pred_pc, 3, 8)
  }
  
  # DBSCAN
  opt_res <- dbscan::optics(pred_pc, eps = 2.5, minPts = 8)

  db_res <- extractDBSCAN(opt_res, eps_cl = 1.45)$cluster

  xi_res <- extractXi(opt_res, xi = 0.045, minimum = F, correctPredecessors = T)

  # Final predictions
  ens_stack_pred_km <- rowMeans(sapply(unique(pred_clus_res), function(x) colMeans(check.matrix(pred_matrix[pred_clus_res == x,]))))
  
  if(max(unique(db_res)) == 0) {
    ens_stack_pred_db <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_db <- rowMeans(sapply(1:max(unique(db_res)), function(x) apply(check.matrix(pred_matrix[db_res == x,]), 2, median)))
  }
  
  if(is.null(xi_res$cluster)) {
    ens_stack_pred_xi_med <- apply(pred_matrix, 2, median)
  } else {
    ens_stack_pred_xi_med <- apply(sapply(1:max(unique(xi_res$cluster)), function(x) apply(check.matrix(pred_matrix[xi_res$cluster == x,]), 2, median)), 1, median)
  }
  
  ens_medians <- rowMeans(sapply(1:length(predictions), function(i) apply(predictions[[i]], 2, median)))
  
  return(list(RPART.median = apply(predictions$RPART, 2, median),
              CTREE_LAG.median = apply(predictions$CTREE_LAG, 2, median),
              CTREE_FUR.median = apply(predictions$CTREE_FUR, 2, median),
              ARIMA.median = apply(predictions$ARIMA, 2, median),
              ETS.median = apply(predictions$ETS, 2, median),
              ES.median = apply(predictions$ES, 2, median),
              Simple.mean = apply(pred_matrix, 2, mean),
              Simple.median = apply(pred_matrix, 2, median),
              Ens.medians = ens_medians,
              Ensemble_km = ens_stack_pred_km,
              Ensemble_db = ens_stack_pred_db,
              Ensemble_xi_med = ens_stack_pred_xi_med))
}
