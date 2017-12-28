# Moving Block Bootstrap ----
source("kmeanspp2.R")
# source code from *forecast* package - MBB

MBB <- function(x, window_size) {

  # window_size <- 2*seas*7
  # x <- rnorm(seas*7*3)
  #
  bx <- array(0, (floor(length(x)/window_size)+2)*window_size)
  for (i in 1:(floor(length(x)/window_size)+2)){
    c <- sample(1:(length(x)-window_size+1),1)
    bx[((i-1)*window_size+1):(i*window_size)] <- x[c:(c+window_size-1)]
  }
  start_from <- sample(0:(window_size-1),1) + 1
  bx[start_from:(start_from+length(x)-1)]
}

# My smoothed bootstrap version of previous bld.mbb.bootstrap function ----
smo.bootstrap <- function(x, num, block_size=NULL, alpha = 0.05)
{
  # x <- ts(some_load, freq = seas*7)
  freq <- frequency(x)
  if(is.null(block_size))
  {
    block_size <- ifelse(freq > 1, 2*freq, min(8, floor(length(x)/ 2)))
  }

  xs <- list()
  xs[[1]] <- x # the first series is the original one

  if (num>1){
    # Box-Cox transformation
    lambda <- BoxCox.lambda(x, lower=0.01, upper=1)
    x.bc <- BoxCox(x, lambda)

    if (freq>1){
      # STL decomposition
      x.stl <- stl(ts(x.bc, frequency=freq), "per", robust = TRUE)$time.series
      seasonal <- x.stl[,1]
      trend <- x.stl[,2]
      remainder <- ts(c(0, HoltWinters(x.stl[,3], alpha = alpha, beta = FALSE, gamma = FALSE)$fitted[,1]), freq = freq)
    
    } else {
      # Loess
      trend <- 1:length(x)
      suppressWarnings(x.loess <- loess(x.bc ~ trend, span=6/length(x), degree=1))
      seasonal <- rep(0, length(x))
      trend <- x.loess$fitted
      remainder <- x.loess$residuals
    }

    # Bootstrap some series, using MBB
    for (i in 2:num){
      xs[[i]] <- InvBoxCox(trend + seasonal + MBB(remainder, block_size), lambda)
    }
  }

  xs
}

# K-means based bootstrap ----

clusterOptimKmeans <- function(matrixOOM, k_low, k_high){
  
  mat <- data.matrix(matrixOOM)
  
  clusterings <- sapply(c(k_low:k_high), function(x) kmeanspp(mat, x, nstart = 10, iter.max = 20)$cluster)
  
  DB_values <- sapply(1:ncol(clusterings), function(x) intCriteria(mat, as.integer(clusterings[,x]), c("Davies_Bouldin")))
  
  return(as.integer(clusterings[, which.min(DB_values)]))
}

KMboot <- function(x, num = 100) {
  
  freq <- frequency(x)
  datas <- data.matrix(x)
  
  km_res <- clusterOptimKmeans(datas, 12, 20)

  xs <- list()
  xs[[1]] <- ts(x, freq = freq)
  
  for(j in 2:num){
    xs[[j]] <- vector(length = length(x))
    for(i in 1:length(x)) {
      xs[[j]][i] <- sample(x[km_res %in% km_res[i]], 1)
    }
    xs[[j]] <- ts(xs[[j]], freq = freq)
  }
  
  return(xs)
}
