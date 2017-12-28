source("kmeanspp2.R")

clusterOptim <- function(matrixOOM, k_low, k_high){
  
  mat <- data.matrix(matrixOOM)
  
  cl <- parallel::makeCluster(detectCores())
  parallel::clusterExport(cl, varlist = c("mat", "kmeanspp"), envir = environment())
  clusterings <- parallel::parSapply(cl, c(k_low:k_high), function(x) kmeanspp(mat, x, nstart = 5, iter.max = 20)$cluster)
  if(!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  cl <- parallel::makeCluster(detectCores())
  parallel::clusterExport(cl, varlist = c("clusterings", "mat", "intCriteria"), envir = environment())
  DB_values <- parallel::parSapply(cl, 1:ncol(clusterings), function(x) intCriteria(mat, as.integer(clusterings[,x]), c("Davies_Bouldin")))
  if(!is.null(cl)){
    parallel::stopCluster(cl)
    cl <- c()
  }
  
  return(as.integer(clusterings[, which.min(DB_values)]))
}
