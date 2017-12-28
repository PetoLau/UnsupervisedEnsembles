# kmeanspp <- function(x, k, iter.max = 10, nstart = 1, ...) {
#   n <- nrow(x) # number of data points
#   centers <- numeric(k) # IDs of centers
#   distances <- matrix(numeric(n * (k - 1)), ncol = k - 1) # distances[i, j]: The distance between x[i,] and x[centers[j],]
#   res.best <- list(tot.withinss = Inf) # the best result among <nstart> iterations
#   for (rep in 1:nstart) {
#     pr <- rep(1, n) # probability for sampling centers
#     for (i in 1:(k - 1)) {
#       centers[i] <- sample.int(n, 1, prob = pr) # Pick up the ith center
#       distances[, i] <- colSums((t(x) - x[centers[i], ])^2) # Compute (the square of) distances to the center
#       pr <- distances[cbind(1:n, max.col(-distances[, 1:i, drop = FALSE]))] # Compute probaiblity for the next sampling
#     }
#     centers[k] <- sample.int(n, 1, prob = pr)
#     data_centers <- as.matrix(x[centers,])
#     variations <- sum(seq(from = (k-1), to = 1))
#     indexes <- vector(length = variations)
#     k_iter <- 1
#     for(i in 1:(nrow(data_centers)-1)) {
#       for(j in (i+1):nrow(data_centers)){
#         indexes[k_iter] <- identical(data_centers[i,], data_centers[j,])
#         k_iter <- k_iter + 1
#       }
#     }
#     ## Perform k-means with the obtained centers
#     if(sum(indexes) > 0) {
#       res <- kmeans(x, k, iter.max = iter.max, nstart = 1, ...)
#     } else {
#       res <- kmeans(x, data_centers, iter.max = iter.max, nstart = 1, ...)
#       res$inicial.centers <- data_centers
#     }
#     
#     ## Store the best result
#     if (res$tot.withinss < res.best$tot.withinss) {
#       res.best <- res
#     }
#   }
#   res.best
# }

kmeanspp <- function(x, k, iter.max = 10, nstart = 1, ...) {
  
  res <- kmeans(x, k, iter.max = iter.max, nstart = nstart)
  
  res
}
