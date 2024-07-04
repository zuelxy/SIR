####### global function 
resampling <- function(pool, weight, n){
  ## pool: the sample pool
  ## weight: the weight corresponding to each sample in the pool
  ## n: the size of the resampling size
  pool <- as.matrix(pool)
  u <- runif(n)
  cumweight <- c(0, cumsum(weight))
  index <- findInterval(u, cumweight)
  newsample <- pool[index,]
  return(newsample)
}

an_resampling <- function(pool, weight, n){
  ## pool: the sample pool
  ## weight: the weight corresponding to each sample in the pool
  ## n: the size of the resampling size
  pool <- as.matrix(pool)
  u <- runif(n/2)
  an_u <- 1-u
  cumweight <- c(0, cumsum(weight))
  index <- findInterval(u, cumweight)
  an_index <- findInterval(an_u, cumweight)
  newsample <- pool[c(index,an_index),]
  return(newsample)
}


lhs_resampling <- function(pool, weight, n){
  pool <- as.matrix(pool)
  
  ran <- runif(n)
  u <- (1:n-ran)/n
  
  cumweight <- c(0, cumsum(weight))
  index <- findInterval(u, cumweight)
  newsample <- pool[index,]
  
  return(newsample)
}
