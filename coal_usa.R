source("resampling.R")
coal <- read.table("coal_usa.txt", header = T)

coal_data <- coal[which(coal$year>1922), ]

likelihood <- function(Theta, x){
  theta <- Theta[1]
  lambda_1 <- Theta[2]
  lambda_2 <- Theta[3]
  
  d <- length(x)
  x1_t <- x[1:theta]
  xt_d <- x[(theta+1):d]
  
  term1 <- (lambda_1^x1_t)*exp(-lambda_1)/factorial(x1_t)
  term1_value <- prod(term1)
  
  term2 <- (lambda_2^xt_d)*exp(-lambda_2)/factorial(xt_d)
  term2_value <- prod(term2)
  
  return(term1_value*term2_value)
}

############################ case 1 ########################

N = 5000
n = 2000
K = 500

sir_result <- c()
antisir_result <- c()
lhssir_result <- c()

set.seed(898)
for (k in 1:K) {
  theta_s <- sample(1:99, N, replace = T)
  
  lambda_s <- c()
  as_1 <- rgamma(N, 10, 10)
  as_2 <- rgamma(N, 10, 10)
  for (i in 1:N) {
    lambdas_1 <- rgamma(1, 3, as_1[i])
    lambdas_2 <- rgamma(1, 3, as_2[i])
    lambda_s <- rbind(lambda_s, c(lambdas_1,lambdas_2))
  }
  
  Theta_1 <- cbind(theta_s, lambda_s)
  colnames(Theta_1) <- NULL
  
  weight_s <- apply(Theta_1, 1, likelihood, x=coal_data$incidents)
  quanzhong_s <- weight_s/sum(weight_s)
  
  Theta_post <- resampling(Theta_1, quanzhong_s, n)
  sir <- apply(Theta_post, 2, mean)
  
  sir_result <- rbind(sir_result, sir)
}


set.seed(112)
for (k in 1:K) {
  theta2 <- sample(1:99, N, replace = T)
  
  lambda2 <- c()
  a2_1 <- rgamma(N, 10, 10)
  a2_2 <- rgamma(N, 10, 10)
  for (i in 1:N) {
    lambda2_1 <- rgamma(1, 3, a2_1[i])
    lambda2_2 <- rgamma(1, 3, a2_2[i])
    lambda2 <- rbind(lambda2, c(lambda2_1,lambda2_2))
  }
  
  Theta_2 <- cbind(theta2, lambda2)
  colnames(Theta_2) <- NULL
  
  weight2 <- apply(Theta_2, 1, likelihood, x=coal_data$incidents)
  quanzhong2 <- weight2/sum(weight2)
  
  Theta_post_anti <- an_resampling(Theta_2, quanzhong2, n)
  anti_sir <- apply(Theta_post_anti, 2, mean)
  
  antisir_result <- rbind(antisir_result, anti_sir)
}


set.seed(192)
for (k in 1:K) {
  theta <- sample(1:99, N, replace = T)
  
  lambda <- c()
  a_1 <- rgamma(N, 10, 10)
  a_2 <- rgamma(N, 10, 10)
  for (i in 1:N) {
    lambda_1 <- rgamma(1, 3, a_1[i])
    lambda_2 <- rgamma(1, 3, a_2[i])
    lambda <- rbind(lambda, c(lambda_1,lambda_2))
  }
  
  Theta_3 <- cbind(theta, lambda)
  colnames(Theta_3) <- NULL
  
  weight <- apply(Theta_3, 1, likelihood, x=coal_data$incidents)
  quanzhong <- weight/sum(weight)
  
  Theta_post_lhs <- lhs_resampling(Theta_3, quanzhong, n)
  lhs_sir <- apply(Theta_post_lhs, 2, mean)
  
  lhssir_result <- rbind(lhssir_result, lhs_sir)
}


apply(sir_result, 2, sd)
# 1.1854365 0.7475860 0.1562243
apply(antisir_result, 2, sd)
# 1.1732644 0.7358252 0.1539009
apply(lhssir_result, 2, sd)
# 1.1608661 0.7093401 0.1491981

apply(sir_result, 2, mean)
# 26.6622270  6.6561927  0.8323431
apply(antisir_result, 2, mean)
# 26.6400410  6.6426237  0.8323192
apply(lhssir_result, 2, mean)
# 26.6529660  6.6722731  0.8229999

###################### case 2 ##############################
N = 5000
n = 2000
K = 500

sir_result <- c()
antisir_result <- c()
lhssir_result <- c()

set.seed(69)
for (k in 1:K) {
  thetas <- sample(1:99, N, replace = T)
  
  lambdas <- c()
  as <- rgamma(N, 10, 10)
  log_alphas <- runif(N, log(1/8), log(2))
  alphas <- exp(log_alphas)
  for (i in 1:N) {
    lambdas_1 <- rgamma(1, 3, as[i])
    lambdas_2 <- alphas[i]*lambdas_1
    lambdas <- rbind(lambdas, c(lambdas_1,lambdas_2))
  }
  
  Theta_1 <- cbind(thetas, lambdas)
  colnames(Theta_1) <- NULL
  
  weights <- apply(Theta_1, 1, likelihood, x=coal_data$incidents)
  quanzhongs <- weights/sum(weights)
  
  Theta_post <- resampling(Theta_1, quanzhongs, n)
  sir <- apply(Theta_post, 2, mean)
 
  sir_result <- rbind(sir_result, sir)
}


set.seed(39)
for (k in 1:K) {
  thetaa <- sample(1:99, N, replace = T)
  
  lambdaa <- c()
  aa <- rgamma(N, 10, 10)
  log_alphaa <- runif(N, log(1/8), log(2))
  alphaa <- exp(log_alphaa)
  for (i in 1:N) {
    lambdaa_1 <- rgamma(1, 3, aa[i])
    lambdaa_2 <- alphaa[i]*lambdaa_1
    lambdaa <- rbind(lambdaa, c(lambdaa_1,lambdaa_2))
  }
  
  Theta_2 <- cbind(thetaa, lambdaa)
  colnames(Theta_2) <- NULL
  
  weighta <- apply(Theta_2, 1, likelihood, x=coal_data$incidents)
  quanzhonga <- weighta/sum(weighta)
  
  Theta_post_anti <- an_resampling(Theta_2, quanzhonga, n)
  
  anti_sir <- apply(Theta_post_anti, 2, mean)
  antisir_result <- rbind(antisir_result, anti_sir)
}

set.seed(91)
for (k in 1:K) {
  theta <- sample(1:99, N, replace = T)
  
  lambda <- c()
  a <- rgamma(N, 10, 10)
  log_alpha <- runif(N, log(1/8), log(2))
  alpha <- exp(log_alpha)
  for (i in 1:N) {
    lambda_1 <- rgamma(1, 3, a[i])
    lambda_2 <- alpha[i]*lambda_1
    lambda <- rbind(lambda, c(lambda_1,lambda_2))
  }
  
  Theta_3 <- cbind(theta, lambda)
  colnames(Theta_3) <- NULL
  
  weightl <- apply(Theta_3, 1, likelihood, x=coal_data$incidents)
  quanzhongl <- weightl/sum(weightl)
  
  Theta_post_lhs <- lhs_resampling(Theta_3, quanzhongl, n)
  lhs_sir <- apply(Theta_post_lhs, 2, mean)
  
  lhssir_result <- rbind(lhssir_result, lhs_sir)
}


apply(sir_result, 2, sd)
# 1.0337767 0.6004933 0.1064548
apply(antisir_result, 2, sd)
# 1.0185208 0.5667787 0.1039985
apply(lhssir_result, 2, sd)
# 0.96216655 0.57541682 0.09545736

apply(sir_result, 2, mean)
# 26.5327255  6.2819503  0.9238379
apply(antisir_result, 2, mean)
# 26.5477630  6.2497285  0.9235687
apply(lhssir_result, 2, mean)
# 26.5265620  6.2858467  0.9246661
