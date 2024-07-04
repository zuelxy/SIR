Kotz_densty <- function(xx, mu, Sigma, r=0.5, s=2, M=3){
  ## xx, mu should be d-dimensional column vector
  xx <- as.matrix(xx)
  mu <- as.matrix(mu)
  invsigma <- solve(Sigma)
  term1 <- t(xx-mu)%*%invsigma%*%(xx-mu)
  term2 <- (det(Sigma))^(-0.5)*term1^(M-1)
  term3 <- exp(-r*(term1^s))
  dKotz <- term2*term3
  return(dKotz)
}

mu <- rep(0,4)
sigma_kotz <- matrix(c(5.3, 0, 0, -0.2, 0, 4, -0.4, 0.3,
                       0, -0.4, 6.8, 0, -0.2, 0.3, 0, 9), ncol=4, byrow = T)
N <- 2000
n <- 400

library(mvtnorm)

kotzmean_sir <- kotzmean_ansir <- kotzmean_lhssir <- matrix(ncol = 4, nrow = 1000)

set.seed(947)
for (k in 1:1000) {
  x <- rmvnorm(N, mu, sigma_kotz)
  dkotz_x <- apply(x, 1, Kotz_densty, mu=mu, Sigma=sigma_kotz)
  dnorm_x <- apply(x, 1, dmvnorm, mean=mu, sigma=sigma_kotz)
  
  ratio <- dkotz_x/dnorm_x
  quanzhong <- ratio/sum(ratio)
  
  x_sir <- resampling(x, quanzhong, n)
  x_ansir <- an_resampling(x, quanzhong, n)
  x_lsir <- lhs_resampling(x, quanzhong, n)
  
  kotzmean_sir[k, ] <- apply(x_sir, 2, mean)
  kotzmean_ansir[k, ] <- apply(x_ansir, 2, mean)
  kotzmean_lhssir[k, ] <- apply(x_lsir, 2, mean)
}

mse_sir <- apply(kotzmean_sir^2, 2, mean)
### 0.00889729 0.00681040 0.01227395 0.01560597
### 0.04358761 
mse_ansir <- apply(kotzmean_ansir^2, 2, mean)
### 0.008949697 0.006706782 0.011696862 0.015264945
### 0.04261829  (1.022744)
mse_lhssir <- apply(kotzmean_lhssir^2, 2, mean)
### 0.006330384 0.004795368 0.009019552 0.011485558
### 0.03163086  (1.378009)


