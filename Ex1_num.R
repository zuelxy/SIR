N <- 20000
n <- 1000

########
# U(0,1) Be(2,3)
set.seed(972)
case1_sir <- case1_ansir  <- case1_lhssir <- numeric()
for(k in 1:1000){
  x <- runif(N)
  ratio <- dbeta(x, 2, 3)/dunif(x)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case1_sir[k] <- mean(x_star)
  case1_ansir[k] <- mean(x_astar)
  case1_lhssir[k] <- mean(x_lstar)
}

mean((case1_sir-0.4)^2)
#  [1] 4.189533e-05
mean((case1_ansir-0.4)^2)
# 【1] 4.006772e-05  (1.045613)
mean((case1_lhssir-0.4)^2)
#  [1] 3.989958e-05  (1.050019)


##########
# U(0,1) Be(0.9, 0.9)
set.seed(81)
case2_sir <- case2_ansir  <- case2_lhssir  <- numeric()
for(k in 1:1000){
  x <- runif(N)
  ratio <- dbeta(x, 0.9, 0.9)/dunif(x)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case2_sir[k] <- mean(x_star)
  case2_ansir[k] <- mean(x_astar)
  case2_lhssir[k] <- mean(x_lstar)
}

mean((case2_sir-0.5)^2)
# [1] 0.0001036633 
mean((case2_ansir-0.5)^2)
# [1] 9.612684e-05 (1.078401)
mean((case2_lhssir-0.5)^2)
# [1] 9.022977e-05 (1.148881)

############
# log(0, 1) N(0,1)
set.seed(86)
case3_sir <- case3_ansir  <- case3_lhssir  <- numeric()
for(k in 1:1000){
  x <- rlogis(N)
  ratio <- dnorm(x, 0, 1)/dlogis(x)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case3_sir[k] <- mean(x_star)
  case3_ansir[k] <- mean(x_astar)
  case3_lhssir[k] <- mean(x_lstar)
}

mean((case3_sir-0)^2)
# [1] 0.001144058
mean((case3_ansir-0)^2)
# [1] 0.001070271 (1.068942)
mean((case3_lhssir-0)^2)
# [1] 0.001053291 (1.086175)


###########
# C(0,1) N(0,1)
set.seed(35)
case4_sir <- case4_ansir  <- case4_lhssir  <- numeric()
for(k in 1:1000){
  x <- rcauchy(N)
  ratio <- dnorm(x, 0, 1)/dcauchy(x)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case4_sir[k] <- mean(x_star)
  case4_ansir[k] <- mean(x_astar)
  case4_lhssir[k] <- mean(x_lstar)
}

mean((case4_sir-0)^2)
# [1] 0.001136371
mean((case4_ansir-0)^2)
# [1] 0.0009898121 (1.148067)
mean((case4_lhssir-0)^2)
# [1] 0.001068776  (1.063245)

########### 
# C(0,1) t(2,0,1)
set.seed(819)
case5_sir <- case5_ansir  <- case5_lhssir  <- numeric()
for(k in 1:1000){
  x <- rcauchy(N)
  ratio <- dt(x, 2)/dcauchy(x)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case5_sir[k] <- mean(x_star)
  case5_ansir[k] <- mean(x_astar)
  case5_lhssir[k] <- mean(x_lstar)
}

mean((case5_sir-0)^2)
# [1] 0.02377132
mean((case5_ansir-0)^2)
# [1] 0.01478423  (1.607884)
mean((case5_lhssir-0)^2)
# [1] 0.01544801  (1.538795)


########### 
# IG(1,1) F(10, 6)
library(actuar)
set.seed(191)
case6_sir <- case6_ansir  <- case6_lhssir  <- numeric()
for(k in 1:1000){
  x <- rinvgamma(N, shape = 1, scale = 1)
  ratio <- df(x, 10, 6)/dinvgamma(x, shape = 1, scale = 1)
  quanzhong <- ratio/sum(ratio)
  
  x_star <- resampling(x, quanzhong, n)
  x_astar <- an_resampling(x, quanzhong, n)
  x_lstar <- lhs_resampling(x, quanzhong, n)
  
  case6_sir[k] <- mean(x_star)
  case6_ansir[k] <- mean(x_astar)
  case6_lhssir[k] <- mean(x_lstar)
}

mean((case6_sir-1.5)^2)
# [1] 0.003511655
mean((case6_ansir-1.5)^2)
# [1] 0.002962995  (1.185171)
mean((case6_lhssir-1.5)^2)
# [1] 0.00309055  (1.136256)
