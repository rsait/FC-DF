# Data generation model from  P. Ashtari, F. N. Haredasht, H. Beigy, Supervised fuzzy partitioning, Pattern Recognition, (2020), 97, Article 107013,
generate <- function(n, seed=1366)
{
  # n: number of individual sampled
  mu1 <- c(0,0)
  S1 <- diag(c(15, 0.05))
  mu2 <- c(-12, 0)
  S2 <- diag(c(1, 1))
  mu3.1 <- c(0,8)
  S3.1 <- 4*S2
  mu3.2 <- c(0, -4)
  S3.2 <- diag(c(1, 1))
  
  n1 <- floor(n/4)
  n2 <- floor(n/4)
  n3.1 <- floor(n/3)
  n3.2 <- n- n1 - n2 -n3.1 
  set.seed(seed)
  x1 <- mnormt::rmnorm(n = n1, mean = mu1, varcov=S1, sqrt=NULL) 
  x2 <- mnormt::rmnorm(n = n2, mean = mu2, varcov=S2, sqrt=NULL) 
  x3.1 <- mnormt::rmnorm(n = n3.1, mean = mu3.1, varcov=S3.1, sqrt=NULL) 
  x3.2 <- mnormt::rmnorm(n = n3.2, mean = mu3.2, varcov=S3.2, sqrt=NULL) 
  x <- rbind(x1, x2, x3.1, x3.2)
  
  y <- rep(1:3, c(n1,n2, n3.1+n3.2))
  return(list(x, y))
}
