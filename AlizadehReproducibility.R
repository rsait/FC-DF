# Aiming reproducibility purposes, we publish the code we used with Alizadeh data set
# Tuning for clustering
#To understand the use of the available R functions, the reader should revise example1.R and example2.R


source("functions.R")
library(dplyr)
library(fclust)

# Load and order data set
dat <- read.table("Alizadeh_v3.csv", sep=";", dec=",", header=TRUE)
# Class label
y <- factor((dat[, 1]))
ynum <- as.numeric(y)

# Remaining data
datred <- dat[, -1] 
p <- dim(datred)[2]
for (i in 1:p)
{
  datred[, i] <- as.numeric(unlist(datred[, i]))
}
x <- as.matrix(datred)

# Correlation distance
dtrain <-  as.matrix(sqrt(1 - cor(t(x), use="na.or.complete")))
n <- dim(dtrain)[1]
#-----------------------------


# Set grid for hyper-parameters (alpha=0 for clustering)
valoresgamma <- c(2^seq(-9, 2, by=1)) 
valoresK <- 2:15
dfvalores <- expand.grid(gamma = valoresgamma , K=valoresK)
dim(dfvalores)  
dimhyper <- 2

logtol <- function(x, eps)
{
  x[x < eps] <- NA
  log(x)
}

# Set number of bootstrap replications. Preliminary results are gathered in "results". Number of repetitions in "nrep"
B <- 100
init <- matrix(NA, nrow=dim(dfvalores)[1], ncol=2 + max(valoresK)) 
results <- data.frame(dfvalores, init)
names(results) <- c("gamma",  "K", paste("a", 1:max(valoresK), sep=""), "G", "meanLogGb")
nrep <- 5
for (ih in 1:dim(dfvalores)[1])
{
  cat("\nih = ", ih, "out of ", dim(dfvalores)[1])
  set.seed(13)
  haziak <- sample(1:10000, nrep, replace=FALSE)
  aux <- rep(NA, nrep)
  prov <- vector("list", nrep)
  objetivevalues <- rep(NA, nrep)
  for (reprandom in 1:nrep)
  {
    prov[[reprandom]] <- dbfsc(x=dtrain, y=NULL, K=dfvalores[ih, 2], gamma=dfvalores[ih, 1],  
                               diss=TRUE, eps=1e-8, random=TRUE, seed=haziak[reprandom])
    objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
  }
  aux <- which.min(objetivevalues)
  out <- prov[[aux]]
  a <- out$a
  
  results[ih, (dimhyper+1):(dimhyper + dfvalores[ih, 2] )] <- out$a
  results[ih, "G"] <- min(out$objG)
  
  # Permutation distances
  Gb <- rep(NA, B)
  Deltab <- rep(NA, B)
  d_vec <- as.dist(dtrain)
  N <- n*(n-1)/2
  set.seed(258)
  haziak <- sample(1:10000, B, replace=FALSE)
  for (b in 1:B)
  {
    set.seed(haziak[b])
    permutation <- sample(1:N, N, replace=FALSE)
    db_vec <- d_vec[permutation]
    attr(db_vec, "class") <- "dist"
    attr(db_vec, "Size") <- n
    db <- as.matrix(db_vec)
    aux <- dbfsc(x=db, y=NULL, a= a, K=dfvalores[ih, 2], gamma=dfvalores[ih, 1],  diss=TRUE, eps=1e-8, random=FALSE)
    Gb[b] <- min(aux$objG)
  }
  results[ih, "meanLogGb"] <- mean(logtol(Gb, eps=1e-8))
}

# The preliminary results might be saved in case they will be analysed later
save(results, file="PreliminaryResults.RData") 
#--------------------------------------------------

# To select hyper-parameters for each K
results$Gap <-  results$meanLogGb - logtol(results$G, eps=1e-8) 

gd <- results %>% group_by(K) %>% summarise(gamma = gamma[which.max(Gap)])
nR <- dim(results)[1]
gd$index <- NA
for (k in valoresK)
{
  i <- k-1
  gd$index[i] <- (1:nR)[(results$K == k) & ifelse(abs(results$gamma - gd$gamma[i])<=1e-9, TRUE, FALSE)]
}
gd

# Compute Fuzzy silhouette
gd$silF <- NA
for (k in valoresK)
{
  i <- k-1
  fila <- gd$index[i]
  gamma.chosen <- results$gamma[fila]
  K.chosen <- k
  a <- unlist(results[fila, 3:(2+K.chosen)])
  out <- dbfsc(x=dtrain,  K= K.chosen, gamma=gamma.chosen, diss=TRUE, eps=1e-8, a=a)
  gd$silF[i] <- SIL.F(dtrain, out$u, distance=TRUE)
}
# Note: warning note comes from SIL.F

gd$silF

#---------------------------------------------
# Partition in K=2 clusters
  k <- 2
  i <- k-1
  selectedrow <- gd$index[i]
  gamma.chosen <- results$gamma[selectedrow]
  K.chosen <- k
  a <- unlist(results[selectedrow, 3:(2+K.chosen)])
  out <- dbfsc(x=dtrain,  K= K.chosen, gamma=gamma.chosen, diss=TRUE, eps=1e-8, a=a)
  pred <- apply(out$u, 1, which.max)
  table(y, pred)
  
  # Partition in K=4 clusters
  k <- 4
  i <- k-1
  selectedrow <- gd$index[i]
  gamma.chosen <- results$gamma[selectedrow]
  K.chosen <- k
  a <- unlist(results[selectedrow, 3:(2+K.chosen)])
  out <- dbfsc(x=dtrain,  K= K.chosen, gamma=gamma.chosen, diss=TRUE, eps=1e-8, a=a)
  pred <- apply(out$u, 1, which.max)
  table(y, pred)
  cbind(a, out$u[a,]) # Prototypes and their membership values
  deltaF(out$u, dtrain) # Fuzzy Delta values
  #------------------------------------------------
  

