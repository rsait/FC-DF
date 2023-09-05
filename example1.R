rm(list=ls())
#------ Simple example to execute FC-DF-------
source("functions.R")
# 
# Generate simulated data
source("generatedata.R")
out <- generate(n=500, seed=3663)
x <- out[[1]]
y <- out[[2]]
rm(out)

# Set hyperparameters,
K <- 5
gamma <- 0.5
alpha <-  2^1.5

# Repetition because of random initialization  and seed setting
nrep <- 5
set.seed(3201)
haziak <- sample(1:10000, nrep, replace=FALSE)

# -----------FC-DF ---------------------------

# Euclidean distance matrix   
d <- as.matrix(dist(x))
aux <- rep(NA, nrep)
prov <- vector("list", nrep)
objetivevalues <- rep(NA, nrep)
for (reprandom in 1:nrep)
{
  prov[[reprandom]] <- dbfsc(x=d, y=y, K=K, gamma=gamma, alpha=alpha, 
                             diss=TRUE, eps=1e-8, random=TRUE, seed=haziak[reprandom])
  objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
  aux <- which.min(objetivevalues)
  out <- prov[[aux]]
}
# Final model 
Finalmodel <- out

# ------------Prediction-------------------
# Generate test set
test <- generate(n=500, seed=325)
xtest <- test[[1]]
ytest <- test[[2]]
rm(test)
a <- Finalmodel$a # 
z <- Finalmodel$labelprot
aux <- rbind(x[a,], xtest)
deval <- as.matrix(dist(aux))
dt <- deval[(length(a)+1):(length(a)+dim(xtest)[1]), 1:length(a)] # ntest x K matrix
dt <- dt^2
u <- update.u(dt, gamma=out$gamma)
# Final prediction
pred <- apply(u%*%t(z), 1, which.max)
# Confusion matrix
table(pred, ytest)   
