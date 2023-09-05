rm(list=ls())
#------ Example to execute FC-DF with search for hyperparameters-------
source("functions.R")
# 
# Generate simulated data
source("generatedata.R")
out <- generate(n=500, seed=5231)
x <- out[[1]]
y <- out[[2]]
n <- dim(x)[1]
rm(out)

# Set grid for hyperparameters,
values <- c(2^seq(-1, 4, by=0.5))
valuesK <- 4#3:8
dfval <- expand.grid(gamma = values , alpha= values, K=valuesK)
dim(dfval) 
dimhyper <- 3

# Repetition because of random initialization  and seed setting
nrep <- 5
set.seed(3201)
haziak <- sample(1:10000, nrep, replace=FALSE)

# -----------FC-DF ---------------------------
# Stratified, 10 fold cross validation for hyper-parameter selection
folds <- 10
init <- matrix(NA, nrow=dim(dfval)[1], ncol=n + max(valuesK)) # 
# Provisional results will be gathered in "results" so that the best choice of hyper-parameters can be choosen
results <- data.frame(dfval, init)
names(results) <- c("gamma", "alpha", "K", paste("a", 1:max(valuesK), sep=""), paste("pred", 1:n, sep=""))

d <- as.matrix(dist(x))
# Seed set inside the function
splits <- splitTools::create_folds(y, k=folds, type="stratified", seed=2562)

for (ih in 1:dim(dfval)[1])
{
  cat("\nih = ", ih, "out of ", dim(dfval)[1])
  for (i in 1:folds)
  {
    # No need to recalcutate the distance matrix
    indextrain <- splits[[i]]
    indextest <- (1:n)[-indextrain]
    dtrain <- d[indextrain, indextrain]
    ytrain <- y[indextrain]
    deval <- d[-indextrain, indextrain]
    set.seed(364)
    haziak <- sample(1:10000, nrep, replace=FALSE)
    aux <- rep(NA, nrep)
    prov <- vector("list", nrep)
    objetivevalues <- rep(NA, nrep)
    for (reprandom in 1:nrep)
    {
      prov[[reprandom]] <- dbfsc(x=dtrain, y=ytrain, K=dfval[ih, 3], gamma=dfval[ih, 1], alpha=dfval[ih, 2], 
                                 diss=TRUE, eps=1e-8, random=TRUE, seed=haziak[reprandom])
      objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
      aux <- which.min(objetivevalues)
      out <- prov[[aux]]
    }
  
    a <- out$a # 
    z <- out$labelprot
    dt <- matrix(deval[, a], nrow=length(indextest)) # ntest x K matrix
    dt <- dt^2
    u <- update.u(dt, gamma=dfval[ih, 1])
    pred <- apply(u%*%t(z), 1, which.max)
    results[ih, (dimhyper+1):(dimhyper + dfval[ih, 3] )] <- a
    results[ih, dimhyper+ max(valuesK) +indextest] <- pred
  }
}

propacc <- function(x, ytest)
{
  aux <- ifelse(x==ytest, 1, 0)
  return(mean(aux))
}
red <- as.matrix(results[, -(1:(dimhyper+ max(valuesK)))])
results$accuracy <- apply(red, 1, propacc, ytest=y)

# Choose best hyper-parameters
selecthyp <- results$accuracy==max(results$accuracy)
sum(selecthyp) # several combinations lead best accuracies
results[selecthyp, c("gamma", "alpha", "accuracy")]
# We keep the first one
best <- which.max(results$accuracy) 

# Fit the model with the selected hyper-parameters
selectedgamma <- results$gamma[best]
selectedalpha <- results$alpha[best]
selectedK <- results$K[best]

d <- as.matrix(dist(x))
aux <- rep(NA, nrep)
prov <- vector("list", nrep)
objetivevalues <- rep(NA, nrep)
for (reprandom in 1:nrep)
{
  prov[[reprandom]] <- dbfsc(x=d, y=y, K=selectedK, gamma=selectedgamma, alpha=selectedalpha, 
                             diss=TRUE, eps=1e-8, random=TRUE, seed=haziak[reprandom])
  objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
  aux <- which.min(objetivevalues)
  out <- prov[[aux]]
}
# Final model 
Finalmodel <- out



# ------------Prediction-------------------
# Generate test set
test <- generate(n=500, seed=2575)
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