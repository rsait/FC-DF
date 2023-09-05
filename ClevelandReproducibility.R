# With reproducibility purposes, we publish the code we used with Alizadeh data set
# Tuning for classification
# To understand the use of the available R functions, the reader should revise example1.R and example2.R


library(ICGE)
library(splitTools)
library(MASS)
library(dplyr)
source("functions.R")
source("RMSdistance.R")



dat <- read.table("datos/Cleveland_data.csv", header=TRUE, sep=";", dec=",")
# Omision of NA values
dat <- na.omit(dat)
names(dat) 

# Class variable
y <- factor(dat$Classes)
ycuant <- as.numeric(y)

# Explicative variables
x <- dat[, 1:13] # 
n <- dim(x)[1]
row.names(x) <- 1:n
nombres <- 1:n



#------------------------------------------

# Seed set inside the function
splits <- splitTools::create_folds(y, k=3, type="stratified", seed=202852)
#-----------------------------------------------------------
selectedFold <- 1
indextest <- (1:n)[-splits[[selectedFold]]] # Set test set
indextrain <- splits[[selectedFold]] # Set training set
ntrain <- length(indextrain)
table(y[indextrain])
table(y[indextest])
# Set type of varibles: quantitatives, binaries, qualitatives
quant <- c(1:5, 12)
bin <- 6:7
quali <- c(8:11, 13)
# Compute standard deviations within training set for quantitative variables
SDs <- apply(x[indextrain, quant], 2, sd, na.rm=TRUE) 
SDs


# Gower distance for qualitative variables
qbin <- length(bin) # 2
qquali <- length(quali) # 4
dAtrain <- dgower(x[indextrain, c(bin, quali)], type=list(bin=1:qbin, nom=(qbin+1):(qbin+qquali)))

# Standarize quantitative variables: train and test individuals are standarized according to training values
# Euclidean distance
xT <- t(apply(x[, quant], 1, function(x, sd=SDs){x/SDs}))
dim(xT)
apply(xT, 2, sd)
apply(xT[indextrain, ], 2, sd) # Ok
dTtrain <- dist(xT[indextrain, ])

# Related distance for train set
d <- list(as.matrix(dAtrain), as.matrix(dTtrain))
dm <- RMS.dist(d)

# Related distance for all individuals, train + test
dAcomp <- dgower(x[, c(bin, quali)], type=list(bin=1:qbin, nom=(qbin+1):(qbin+qquali)))
dTcomp <-  dist(xT[ , ]) 
daux <- list(as.matrix(dAcomp), as.matrix(dTcomp))
dcomp <- RMS.dist(daux)
#plot(cmdscale(dcomp), col=y, asp=1 )
#points(cmdscale(dcomp)[indextrain,], bg=y, pch=21)
 
#-----------------------------------------------------------
# # Hyper-parameters search (It takes about 10.6 min)
nrep <- 5
# To save preliminary results
izena <- "PreliminaryCleveland.RData"
valoresgamma <- c(2^seq(-4, 1, by=0.5)) 
valoresalpha <- c(2^seq(-6, 1, by=0.5))
valoresK <- 2:10
dfvalores <- expand.grid(gamma = valoresgamma , alpha= valoresalpha, K=valoresK)
dim(dfvalores) 
dimhyper <- 3

folds <- 10
ptm <- proc.time()
ntrain <- dim(dm)[1]
init <- matrix(NA, nrow=dim(dfvalores)[1], ncol=ntrain) 
results <- data.frame(dfvalores, init)
names(results) <- c("gamma", "alpha", "K", paste("pred", 1:ntrain, sep=""))
splitsINtrain <- splitTools::create_folds(y[indextrain], k=10, type="stratified", seed=62)
dreduced <- dm #  dm contains distance between individuals in training set
yreduced <- y[indextrain]
ycuantreduced <- ycuant[indextrain]

for (ih in 1:dim(dfvalores)[1])
{
  cat("\nih = ", ih, "out of ", dim(dfvalores)[1])
  for (i in 1:folds)
  {
    # It is not necessary to compute again the whole distance matrix
    indextrainreduced <- splitsINtrain[[i]]
    indexvalidation <- (1:ntrain)[-indextrainreduced]
    dtrainreduced <- dreduced[indextrainreduced, indextrainreduced]
    ytrainreduced <- yreduced[indextrainreduced]
    deval <- dreduced[indexvalidation, indextrainreduced]
    set.seed(2164)
    haziak <- sample(1:10000, nrep, replace=FALSE)
    aux <- rep(NA, nrep)
    prov <- vector("list", nrep)
    objetivevalues <- rep(NA, nrep)
    for (reprandom in 1:nrep)
    {
      prov[[reprandom]] <- dbfsc(x=dtrainreduced, y=ytrainreduced, K=dfvalores[ih, 3], gamma=dfvalores[ih, 1], alpha=dfvalores[ih, 2], 
                                 diss=TRUE, eps=1e-8, 
                                 random=TRUE, seed=haziak[reprandom])
      objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
    }
    aux <- which.min(objetivevalues)
    out <- prov[[aux]]
    
    a <- out$a # It is a number between 1 and ntrain
    z <- out$labelprot
    dt <- matrix(deval[, a], nrow=length(indexvalidation)) # ntest x K matrix
    dt <- dt^2
    u <- update.u(dt, gamma=dfvalores[ih, 1])
    pred <- apply(u%*%t(z), 1, which.max)
    results[ih, dimhyper + indexvalidation] <- pred
  }
}
proc.time() - ptm # elapsed 


propacc <- function(x, ytest)
{
  aux <- ifelse(x==ytest, 1, 0)
  return(mean(aux))
}

red <- as.matrix(results[, -(1:dimhyper)])
results$accuracy <- apply(red, 1, propacc, ytest=ycuantreduced)

save(results, file=izena) 

#-----------------------------------------------------------------------------
# Estudy of preliminary results to choose suitable hyper-parameters
# Maximum accuracy value can be achieved by several combination of hyper-parameters. It is shown the first one.
load(izena)

gd <- results %>% 
  group_by(K) %>% 
  summarise(
    gamma = gamma[which.max(accuracy)],
    alpha = alpha[which.max(accuracy)], 
    EvalSetAcc = accuracy[which.max(accuracy)]
  )
gd

##-----------------------------------------------------------------------------

# Evaluation on test set
nrep <- 5 
dtrain <- dm 
ytrain <- y[indextrain]
ycuanttrain <- ycuant[indextrain]
dtest <- dcomp[indextest, indextrain]
objetos <- vector("list", dim(gd)[1])
gd$accuracy <- NA

for (i in 1:dim(gd)[1])
{
  set.seed(12568)
  haziak <- sample(1:10000, nrep, replace=FALSE)
  gammachosen <- gd$gamma[i]
  alphachosen <- gd$alpha[i]
  Kchosen <- gd$K[i]
  prov <- vector("list", nrep)
  objetivevalues <- rep(NA, nrep)
  for (reprandom in 1:nrep)
  {
    prov[[reprandom]] <- dbfsc(x=dtrain, y=ytrain, K=Kchosen, gamma=gammachosen, alpha=alphachosen, 
                               diss=TRUE, eps=1e-8,
                               random=TRUE, seed=haziak[reprandom])
    objetivevalues[reprandom] <- min(prov[[reprandom]]$objG)
  }
  aux <- which.min(objetivevalues)
  out <- prov[[aux]]
  a <- out$a # It is a number between 1 and dim(dtrainreduced), i.e, ntrain
  ycuant_a <- ycuanttrain[a]# Quantitative values of prototypes
  z <- out$labelprot
  dt <- matrix(dtest[, a], nrow=length(indextest)) # ntest x K matrix
  dt <- dt^2
  u <- update.u(dt, gamma=gammachosen)
  p <- u%*%t(z)
  pred <-  apply(p, 1, which.max)
  ytrue <- y[indextest]
  ycuanttrue <- ycuant[indextest]
  acc <- sum(pred==ycuanttrue)/length(ycuanttrue)
  gd$accuracy[i] <- acc
  objetos[[i]] <- list(out=out, u=u, p=p, pred=pred, ytrue=ytrue, ycuanttrue=ycuanttrue, accuracy=acc)
}

gd
#----------------------------------------------------------------------
