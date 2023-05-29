#---------------------------------------------------------------------------
# Fuzzy classification with distance-based depth prototypes: High-
# dimensional Unsupervised and/or Supervised problems
#---------------------------------------------------------------------------

# Code
# Main functions

dbfsc <- function(x, y=NULL, K=2, a=NULL, gamma=0.05, alpha=1, maxit=15, delta=0.01, diss=FALSE, seed=13, eps=1e-10, 
                  random=TRUE, ndim=2)
{
  # Input:
  #       x: data set
  #       y: labels of individuals (optional), numbers 1 to M
  #       K: number of clusters
  #       a: vector of length K with initial prototypes. If NULL, they are set at random
  #       maxit: maximum of iterations
  #       train: if y is given, train establishes which elements or what percentage of individuals conform the train set
  # Output:
  #       u: memberships of individuals to clusters (nxK)
  #       a: final deepest individuals
  #       atraj: records of deepest individuals
  x <- as.matrix(x)
  if (!is.null(a))
  {
    if(any(a > dim(x)[1]))
    {
      stop("elements in a must be between 1 and ", dim(x)[1])
    }
  }
  if (is.null(y))
  {
    out <- subdbfsc1(x=x, K=K, a=a, gamma=gamma, maxit=maxit, diss=diss, seed=seed, eps=eps, random=random)
  }else
  {
    out <- subdbfsc2(x=x, y=y, K=K, a=a, gamma=gamma, alpha=alpha, maxit=maxit, delta=delta, diss=diss, seed=seed, eps=eps, random=random, ndim=ndim)
  }
  return(out)
}


subdbfsc1 <- function(x, K, a, gamma, maxit, diss, seed, eps, random)
{
  # Input:
  #       x: data set or distance matrix
  #       y: labels of individuals (optional)
  #       K: number of groups
  #       a: vector of length K with initial prototypes
  #       maxit: maximum of iterations
  # Output: NON SUPERVISED VERSION
  #       u: memberships of individuals to groups (nxK)
  #       a: final deepest individuals
  #       atraj: records of deepest individuals
  if(diss)
  {
    dij <- as.matrix(x)
  }else{
    dij <- as.matrix(dist(x))
  }
  
  n <- dim(x)[1]
  if (is.null(a))
  {
    if (random)
    {
      set.seed(seed)
      a <- sample(1:n, K, replace=FALSE)
    }else
      stop("For clustering set elements a")
  }
  if (length(a)!= K)
  {
    stop("There are not K=", K, " elements in a" )
  }
  
  atraj <- matrix(a, nrow=1)
  keep <- TRUE
  count <- 0
  objG <- NULL
  objF <- NULL
  while(keep)
  {
    count <- count + 1
    ddeep <- dij[, a]#matrix(0, n, K)
    ddeep <- ddeep^2
    
    u <- update.u(ddeep, gamma=gamma, eps=eps)
    vik <- apply(u, 2, fpartialPhi, d=dij)
    # 
    a <- apply(vik, 2, update.a, n=n)
    errep <- check.a(a)
    ikusi <- 0
    while(!is.null(errep) & (ikusi <= K))
    {
      aold <- a
      for (f in 1:dim(errep)[1])
      {
        changeind <- which(aold == errep[f, 1]) 
        a[changeind] <- update.a(vik[, changeind[1]], n=n, amount=errep[f, 2])
      }
      errep <- check.a(a)
      ikusi <- ikusi +1
    }
    names(a) <- paste("a", 1:K, sep="")
    atraj <- rbind(atraj, a)
    row.names(atraj) <- NULL
    if ((count > maxit)|(sum(a - atraj[count,])==0))
    {
      keep <- FALSE
    }
    objG <- c(objG, G(ddeep, u, gamma, tol=eps))
    objF <- c(objF, sum(apply(u, 2, Fk, d=dij)))
  } # while
  return(list(u=u, a=a, atraj=atraj, count=count, gamma=gamma, objG=objG, objF=objF, diss=diss, seed=seed, eps=eps))
}


subdbfsc2 <- function(x, y, K, a, gamma, alpha, maxit, delta, seed, diss, eps, random, ndim=ndim)
{
  # Input:
  #       x: data set
  #       y: labels of individuals, numbers 1 to M
  #       K: number of prototypes
  #       a: vector of length K with initial prototypes. If NULL, they are 
  #          assigned by deepest individuals in classes. When K>M, the K-M left
  #          deep individuals are assigned according to the variability in classes. 
  #       maxit: maximum of iterations
  # Output: SUPERVISED VERSION
  #       u: memberships of individuals to groups (nxK)
  #       a: final deepest individuals
  #       atraj: records of deepest individuals
  
  if(diss)
  {
    dij <- as.matrix(x)
  }else{
    dij <- as.matrix(dist(x))
  }
  
  ynum <- as.numeric(y)
  
  
  n <- dim(dij)[1]
  ny <- tabulate(ynum)
  M <- max(ynum)
  
  indicator <- t(sapply(ynum, function(x){x==1:M}))
  if (is.null(a))
  {
    if (random)
    {
      aux <- initialize2.az(n=n, ynum=ynum, K=K, delta=delta, seed=seed)
    }else
    {
      aux <- initialize3.az(ynum=ynum, K=K, delta=delta, dtrain=dij, ndim=ndim)
    }
    a <- aux$a
    z <- aux$z
  }else
  {
    z <- initialize.z(n=n, ynum=ynum, delta=delta, a=a)
    names(a) <- paste("a", 1:K, sep="")
  }
  
  atraj <- a #matrix(a, nrow=1)
  keep <- TRUE
  count <- 0
  objG <- NULL
  objF <- NULL
  
  while(keep)
  {
    count <- count + 1
    ddeep <- dij[, a]
    ddeep <- ddeep^2
    ddeep_ikusi <- ddeep
    
    aux <- t(apply(indicator, 1, function(x){z[x, ]}))
    aux <- apply(aux, 2, puttol, tol=eps)
    ddeep <- ddeep + alpha*(-log(aux))
    
    u <- update.u(ddeep, gamma=gamma, eps=eps)
    vik <- apply(u, 2, fpartialPhi, d=dij)
    
    # 
    a <- apply(vik, 2, update.a, n=n)
    errep <- check.a(a)
    if(!is.null(errep))
    {
      aold <- a
      for (f in 1:dim(errep)[1])
      {
        changeind <- which(aold == errep[f, 1]) 
        a[changeind] <- update.a(vik[, changeind[1]], n=n, amount=errep[f, 2])
      }
      
    }
    names(a) <- paste("a", 1:K, sep="")
    atraj <- rbind(atraj, a)
    row.names(atraj) <- NULL
    for (k in 1:K)
    {
      z[, k] <- apply(u[, k]*indicator, 2, sum)/sum(u[, k])
    }
    objG <- c(objG, G(ddeep, u, gamma, tol=eps))
    objF <- c(objF, sum(apply(u, 2, Fk, d=dij)))
    
    if ((count > maxit)|(sum(a - atraj[count,])==0))
    {
      keep <- FALSE
    }
  } # while
  
  
  return(list(u=u, a=a, atraj=atraj, count=count, labelprot=z, 
              xtrain=x,  gamma=gamma, alpha=alpha, objG=objG, objF=objF, diss=diss, seed=seed, eps=eps, ddeep2=ddeep_ikusi, ddeep2Pl =ddeep))
}

# Auxiliary functions
vgF <- function(uk, d)
{
  # Fuzzy geometric variability of Ck 
  # Input:
  #    uk: membership vector with respect to Ck 
  #    d: pairwise distances
  d <- d^2
  d <- as.matrix(d)
  uk <- unlist(uk)
  uk <- matrix(uk, ncol=1)
  vgeoF <- t(uk)%*%d%*%uk
  nk <- sum(uk)
  vgeoF <- vgeoF/(2*nk^2) 
}

deltaF_kl <- function(uk, ul, d)
{
  # Fuzzy squared distance between Ck  and Cl
  # Input:
  #    uk: membership vector with respect to Ck
  #    ul: membership vector with respect to Cl
  #    d: pairwise distances
  uk <- matrix(uk, ncol=1)
  ul <- matrix(ul, ncol=1)
  vFk <- vgF(uk, d)
  vFl <- vgF(ul, d)
  d <- d^2
  d <- as.matrix(d)
  aux <- t(uk)%*%d%*%ul
  nk <- sum(uk)
  nl <- sum(ul)
  delta <- aux/(nk*nl) - vFk - vFl 
  return(delta)
}


deltaF <- function(u, d)
{
  # Fuzzy distance between Ck  and Cl, k, l=1, ..., K
  # Input:
  #     u: matrix of memberships
  #     d: pairwise distances
  # Output:
  #     delta: Fuzzy distance between Ck and Cl
  K <- dim(u)[2]
  protnames <- colnames(u)
  delta <- matrix(0, K, K, dimnames=list(protnames, protnames))
  for (k in 1:(K-1))
  {
    for (l in (k+1):K)
    {
      uk <- u[, k]
      ul <- u[, l]
      delta[k, l] <- deltaF_kl(uk, ul, d=d)
      delta[l, k] <- delta[k, l]
    }
  }
  return(delta)
}

Fk <- function(uk, d)
{
  # Auxiliary function to build the objective function:
  # Fuzzy geometric variability * sum(u_ik) 
  # Input:
  #    uk: membership vector with respect to Ck
  #     d: pairwise distances
  d <- d^2
  d <- as.matrix(d)
  uk <- matrix(uk, ncol=1)
  Fk <- t(uk)%*%d%*%uk
  Fk <- Fk/(2*sum(uk)) 
}


fpartialPhi <- function(uk, d)
{
  # Input: 
  #      uk: membership to the k-th group
  #      d: pairwise distances between individuals
  # Output: first terms of the adapted phi values phi^2(i, C_k) (sum_j u_jk d^2(i,l)/sum_j u_jk)
  #         proportional to sum_j u_jk d^2(i,l) within class u_jk
  d <- as.matrix(d)
  d <- d^2
  aux <- apply(d, 1, function(x){sum(x*uk)})
  return(aux)
}


G <- function(dprot, u, gamma, tol)
{
  # Global objective function
  # Input:
  #     dprot: distances (squared) to each prototype + alpha*LOSS()
  #     u: matrix of memberships
  #     gamma: hyperparameter in objective function
  #     tol: tolerance for small values
  u <- puttol(u, tol)
  value <- sum(u*dprot) + gamma*sum(u*log(u))
  return(value)
}


initialize2.az <- function(n, ynum, K, delta, seed)
{
  # Initialization of prototypes at random and label-prototypes
  # Input:
  #    n: number of individuals
  #    ynum: labels of units, numerically coded
  #    delta: fuzzy squared distance between C1, ..., CK
  #    seed: integer to set the seed
  # Output:
  #    a: prototypes
  #    z: matrix of label-prototypes
  
  M <- max(ynum)
  z <- matrix(delta, M, K)
 
  # Initialize prototypes
  set.seed(seed)
  a <- sample(1:n, K, replace=FALSE)
  names(a) <- paste("a", 1:K, sep="")

  for (k in 1:K)
  {
    cl <- ynum[a[k]]
    z[cl, k] <- 1
  }
  
  z <- apply(z, 2, function(x){x/sum(x)})
  return(list(a=a, z=z))
}


initialize.z <- function(n, ynum, delta, a)
{
  # Initialization of prototypes by the user and label-prototypes
  # Input:
  #    n: number of individuals
  #    ynum: labels of units, numerically coded
  #    delta: fuzzy squared distance between C1, ..., CK
  #    a: prototypes initialized by user
  # Output:
  #    z: matrix of label-prototypes
  K <- length(a) 
  M <- max(ynum)
  z <- matrix(delta, M, K)
  
  for (k in 1:K)
  {
    cl <- ynum[a[k]]
    z[cl, k] <- 1
  }
  
  z <- apply(z, 2, function(x){x/sum(x)})
  return(z)
}

update.u <- function(dd, gamma, eps)
{
  # Updating of membership matrix
  # Input:
  #     dd: distances to each of the K deepest units (n x K)
  #     gamma: hyperparameter
  #     eps: eps
  u <- exp(-dd/gamma)
  if(nrow(dd)>1)
  {
    pond <- apply(u, 1, sum)
  }else
  {
    pond <- sum(u)
  }
  if(any(pond==0))
  {
    warning("Value of gamma is too small, is causing problems in membership calculations. Set to eps")
    pond[pond==0] <- eps
  }
  u <- u/pond
  return(u)
}


update.a <- function(Iinvk, n=n, amount=1)
{
  # Updating prototypes
  # Input:
  #       Ik: Depth index of group k
  #       amount: number of deepest that are choosen
  # Output: index of the deepest indiviual
  o <- order(Iinvk, decreasing = FALSE)
  a <- (1:n)[o[1:amount]]
  return(a)
}

puttol <- function(v, tol)
{
  # Auxiliary function to to avoid membership values below tolerance
  selec <- (v < tol)
  v[selec] <- tol
  selec <- (v > 1- tol)
  v[selec] <- 1 - tol
  return(v)
}



settrain <- function(N, train, y)
{
  # Input: 
       # N: total individuals in data set
       # train
  # Output: returns indexes of individuals in the train set
  if(length(train)==1)
  {
    if ((train >0)&(train <= 1))
    {
      set.seed(seed)
      selected <- caret::createDataPartition(y, p = train, 
                                             list = FALSE, 
                                             times = 1)
    }
  }else
  {
      if (any(!train%in%(1:N) ))
      {
        stop("Train elements must individual in data") 
      }else
      {
        selected <- train
      }
  }
  
  return(selected)
}

check.a <- function(a)
{
  # Checks whether prototypes are repeated. If yes, the repeated prototype and the amount it is repeated is returned
  kontatu <- table(a)
  cond <- kontatu > 1
  if(sum(cond)>=1)
  {
    marks <- kontatu[cond]
    repeated <- data.frame(el=as.numeric(names(marks)), el.n=c(marks))
    row.names(repeated) <- NULL
  }else
  {
    repeated <- NULL
  }
  return(repeated)
}




# Para introducir la version supervisada tipo L1
subdbfsc2.L1 <- function(x, y, K, a, gamma, alpha, maxit, delta, seed, diss, eps, random)
{
  # Input:
  #       x: data set
  #       y: labels of individuals, numbers 1 to M
  #       K: number of prototypes
  #       a: vector of length K with initial prototypes. If NULL, they are 
  #          assigned by deepest individuals in classes. When K>M, the K-M left
  #          deep individuals are assigned according to the variability in classes. 
  #       maxit: maximum of iterations
  # Output: SUPERVISED vERSION
  #       u: memberships of individuals to groups (nxK)
  #       a: final deepest individuals
  #       atraj: records of deepest individuals
  
  if(diss)
  {
    dij <- as.matrix(x)
  }else{
    dij <- as.matrix(dist(x))
  }
  
  ynum <- as.numeric(y)
  
  
  n <- dim(dij)[1]
  ny <- tabulate(ynum)
  M <- max(ynum)
  
  indicator <- t(sapply(ynum, function(x){x==1:M}))
  if (is.null(a))
  {
    if (random)
    {
      aux <- initialize2.az(n=n, ynum=ynum, K=K, delta=delta, seed=seed)
    }else
    {
      aux <- initialize.az(d=dij, ynum=ynum, K=K, delta=delta, seed=seed)
    }
    a <- aux$a
    z <- aux$z
  }else
  {
    z <- initialize.z(n=n, ynum=ynum, delta=delta, a=a)
    names(a) <- paste("a", 1:K, sep="")
  }
  

  atraj <- a #matrix(a, nrow=1)
  #colnames(atraj) <- names(a)
  keep <- TRUE
  count <- 0
  objG <- NULL
  objF <- NULL
  
  while(keep)
  {
    count <- count + 1
    ddeep <- dij[, a]#matrix(0, n, K)
   # ddeep <- ddeep^2
    aux <- t(apply(indicator, 1, function(x){z[x, ]}))
    aux <- apply(aux, 2, puttol, tol=eps)
    ddeep <- ddeep + alpha*(-log(aux))
    
    u <- update.u(ddeep, gamma=gamma, eps=eps)
    #vg <- apply(u, 2, vgF.L1, d=dij)
    vik <- apply(u, 2, fpartialPhi.L1, d=dij)
    # 
    a <- apply(vik, 2, update.a, n=n)
    names(a) <- paste("a", 1:K, sep="")
    atraj <- rbind(atraj, a)
    row.names(atraj) <- NULL
    for (k in 1:K)
    {
      z[, k] <- apply(u[, k]*indicator, 2, sum)/sum(u[, k])
    }
    objG <- c(objG, G(ddeep, u, gamma, tol=eps))
    objF <- c(objF, sum(apply(u, 2, Fk.L1, d=dij)))
    
    if ((count > maxit)|(sum(a - atraj[count,])==0))
    {
      keep <- FALSE
    }
  } # while
  
  
  return(list(u=u, a=a, atraj=atraj, count=count, labelprot=z, 
              xtrain=x,  gamma=gamma, alpha=alpha, objG=objG, objF=objF, diss=diss))
}


calcpred2 <- function(dbfscobj, ddeep)
{
 
  n <- dim(xtest)[1]
  gamma <- dbfscobj$gamma
 # a <- dbfscobj$a
  z <- dbfscobj$labelprot
  K <- dim(z)[2]
  M <- dim(z)[1]
  ddeep <- ddeep^2
  u <- update.u(ddeep, gamma=gamma, eps=eps)
  pred <- apply(u%*%t(z), 1, which.max)
  return(pred)
}

 calcpred <- function(dbfscobj, xtest, diss=FALSE)
 {
   x <- dbfscobj$xtrain
   xtest <- as.matrix(xtest)
   n <- dim(xtest)[1]
   gamma <- dbfscobj$gamma
   a <- dbfscobj$a
   z <- dbfscobj$labelprot
   K <- dim(z)[2]
   M <- dim(z)[1]
   ddeep <- matrix(0, n, K)
   for (k in 1:K)
   {
     ddeep[, k] <- apply(xtest, 1, fd, xa=x[a[k], ], diss=diss) 
   }
   u <- update.u(ddeep, gamma=gamma, eps=eps)
   pred <- apply(u%*%t(z), 1, which.max)
   return(pred)
 }


