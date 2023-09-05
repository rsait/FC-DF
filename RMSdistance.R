#------------------ RELATED METRIC SCALING ------------------
#
# In the following, there is the code to calculate the Related Metric
# Scaling distance given d1, ..., dk distance matrices.
# The code contains 3 auxiliary functions (doublecentered, sqrtB, vg2unit) and one principal
# function (RMS.dist).
#
# Details:
# - In the procedure, all distances are scaled so that their geometric variability is equal to 1.
# - Eigen values of the doubly centered matrix are set to 0 if their numerical value 
# is in (-eps, eps). By default, eps is set as eps=10^-10.
# - When a considered distance is not Euclidean, that means that the lowest eigen value of its
# doubly centered matrix is negative (l < -eps), then a transformation is made so that 
# the trasnformed distance matrix is Euclidean. The transformed distances are d* = sqrt(d^2 - 2*l), 
# where d^2 stands for the squared values of the distance matrix entries and l stands for the lowest eigen value.
#
# Example of use: 
# Consider distance matrices d1 and d2 calculated on the same n individual
# To calculate the RMS distance:
# d <- list(d1, d2)
# d.rms <- RMS.dist(d)
#----------------------------------------------------------------------------

# Auxiliar functions
doublecentered <- function(d)
    {
        #----- Calculates B, to doubly center distance matrix
        # Input:
        #     d: distance matrix
        # Output:
        #     B:  matrix to doubly center matrix
        #-----------------------------------------
        d <- as.matrix(d)
        n <- dim(d)[1]
        unos <- matrix(1, nrow=n, ncol=1)
        I <- diag(rep(1, n))
        H <- I - 1/n*unos%*%t(unos)
        d2 <- d^2
        B <- -1/2*H%*%d2%*%H
        return(B)
    }
sqrtB <- function(B, eps=10^-10)
    {
        #----- Calculates the square root of a semidefinite positive matrix 
        #  
        # Input:
        #     B: a  matrix
        #     eps: accuracy limit for numerical 0
        # Output:
        #     B.sqrt:  B^{1/2}
        #     l: minumum eigen value
        #-----------------------------------------  
      vp <- eigen(B)
      lambda <- vp$values
      s <- (lambda > -eps) & (lambda < eps)
      lambda[s] <- 0
      L <- diag(lambda)
      U <- vp$vectors
      L.sqrt <- sqrt(abs(L))
      B.sqrt <- U%*%L.sqrt%*%t(U) #Ok
      lmin <- min(lambda)
      out <- list(B.sqrt=B.sqrt, l=lmin)
      return(out)
    }

vg2unit <- function(d)
    {
        #----- Scale a distance matrix to have unit geometric variability
        # Input:
        #     d: original distance matrix
        # Output:
        #     d: scaled distance matrix so that vgeo(d)=1
        #-----------------------------------------                
        d <- as.matrix(d)
        n <- dim(d)[1]
        vg <- sum(d^2)/(2*n*n)
        d <- d/sqrt(vg)
    }


# Principal function
RMS.dist <- function(d, eps=10^-10)
{
        #----- To calculate Related MEtric Scaling
        # Input:
        #     d: a list with k>=2 distance matrices 
        #     eps: accuracy for numerical 0. Values in (-eps, eps) are considered as 0
        # Output:
        #     drel: RMS distance
        #-----------------------------------------                    

    k <- length(d)
    n <- dim(d[[1]])[1]
    du <- lapply(d, vg2unit)
    Bs <- lapply(du, doublecentered)
    aux <- lapply(Bs, sqrtB, eps=eps)
    Bsqrts <- list()
    for (i in 1:k)
    {
      if(aux[[i]]$l < -eps)
      {
        cat("Distance matrix ", i, " is not Euclidean. It has been transformed\n")
        d[[i]] <- as.matrix(d[[i]])
        dt <- sqrt(d[[i]]^2 - 2*aux[[i]]$l)
        diag(dt) <- 0
        dt <- vg2unit(dt)
        Bs[[i]] <- doublecentered(dt)
        Bsqrts[[i]] <- sqrtB(Bs[[i]])$B.sqrt
      }else
      {
        Bsqrts[[i]] <- aux[[i]]$B.sqrt
      }
    }
    Brel <- matrix(0, n, n)
    aux <- matrix(0, n, n)
    for (i in 1:k)
        {
           Brel <- Brel + Bs[[i]]
           for (j in 1:k)
               {
                   if (j != i)
                       {
                           aux <- aux + Bsqrts[[i]]%*%Bsqrts[[j]]
                       }
               }
        }
    Brel <- Brel - 1/k*aux
    b <- matrix(diag(Brel), ncol=1)
    unos <- matrix(1, ncol=1, nrow=n)
    drel.2 <- b%*%t(unos) + unos%*%t(b) - 2*Brel
    drel <- sqrt(drel.2)
    return(drel)
}


