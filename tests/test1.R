library(msm)  # for rtnorm
library(mvtnorm)
library(Matrix)
library(mgcv)
library(clv)
library(tmvtnorm)
library(MCMCpack)
library(devtools)
library(MNP) 

source('~/src/simCrop/MNPSARTools.R')



########################## FUNCTIONS ##########################################

# get category from a latent variable using probit rules
getLatentState <- function( z ) {
  return(
    apply(z,1,function(x){ if( sum(x > 0) == 0 ){return(1)} else{return(1+which.max(x))} } )  
  )
}

# this is a function that takes a list of indexes and creates
# a block diagonal matrix out of them where each block is of
# size 1 by x.max
#
# | -- j=1 -- | -- j=2 -- | -- ... -- | -- j=J -- |  as diag
# | p1,...,pP | 

designElement <- function(x,x.max) {
  # create a blank block 
  blankBlock <- matrix(0,nrow=length(x),ncol=length(x)*x.max)
  for( i in 1:length(x))  blankBlock[ i, (i-1)*x.max + x[i] ] <- 1
  
  return( blankBlock )
}


# it is assumed that the matrix x has columns that are increasing with respect to year
# the result is a list of prior years from left to right t-1 t-2 ... t-k
createDesignMatrix <- function(x, priorStates) {

  n <- NROW(x) 
  m <- length(x)/n
  J <- max(x) 
  P <- (J+1)*(J+1)
  
  X.new <- c()
  for( i in 1:(m - priorStates) ) {
    
    start.index <- (i-1)*n + 1
    stop.index  <- (i+priorStates)*n
    
    X.new <- rbind(
      X.new,
      cbind(
        i,
        matrix(x[start.index:stop.index],nrow=n)[,(1+priorStates):1]
      )
    )
  }

  colnames(X.new) <- c( 'i', sprintf("%d", -1 * 1:(priorStates + 1)) ) 

  # convert prior states into sequences
#  Y10 <- (Y1*3 + Y0)-3
#
#    X.design0 <- c()
#    for( i in 1:length(Y10)) {
#      de <- designElement( rep( Y10[i] ,J), P) 
#      X.design0 <- rbind( X.design0,de  )
#    }
  
  # get prior states 
#  priorStateDesignMatrix(y) 

  return(X.new)
}



######################################
# hyper parameters 
######################################

#Size of Transition matrix:
#  prior states
#   nonAg <- nonAg
#   nonAg <- corn
#   nonAg <- soybeans
#   corn <- nonAg
#   corn <- corn
#   corn <- soybeans
#   soybeans <- nonAg
#   soybeans <- corn
#   soybeans <- soybeans
J <- 2
P <- (J+1)*(J+1) 
# Sigma for transition
BetaStar  <- matrix( 0 , nrow= n*J*P) 
Sigma0  <- diag(P)   # Sigma for the prior states
SigmaJ0 <- diag(J)   # Sigma for the latent variable dim in Beta prior

# hyper parameters for Sigma.inv
v <-  10 
S <- diag( J )


######## design Matrix ##############


######## initial parameters ##############
Beta      <- rep(0, n*J*P)
Sigma.inv <- diag(J)
muZ       <- X.design %*% Beta

######## constants ########################
m <- 10 # thinning for multivariate normal
J <- 2  # number of alternatives



##############################################################

# Imai 2005
set.seed(400)


# initial value for the random variable alpha^2
v <- 3 
pMinusOne <- ncol(X.design)
alpha0 <- 1

Sigma <- diag(2)
muZ <- X.design %*% c(0,0) 

# get initial z
z <- rnorm(n*J)

SigInv <- Sigma.inv 

#### start iter here ####

for( iter in 1:10 ) {

  ############### ALPHA for Z ######################### 
    
  ## get scale parameter 
  chiScale <- alpha0^2 * sum(diag( S %*% SigInv   ))
  
  # shape, rate
  alpha.z <- chiScale / rchisq(1, df=(v*(pMinusOne)) ) 
  print(sprintf("alpha.z = %f", alpha.z)) 
  
  ############### Z ########################## 

  zMaxArray <- z
  XBeta <- X.design %*% Beta 
  
  # guess for cmean and cvar
  #cvars <-c()
  #cmeans <-c()

  for( i in 1:n) {
      
    #cmean <- z.central[(i*J-1):(i*J)][-p]
   
    for( p in 1:J) {

      cvar <- 1/SigInv[p,p]
     
      cmean <- XBeta[(i*J-1):(i*J)][p] -  (  z[(i*J-1):(i*J)] - XBeta[(i*J-1):(i*J)])[-p] * SigInv[p,-p] * cvar
  
      zMax = max( max( z[(i*J-1):(i*J)][-p] ),0)
                      
 
#      for( p2 in 1:J) {
#        cmean <- cmean - cmean[p2] * cvar;
#      }

      zMaxArray[i*J-2 + p] = zMax
      if( Y[i]  == p +1 ) {
        z[i*J -2 + p] <- TruncNorm( zMax, cmean + 1000 * sqrt(cvar), cmean, cvar, 0 ) 
      } else {
        z[i*J -2 + p] <- TruncNorm( cmean - 1000 * sqrt(cvar), zMax, cmean, cvar, 0 ) 
      }
 
#      cvars <- c(cvars,cvar)
#      cmeans <- c(cmeans,cmean)

    }
  }
  
  zSquiggle <- z * sqrt(alpha.z) 
  
    
  ############### ALPHA for Beta ######################### 
  XStar = cbind( X.design %*% chol(SigInv), c(chol(SigInv) %*% matrix(zSquiggle,nrow=2)))
  SS = t(XStar) %*% XStar

  # sweep to get beta hat etc...
  SS.sweep <- SWP(SWP(SS,0),1)
  
  chiScale2 <- SS.sweep[J+1,J+1] + chiScale 
  
  alpha.Beta <- c(chiScale2 / rchisq(1,df=(n + v)*(pMinusOne ) ) )
  print(sprintf("alpha.Beta = %f", alpha.Beta)) 
  
  ############ Beta ######################
  
  beta.var <- SS.sweep[1:J,1:J] * -1 * alpha.Beta
  beta.mu  <- SS.sweep[1:J,J+1]
  
  BetaSquiggle <- RMVN( beta.mu, beta.var ) 
  
  ############ Wishart ######################
  
  epsilon <- zSquiggle - X.design %*% BetaSquiggle
  epsilon <- matrix(epsilon,ncol=2,byrow=TRUE)
  
  R <- matrix(0,nrow=J,ncol=J)
  for(i in 1:n) R <- R +  t(epsilon[i,,drop=FALSE]) %*% epsilon[i,,drop=FALSE]
  
  
  SigInv <- RWISH( solve(R + A), v + n)
  Sigma <- solve(SigInv)
  
  ############ Rescale ######################
  
  Beta <- BetaSquiggle/sqrt(alpha.Beta)
  
  alpha.Sigma <- mean(diag(Sigma))
  Sigma <- Sigma / alpha.Sigma
  SigInv <- SigInv * alpha.Sigma 
  z <- zSquiggle / sqrt(alpha.Sigma)

  z.save <- z
 
} 






#set.seed(400)
###############################################################
## run model
#res <- mnp( Y ~ 1,
#           n.draws = 10,
#           burnin = 0,
#           thin = 0,
#           verbose = TRUE,
#           latent = TRUE
#           )
###############################################################
#

#
# 
#
# need:  det|(I - pW)|^{-J} det|\Sigma_0|^{-n/2}
#  exp( -1/2 * b_i^{T}A^{-1}b_i )
#  A = (I - pW)^{-1}\Sigma_0(I - pW)^{-1}


Z  <- z
x0 <- .1
Years <- 1
Sigma0.inv <- S
lambda.range <- sort(1/range(eigen(W)$values))

mh.lambda <- function(Z,W,x0,iter,burnIn=50, lambda.range) {
    
  n <- nrow(W)
  J <- length(Z) / n
  nJ <- n*J

  results <- 1:iter

  pdfValue <- function(x,W,Sigma0.inv,n,J,nJ) { 
    B <- diag(n) - x * W 
    B.det <- 1/det(B) 
    B.inv <- solve(B)
    B.inv.kronecker <- kronecker(B.inv, diag(J))
        
    A <- solve( B.inv.kronecker %*% 
               kronecker(diag(n),Sigma0.inv) 
               %*% B.inv.kronecker)
     
    B.linear <- 0
    for(year in Years) {
      Zj <- matrix(Z[(nJ*(year-1)):(nJ*year)],ncol=1) 
      B.linear <- B.linear + t(Zj) %*% A %*% Zj
    }
    return( B.det^J *  exp( - 1/2 * B.linear) )
  }

  x.pdf <- pdfValue(x0,W,Sigma0.inv,n,J,nJ) 

  # iterations for metropolis hastings algorithm
  for(i in 1:(iter + burnIn) ) {
    y <- runif(1, min=lambda.range[1],max=lambda.range[2]) 
    y.pdf <- pdfValue(y,W,Sigma0.inv,n,J,nJ) 
    
    u <- runif(1)

    #rho <- min( foo(y) /  foo(x), 1 ) 
    rho <- min( y.pdf/x.pdf, 1 )

    if (u < rho) {
      x <- y
      x.pdf <- y.pdf
    }

    if( i > burnIn) results[i - burnIn] <- x
  }

  return( results )
}






