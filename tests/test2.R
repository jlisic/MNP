library(msm)  # for rtnorm
library(mvtnorm)
library(Matrix)
library(mgcv)
library(clv)
library(tmvtnorm)
library(MCMCpack)
library(devtools)
library(MNP) 
library(plyr)

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
# the prior states matrix is assumed to be in the order t-1 to t-k.
createDesignMatrix <- function(x, nPriorStates, priorStatesMatrix=NULL, deleteNoMatch=FALSE, autocorrelation=FALSE) {

  n <- NROW(x) 
  m <- length(x)/n
  J <- max(x) - 1 

  X.new <- c()
  for( i in 1:(m - nPriorStates) ) {
    
    start.index <- (i-1)*n + 1
    stop.index  <- (i+nPriorStates)*n
    
    X.new <- rbind(
      X.new,
      cbind(
        i,
        matrix(x[start.index:stop.index],nrow=n)[,(1+nPriorStates):1]
      )
    )
  }

  # at this point we have a matrix with time index, current state, and prior states
  colnames(X.new) <- c( 'i', sprintf("%d", -1 * 0:nPriorStates) ) 

  # if prior state matrix is missing create one
  if( is.null(priorStatesMatrix) ) {
    priorStatesMatrix <- matrix(1:J,J)
    if( nPriorStates > 1 ) {
      for( i in 2:nPriorStates ) 
        priorStatesMatrix <- cbind( kronecker( priorStatesMatrix, matrix(1,J,1) ), rep(1:J, J^(i-1)) )
    }
  }

  # fixes for prior state matrix
  if( is.null( dim(priorStatesMatrix) ) ) priorStatesMatrix <- c(priorStatesMatrix,ncol=1)
  P <- nrow(priorStatesMatrix)
  priorStatesMatrix <- cbind( priorStatesMatrix, 1:P)
  colnames( priorStatesMatrix) <- c(sprintf("%d",-1* 1:nPriorStates), 'p')
  

  # merge data together
  X.new <- join( as.data.frame(X.new), as.data.frame(priorStatesMatrix), by=sprintf("%d",-1*1:nPriorStates) )


  if( is.na( sum(X.new[,'p']) ) ) {
    P <- P + 1
    X.new[ is.na(X.new[,'p']) ,'p'] = P
    cat("Creating prior state for non-matching prior states.\n")
  } 


  # generate design matrix
  X.design <- c()
  for( i in 1:nrow(X.new) ) {
    if( autocorrelation ) { 
      de <- designElement( rep( X.new[i,'p'] ,J), P) 
      de.offset <- matrix(0,ncol=n)
      de.offset[1,ceiling(i/(m-nPriorStates))] <- 1 
      de <- kronecker( de.offset, de )  
    } else {
      de <- designElement( rep( X.new[i,'p'] ,J), P) 
    }
    X.design <- rbind( X.design,de  )
  }

  return(list(map=X.new, design=X.design))
}

# example 
a <- matrix( sample(1:3,size=10,replace=T), ncol=2)
a.design <- createDesignMatrix( a, 1,autocorrelation=T)

# not sure if this is used anymore
#priorStateDesignMatrix <- function(x) {
#  # create design matrix
#  y <- rep(0,max(x))
#  return( t(sapply(x,function(z) {y[z] = 1; return(y)}) ) )
#}



# input X 
# sum_{i=1}^n t(x_i) Sigma.inv x_i
#getXX <- function( X, Sigma.inv ) {
#  
#  X <- matrix(X,nrow=NROW(X)) 
#  J <- sqrt(length(Sigma.inv))
#  n <- nrow(X)/J
#  p <- ncol(X)
#
#  XX <- matrix(0,p,p)
#
#  for( i in 1:n) {
#    x <- X[(J*(i-1) + 1):(J*i),,drop=FALSE]
#    XX <- XX + t(x) %*% Sigma.inv %*% x 
#  }
#
#  return(XX)  
#}
#
#
#getXZ <- function( X, Z, Sigma.inv, outerProduct=FALSE ) {
#
#  X <- matrix(X,nrow=NROW(X)) 
#  Z <- matrix(Z,nrow=NROW(Z)) 
#
#  if( nrow(X) != nrow(Z) ) stop("rows of X do not equal Z")
#
#  J <- sqrt(length(Sigma.inv))
#  n <- nrow(X)/J
#  p1 <- ncol(X)
#  p2 <- ncol(Z)
#
#
#  if( !outerProduct ) {
#  
#    XZ <- matrix(0,p1,p2)
#
#    for( i in 1:n) {
#      x <- X[(J*(i-1) + 1):(J*i),,drop=FALSE]
#      z <- Z[(J*(i-1) + 1):(J*i),,drop=FALSE]
#      XZ <- XZ + t(x) %*% Sigma.inv %*% z 
#    }
#  
#    return(XZ)  
#  }
#    
#  XZ <- matrix(0,J,J)
#
#  for( i in 1:n) {
#    x <- X[(J*(i-1) + 1):(J*i),,drop=FALSE]
#    z <- Z[(J*(i-1) + 1):(J*i),,drop=FALSE]
#    XZ <- XZ + x %*% t(z)
#  }
#
#  return(XZ)  
#}

#############################################################################

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
# Hyperparameter for 
Beta0 <- rep( 0.1, P )
# Sigma for transition
Sigma0  <- diag(P)   # Sigma for the prior states
SigmaJ0 <- diag(J)   # Sigma for the latent variable dim in Beta prior
SigmaE0 <- diag(J)   # Sigma for epsilon

# generate sample size
# note that n is the number of spatial units, not necessarily the number 
n <- 300 
n <- 64

# l is the number of elements, needs to be a sqrt
W <- as.matrix(
       bandSparse(n,n,
       k=list(-1*sqrt(n),-1,1,sqrt(n)))
       )*1

# get range of rho
lambda.range <- sort(1/range(eigen(W)$values))
lambda0 <- 0.20




# for each observation we want the true sigma replicated via the spatial matrix
# Sigma = ( (I - pW)^-1 Sigma (I-pW)^-1 )^-1
# don't do this for anything but small n's

# ((B o Sj) o Ip)
#
#
#  B11*J11*P11 B11*J11*P12 ... B11*J11*P1p
#  B11*J11*P12 B11*J11*P22 ... B11*J11*P2p
#  ...
#  B11*J11*Pp1             ... B11*J11*Ppp
# ---------------------------
#  B11*J21*P11 B11*J21*P12 ... B11*J21*P1p
#  B11*J21*P12 B11*J21*P22 ... B11*J21*P2p
#  ...
#  B11*J21*Pp1         ... B11*J21Ppp
#
#

B.inv      <- kronecker( kronecker( solve( diag(n) - lambda0 * W),diag(J)), diag(P) )  
Sigma0.dense <- kronecker( kronecker( diag(n), SigmaJ0 ), Sigma0) 

Sigma0.spatial <- solve(B.inv %*% Sigma0.dense %*% B.inv)

# generate a deviate from the spatial distribution of Beta
#
# Beta of length P repeated J times, now repeated n times 
Beta <- rmvnorm(1, rep(Beta0,n*J), sigma=Sigma0.spatial,method="chol" ) 
Beta <- matrix(Beta,ncol=J*P,byrow=T)

### Debug ###
# take a look at the correlation of the simulated latent variable
#plot(1:64,Beta[1:64*P],type='l')
#points(1:64,Beta[1:64*P ],cex=0.5,col='red')
### End Debug ###



### loop for a fixed number of years to create the data set 
# Input:
#   Beta            covariate vector of length nJP
#   J               dim of latent variable
#   P               number or prior states 
#   n               number of observations
#   Sigma0.spatial  spatial covariance matrix 

# returned value is t, t-1, ... t-k
# where k is total years
############################ ############################ ############################

generateCropSequence <- function( Beta, J, P, n, Sigma0.spatial, SigmaE0, totalYears ) {
  Y.all <- c()

  for( simYears in 1:(totalYears - 1) ) {

    # now we have a Beta, so we can simulate x
    # we will assume that the initial conditions are independent
    #  Z0 is for the prior state
    #  Y0 is the state
    Z0 <- matrix( rnorm(J * n ),ncol=J)
    Y0 <- apply(Z0,1,function(x){ if( sum(x > 0) == 0 ){return(1)} else{return(1+which.max(x))} } )  


    #  Z1 is for the initial deviate
    Z1 <- matrix( rnorm(J * n ),ncol=J)
    Y1 <- apply(Z1,1,function(x){ if( sum(x > 0) == 0 ){return(1)} else{return(1+which.max(x))} } )  

    #  Create a design matrix according to Y0 and Beta
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
    Y10 <- (Y1*3 + Y0)-3

#print(table(Y10))
#print(table( apply( cbind(Y1,Y0), 1,paste,collapse=" <- "  ) ))

    X.design0 <- c()
    for( i in 1:length(Y10)) {
      de <- designElement( rep( Y10[i] ,J), P) 
      X.design0 <- rbind( X.design0,de  )
    }


    # at this point Z0 is distributed according to the spatial probit latent variable
    # conditioned on the prior states
    Z0 <- rowSums( X.design0 * kronecker(Beta,matrix(1,J))) + rmvnorm(n,sigma=SigmaE0) 

    Y.all <- cbind(Y.all,Y0)

    Y0 <- Y1
    Y1 <- getLatentState( Z0 ) 
  }
    
  Y.all <- cbind( Y.all, Y1) 
  colnames(Y.all) <- sprintf("%d", 1:totalYears) 

  return( Y.all )
}


######################################
#generate Crop Sequence
######################################
Y.all <- generateCropSequence( Beta, J, P, n, Sigma0.spatial, SigmaE0, 7 ) 




######################################## ######################################## ########################################
# Run MNP
######################################## ######################################## ########################################


######################################
# hyper parameters 
######################################

#Size of Transition matrix:
#  prior states
#   nonAg <- nonAg
#   nonAg <- corn
#   nonAg <- soybeans
#   corn <- nonAg
#   * corn <- corn
#   * corn <- soybeans
#   soybeans <- nonAg
#   soybeans <- corn
#   * soybeans <- soybeans
#
# stars indicate what we are interested in 
J <- 2
nPriorStates <- 2
r <- ncol(Y.all) - nPriorStates 
priorStateMatrix <- matrix( c( 1, 1, 
                               1, 2,
                               2, 1 ),ncol=2, byrow=TRUE)

P <- NROW(priorStateMatrix) + 1 

# Sigma for transition
BetaStar  <- matrix( 0 , nrow= n*J*P) 
Sigma0  <- diag(P)   # Sigma for the prior states
SigmaJ0 <- diag(J)   # Sigma for the latent variable dim in Beta prior

# hyper parameters for Sigma.inv
v <-  10 
S <- diag( J )
alpha0 <- 1



######## initial parameters ##############
Beta      <- rep(0, n*J*P)
SigInv <- diag(J)
X <- createDesignMatrix(Y.all, nPriorStates, priorStatesMatrix=priorStateMatrix, autocorrelation=T)
X.design  <- X$design 
X.map     <- X$map
Y <- X.map[,'0']
muZ       <- X.design %*% Beta
rho       <- .20 
pMinusOne <- ncol(X.design) # number of parameters - 1



##############################################################

# Imai 2005
set.seed(400)

# get initial z
z <- rnorm(r*n*J)

#### start iter here ####

for( iter in 1:10 ) {

  ############### ALPHA for Z ######################### 
    
  ## get scale parameter 
  chiScale <- alpha0^2 * sum(diag( S %*% SigInv   ))
  
  # shape, rate
  alpha.z <- chiScale / rchisq(1, df=(v*(pMinusOne)) ) 
  print(sprintf("alpha.z = %f", alpha.z)) 
  
  ############### Z ########################## 

  XBeta <- X.design %*% Beta 
  
  for( k in 1:r ) { # for unit with prior states
    for( i in 1:n ) {
      zStart <- n*(k-1) + (i*J-1)  # start of z indexes of latent vriable for obs i at time r 
      zStop  <- n*(k-1) + (i*J)    # end of z indexes of latent variable for obs i at time r
      zOut   <- n*(k-1) + i*J -2 + p # output for index i parameter p  

      zRange <- zStart:zStop

      for( p in 1:J) {
  
        cvar <- 1/SigInv[p,p]
       
        cmean <- XBeta[zRange][p] -  (  z[zRange] - XBeta[(i*J-1):(i*J)])[-p] * SigInv[p,-p] * cvar
    
        zMax = max( max( z[zRange][-p] ),0)
                        
        if( Y[n*(k-1) + i]  == p +1 ) {
          z[zOut] <- TruncNorm( zMax, cmean + 1000 * sqrt(cvar), cmean, cvar, 0 ) 
        } else {
          z[zOut] <- TruncNorm( cmean - 1000 * sqrt(cvar), zMax, cmean, cvar, 0 ) 
        }
  
      }
    }
  }
  
  zSquiggle <- z * sqrt(alpha.z) 
  

  ################# Generate Lambda ###################### 
  B.inv      <- kronecker( kronecker( solve( diag(n) - lambda0 * W),diag(J)), diag(P) )  
  Sigma0.dense <- kronecker( kronecker( diag(n), SigmaJ0 ), Sigma0) 

  Sigma0.spatial <- solve(B.inv %*% Sigma0.dense %*% B.inv)
    
  ############### ALPHA for Beta ######################### 
  # this is pretty terrible btw and needs to be reworked - jjl

  SigInv.chol <- chol(SigInv)
  X.design.trans <- t(X.design) # get things into column major form
 
  XSigInv <- c() 
  for( k in 1:r ) { # for unit with prior states
    for( i in 1:n ) {
      zStart <- n*(k-1) + (i*J-1)  # start of z indexes of latent vriable for obs i at time r 
      zStop  <- n*(k-1) + (i*J)    # end of z indexes of latent variable for obs i at time r
      zRange <- zStart:zStop

      XSigInv <- rbind( XSigInv , SigInv.chol %*% t(X.design.trans[,zRange] ))
    }
  }

  ZSigInv <-  matrix( c(chol(SigInv) %*% matrix(zSquiggle,nrow=J)), ncol=1 )

  SSX <- t( XSigInv ) %*% XSigInv
  SSY <- t( XSigInv ) %*% ZSigInv

  SSXA.inv <- solve( SSX + Sigma0.spatial )  

  BetaHat <- SSXA.inv %*% SSY
  XSigInvBetaHat <- XSigInv %*% BetaHat  


  chiScale2 <- t(ZSigInv) %*% ZSigInv + t(XSigInvBetaHat) %*% XSigInvBetaHat + t(BetaHat) %*% Sigma0.spatial %*% BetaHat + chiScale 


  #XStar = cbind( XSigInv , c(chol(SigInv) %*% matrix(zSquiggle,nrow=J)) )
  #SS = t(XStar) %*% XStar

  # sweep to get beta hat etc...
  #SS.sweep <- SWP(SWP(SS,0),1)
  #chiScale2 <- SS.sweep[J+1,J+1] + chiScale 
  
  alpha.Beta <- c(chiScale2 / rchisq(1,df=(n + v)*(pMinusOne ) ) )
  print(sprintf("alpha.Beta = %f", alpha.Beta)) 
  
  ############ Beta ######################
  
  #beta.var <- SS.sweep[1:J,1:J] * -1 * alpha.Beta
  #beta.mu  <- SS.sweep[1:J,J+1]

  beta.mu <- BetaHat
  beta.var <- alpha.Beta * SSXA.inv
  
  BetaSquiggle <- RMVN( beta.mu, beta.var ) 
  
  ############ Wishart ######################
  
  epsilon <- zSquiggle - X.design %*% BetaSquiggle
  epsilon <- matrix(epsilon,ncol=J,byrow=TRUE)
  
  R <- matrix(0,nrow=J,ncol=J)
  for(i in 1:(r*n)) R <- R +  t(epsilon[i,,drop=FALSE]) %*% epsilon[i,,drop=FALSE]
  
  
  SigInv <- RWISH( solve(R + S* alpha0^2), v + n*r)
  Sigma <- solve(SigInv)
  
  ############ Rescale ######################
  
  Beta <- BetaSquiggle/sqrt(alpha.Beta)
  
  alpha.Sigma <- mean(diag(Sigma))
  Sigma <- Sigma / alpha.Sigma
  SigInv <- SigInv * alpha.Sigma 
  z <- zSquiggle / sqrt(alpha.Sigma)

  z.save <- z
 
} 



stop('hi2u')



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






