# it is assumed that the matrix x has columns that are increasing with respect to year
# the result is a list of prior years from left to right t-1 t-2 ... t-k
# the prior states matrix is assumed to be in the order t-1 to t-k.
createDesignMatrix <- function(x, nPriorStates, priorStatesMatrix=NULL, deleteNoMatch=FALSE, autocorrelation=FALSE, noMatrix=FALSE) {

  n <- NROW(x) 
  m <- length(x)/n
  J <- max(x) - 1 
  r <- m - nPriorStates

  X.new <- c()
  for( i in 1:r ) {
    
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
    priorStatesMatrix <- matrix(1:(J+1),J+1)
    if( nPriorStates > 1 ) {
      for( i in 2:nPriorStates ) 
        priorStatesMatrix <- cbind( kronecker( priorStatesMatrix, matrix(1,J+1,1) ), rep(1:(J+1), (J+1)^(i-1)) )
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


  ## generate design matrix


  


  # create beta to Y mapping   

  if( autocorrelation ) {
  X.design <- 
        rep(rep(rep( (0:(n-1))*P*J, r)),each=J) + rep( (X.new$p -1)*J ,each=J) + rep(1:J,r*n)      
  } else {
         print("sorry not implemented yet") 
  }

  # create matrix
  if( !noMatrix ) {
    a <- matrix(0,n*m*J,ncol=n*J*P)
    X.design< - a[cbind( 1:length(X.design),X.design)] <- 1
  }

  return(list(map=X.new, design=X.design))
}








