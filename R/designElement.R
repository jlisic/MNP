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


