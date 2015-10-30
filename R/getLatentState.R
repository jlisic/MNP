# get category from a latent variable using probit rules
getLatentState <- function( z ) {
  return(
    apply(z,1,function(x){ if( sum(x > 0) == 0 ){return(1)} else{return(1+which.max(x))} } )  
  )
}


