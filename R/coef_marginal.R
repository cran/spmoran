coef_marginal  <- function( mod ){

  #if(class( mod ) !="resf" ) stop("Error: The model is not an output from the resf fucntion")
  if( !inherits( mod, "resf" ) ) stop("Error: The model is not an output from the resf fucntion")

  n     <- length( mod$other$y )
  dif   <- mod$other$dif

  if( dim(mod$b)[1] > 1 ){
    nc       <- dim(mod$b)[1]
    c2       <- matrix( 0, nrow = n, ncol = nc )

    if( is.null(mod$c_vc) ){
      for(ic in 2:nc){
        c2[,ic]<- mod$b$Estimate[ic] * mod$other$dif
      }
      c2[,1]   <- NA
      c2       <- as.data.frame(c2)
      c2[,1]   <-as.logical(c2[,1])

    } else {
      c2       <- data.frame( cbind( NA, mod$c_vc * mod$other$dif ) )
    }

    names(c2)<-rownames(mod$b)

  } else {
    c2 <-NULL
  }

  result    <- list( b = c2, call = match.call() )
  class( result ) <- "coef_marginal"
  return( result )
}

print.coef_marginal <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Marginal effects from x (dy_i/dx_i) (summary)-------\n")
  bb<- x$b
  print( summary( bb ) )
  cat("\n Note: Medians are recommended summary statistics\n")
  invisible(x)
}


