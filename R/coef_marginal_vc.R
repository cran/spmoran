coef_marginal_vc  <- function( mod ){

  if(class( mod ) !="resf_vc" ) stop("Error: The model is not an output from the resf_vc fucntion")

  #if( class( mod )== "resf_vc" ){
    n     <- length( mod$other$y )
    dif   <- mod$other$dif

    b_vc2 <- mod$b_vc * mod$other$dif
    b_vc2[,1]<-NA
    if( !is.null(mod$B_vc_s) ){
      B_vc_s2 <- mod$B_vc_s[[1]] * mod$other$dif
      B_vc_s2[,1]<-NA
    } else {
      B_vc_s2 <-NULL
    }
    if( !is.null(mod$B_vc_n) ){
      B_vc_n2 <- mod$B_vc_n[[1]] * mod$other$dif
      B_vc_n2[,1]<-NA
    } else {
      B_vc_n2 <- NULL
    }

    if( !is.null(mod$c_vc) ){
      c2       <- mod$c_vc * mod$other$dif
      c2       <- as.data.frame(c2)
      names(c2)<-rownames(mod$c)

    } else if( !is.null(mod$c) ){
      nc       <- dim(mod$c)[1]
      c2       <- matrix( 0, nrow = n, ncol = nc )
      for(ic in 1:nc){
        c2[,ic]<- mod$c$Estimate[ic] * mod$other$dif
      }
      c2       <- as.data.frame(c2)
      names(c2)<-rownames(mod$c)

    } else {
      c2 <-NULL
    }

    other     <- mod$other$s_c
    result    <- list( b_vc = b_vc2, B_vc_s = B_vc_s2, B_vc_n = B_vc_n2, c = c2, other = other,call = match.call() )
    class( result ) <- "coef_marginal_vc"
  #}

  return( result )
}


print.coef_marginal_vc <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)

  cat("\n----Marginal effects from x (dy_i/dx_i) (summary)----\n")
  print( summary( x$b_vc ) )

  if( is.null(x$other) & !is.null(x$c) ){
    cat("\n----Marginal effects from xconst (dy_i/dx_i)(summary)----\n")
    print( summary( x$c ) )

  } else if( !is.null(x$c_vc) ){
    cat("\n----Marginal effects from xconst (dy_i/dx_i)(summary)----\n")
    print( summary( x$c_vc ) )
  }
  cat("\n Note: Medians are recommended summary statistics\n")
  invisible(x)
}
