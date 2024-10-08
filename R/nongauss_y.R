

nongauss_y<- function( y_type = "continuous", y_nonneg = FALSE, tr_num = 0){
  nongauss <-list( y_type = y_type, tr_num = tr_num, y_nonneg = y_nonneg, call = match.call() )
  class( nongauss ) <- "nongauss_y"
  return( nongauss=nongauss )
}

print.nongauss_y <- function(x, ...)
{
  if(x$y_type=="continuous"){
    if( (x$tr_num == 0) & (x$y_nonneg == FALSE) ){
      cat("y ~ N( xb,  sig )\n")
      cat(" \n")
      cat(" - N(): Normal distribution\n")
      cat(" - xb : Regression term with fixed and random coefficients in b\n")
      cat("        which is specified by resf or resf_vc function\n")
      cat(" - sig: Variance parameter\n")
    } else if( (x$tr_num == 0) & (x$y_nonneg == TRUE) ){
      cat("Box-cox transformation is applied to y to estimate\n")
      cat("y ~ P( xb, par )   (or f(y, par)~N(xb, sig) where f() denotes transformation)\n")
      cat(" \n")
      cat(" - P(): Distribution optimized through the transformation\n")
      cat(" - xb : Regression term with fixed and random coefficients in b\n")
      cat("        which is specified by resf or resf_vc function\n")
      cat(" - par: Parameter estimating data distribution\n")
    } else if( x$tr_num > 0 ){
      if( x$y_nonneg == FALSE ){
        mess0<-ifelse(x$tr_num==1,"1 SAL transformation is",
                      paste(x$tr_num," SAL transformation functions are",sep=""))
      } else {
        mess0<-paste("Box-Cox and ",x$tr_num," SAL transformations are",sep="")
      }
      mess <-paste(mess0, " applied to y to estimate",sep="")
      cat(paste(mess,"\n",sep=""))
      cat("y ~ P( xb, par )   (or f(y,par)~N(xb, sig), where f() denotes transformation)\n")
      cat(" \n")
      cat(" - P(): Distribution optimized through the transformations\n")
      cat(" - xb : Regression term with fixed and random coefficients in b\n")
      cat("        which is specified by resf or resf_vc function\n")
      cat(" - par: Parameters estimating data distribution\n")
    }

  } else if( x$y_type == "count" ){
    if( x$tr_num == 0 ){
      cat("Overdispersed Poisson distribution (Log-Gaussian approximation)\n")
      cat("y ~ oPois( mu, sig ), mu = exp( xb )\n")
      cat(" \n")
      cat(" - oPois(): Overdispersed Poisson distribution\n")
      cat(" - xb     : Regression term with fixed and random coefficients in b\n")
      cat("            which is specified by resf or resf_vc function\n")
      cat(" - sig    : Dispersion parameter (overdispersion if sig > 1)\n")

    } else {
      mess0<-paste("Distribution that transformes overdispersed Poisson distribution to adapt to y using ",x$tr_num," SAL trans.",sep="")
      cat(paste(mess0,"\n",sep=""))
      cat(" y ~ P( mu, par ), mu = exp( xb )\n")
      cat(" \n")
      cat(" - P(): Distribution for count data optimized through the transformations\n")
      cat(" - xb : Regression term with fixed and random coefficients in b\n")
      cat("        which is specified by resf or resf_vc function\n")
      cat(" - par: Parameters estimating data distribution\n")
    }
  }
  invisible(x)
}



