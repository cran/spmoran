
resf  	<- function( y, x = NULL, xgroup = NULL, nvc = FALSE, nvc_sel = TRUE,
                    nvc_num = 5, meig, method = "reml", penalty = "bic" ){

  res   <- resf_vc( y = y, xconst = x, xgroup = xgroup,
                    x = NULL, x_sel = FALSE,
                    x_nvc = FALSE, x_nvc_sel = FALSE,
                    xconst_nvc = nvc, xconst_nvc_sel = nvc_sel,nvc_num = nvc_num,
                    meig = meig, method = method, penalty = penalty, maxiter = 30 )

  b     <- res$c
  b_g   <- res$b_g

  c_vc  <- res$c_vc
  cse_vc<- res$cse_vc
  ct_vc <- res$ct_vc
  cp_vc <- res$cp_vc

  s     <- res$s
  s_c   <- res$s_c
  s_g   <- res$s_g
  e     <- res$e
  r     <- res$other$r
  sf    <- res$b_vc
  pred  <- res$pred
  resid <- res$resid
  other <- res$other
  vc    <- res$vc

  sf_alpha<-res$other$sf_alpha[1]
  x_id  <- res$other$xf_id
  par0  <- res$other$res_int$par
  nx    <- length(b[,1])
  df    <- res$other$df
  bias  <- res$other$bias


  other	  <- list( sf_alpha= sf_alpha, x_id = x_id, model = "resf", par0 = par0, nx = nx, df = df, bias=bias, res=res,
                   x = res$other$xconst, coords = meig$other$coords )

  result  <-list( b = b, b_g = b_g, c_vc=c_vc, cse_vc=cse_vc, ct_vc = ct_vc, cp_vc = cp_vc,
                  s = s, s_c=s_c, s_g = s_g, e = e, vc = vc, r = r, sf = sf, pred = pred,
                  resid = resid, other = other, call = match.call() )

  class( result ) <- "resf"
    return( result )
}

print.resf <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  if( !is.null(x$c_vc) ){
    cat("\n----Non-spatially varying coefficients (summary)----\n")
    cat("\nCoefficients:\n")
    xx2 <- data.frame(x$b$Estimate[1],x$c_vc)
    names(xx2)[1]<-"Intercept"
    print( summary( xx2 ) )
    cat("\nStatistical significance:\n")
    cp01<-apply(cbind(x$b$p_value[1],x$cp_vc),2,function(x) sum(x<0.01))
    cp05<-apply(cbind(x$b$p_value[1],x$cp_vc),2,function(x) sum(x<0.05)) - cp01
    cp10<-apply(cbind(x$b$p_value[1],x$cp_vc),2,function(x) sum(x<0.10)) - cp01 - cp05
    cp90<-length(x$cp_vc[,1]) - cp01 - cp05 - cp10
    cpv <-data.frame(rbind( cp90, cp10, cp05, cp01))
    names(cpv)[1]<-"Intercept"
    row.names(cpv) <- c("Not significant", "Significant (10% level)",
                        "Significant ( 5% level)","Significant ( 1% level)")
    print(cpv)

  } else {
    cat("\n----Coefficients------------------------------\n")
    print(x$b)
  }

  cat("\n----Variance parameter------------------------\n")
  cat("\nSpatial effects (residuals):\n")
  print(x$s)
  if( is.null(x$s_c) == FALSE ){
    cat("\nNon-spatially varying coefficients:\n")
    print(x$s_c)
  }
  if( is.null(x$s_g) == FALSE ){
    cat("\nGroup effects:\n")
    print(x$s_g)
  }
  cat("\n----Error statistics--------------------------\n")
  print(x$e)
  invisible(x)
}




