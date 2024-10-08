
resf  	<- function( y, x = NULL, xgroup = NULL, weight = NULL, offset = NULL,
                    nvc = FALSE, nvc_sel = TRUE, nvc_num = 5, meig,
                    method = "reml", penalty = "bic", nongauss = NULL ){

  if( !is.null( x ) ){
    x	  <- as.matrix( x )
    if( is.numeric( x ) == FALSE ){
      mode( x ) <- "numeric"
    }
    x_id0<- apply( x, 2, sd ) != 0
    x	    <- x[ , x_id0 ]
    x     <- as.matrix(cbind( 1, x))
  } else {
    x_id0<- NULL
  }
  res   <- resf_vc( y = y, xconst = x, xgroup = xgroup, weight = weight, offset = offset,
                    x = NULL, x_sel = FALSE, x_nvc = FALSE, x_nvc_sel = FALSE,
                    xconst_nvc = nvc, xconst_nvc_sel = nvc_sel,nvc_num = nvc_num,
                    meig = meig, method = method, penalty = penalty, maxiter = 30,
                    nongauss = nongauss )

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
  r     <- as.matrix(res$other$r)
  sf    <- res$b_vc
  pred  <- res$pred
  pred_quantile<-res$pred_quantile
  resid <- res$resid
  other <- res$other
  vc    <- res$vc

  tr_par   <- res$tr_par
  tr_bpar  <-res$tr_bpar
  tr_y     <-res$tr_y
  pdf      <-res$pdf
  skew_kurt<-res$skew_kurt
  sf_alpha<-res$other$sf_alpha[1]
  x_id  <- res$other$xf_id
  par0  <- res$other$res_int$par
  nx    <- length(b[,1])
  df    <- res$other$df
  bias  <- res$other$bias

  tr_num    <- res$other$tr_num
  y_nonneg  <- res$other$y_nonneg
  y_type    <- res$other$y_type
  y_added   <- res$other$y_added
  y0        <- res$other$y0
  y         <- res$other$y
  jackup    <- res$other$jackup
  offset    <- res$other$offset
  e_NULL    <- res$other$e_NULL
  w_scale   <- res$other$w_scale
  xg_levels <- res$other$xg_levels
  B_covs    <- res$other$B_covs
  sig       <- res$other$sig
  is_weight <- res$other$is_weight
  eevSqrt   <- res$other$eevSqrt
  sig_org   <- res$other$sig_org

  b_s       <-res$other$b_s
  evSqrts   <-res$other$evSqrts
  evSqrts_n <-res$other$evSqrts_n
  evSqrts_t <-res$other$evSqrts_t
  evSqrts_tc<-res$other$evSqrts_tc
  evSqrts_t_int <-res$other$evSqrts_t_int
  evSqrts_tc_int<-res$other$evSqrts_tc_int

  omit_list <-res$other$omit_list
  idd       <-res$other$idd
  interact_sel_all<-res$other$interact_sel_all
  int_ev_sel_list <-res$other$int_ev_sel_list
  id_nsv<-res$other$id_nsv
  id_ntv<-res$other$id_ntv
  id_ntcv<-res$other$id_ntcv
  id_ntv_interact<-res$other$id_ntv_interact
  id_ntcv_interact<-res$other$id_ntcv_interact
  nev0            <-res$other$nev0

  other	  <- list( sf_alpha = sf_alpha, x_id = x_id, x_id0=x_id0, model = "resf", par0 = par0, nx = nx, df = df, bias=bias, res=res,
                   x = res$other$xconst, coords = meig$other$coords, dif=res$other$dif,method=method,
                   tr_num = tr_num, y_nonneg = y_nonneg, y_added = y_added, y_type = y_type,eevSqrt = eevSqrt,
                   xg_levels = xg_levels, is_weight = is_weight, B_covs = B_covs, sig = sig, sig_org=sig_org,
                   y0 = y0, y = y, jackup=jackup, offset=offset, e_NULL = e_NULL, w_scale = w_scale, b_s=b_s,

                   coords_z = meig$other$coords_z,
                   evSqrts=evSqrts,evSqrts_n=evSqrts_n,evSqrts_t=evSqrts_t,evSqrts_tc=evSqrts_tc,
                   evSqrts_t_int=evSqrts_t_int,evSqrts_tc_int=evSqrts_tc_int,omit_list=omit_list,idd=idd,
                   interact_sel_all=interact_sel_all,int_ev_sel_list=int_ev_sel_list,id_nsv=id_nsv,id_ntv=id_ntv,
                   id_ntcv=id_ntcv,id_ntv_interact=id_ntv_interact,id_ntcv_interact=id_ntcv_interact,
                   nev0=nev0)

  result  <-list( b = b, b_g = b_g, c_vc=c_vc, cse_vc=cse_vc, ct_vc = ct_vc, cp_vc = cp_vc,
                  s = s, s_c=s_c, s_g = s_g, e = e, vc = vc, r = r, sf = sf, pred = pred,
                  pred_quantile=pred_quantile, tr_par = tr_par, tr_bpar = tr_bpar, tr_y = tr_y,
                  resid = resid, pdf = pdf, skew_kurt = skew_kurt, other = other,
                  call = match.call() )

  class( result ) <- "resf"
    return( result )
}

print.resf <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  if( !is.null(x$c_vc) ){
    cat("\n----Non-spatially varying coefficients on x (summary) ----\n")
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
    cat("\nNon-spatial effects (coefficients on x):\n")
    print(x$s_c)
  }
  if( is.null(x$s_g) == FALSE ){
    cat("\nGroup effects:\n")
    print(x$s_g)
  }


  #if( !is.null(x$skew_kurt)|!is.null(x$tr_bpar) ){
  if( (x$skew_kurt$Estimates[1]!=0)|(x$skew_kurt$Estimates[1]!=0)|!is.null(x$tr_bpar) ){
    cat("\n----Estimated probability distribution of y--------------\n")
    if( !is.null(x$skew_kurt) ) print(x$skew_kurt)
    #if( !is.null(x$tr_bpar) ){
    #  cat( paste("(Box-Cox parameter: ", format(x$tr_bpar[1], digits=7),")\n",sep="") )
    #}
    if( !is.null(x$tr_bpar) & x$other$y_type =="continuous" ){
      cat( paste("(Box-Cox parameter: ", format(x$tr_bpar[1], digits=7),")\n",sep="") )
    } else if( x$other$y_type =="count" ){
      cat( paste("(dispersion parameter: ", format(x$other$sig_org, digits=7),")\n",sep="") )
    }
  }
  cat("\n----Error statistics--------------------------\n")
  print(x$e)

  loglik_NULL<- x$other$e_NULL[[1]][1,1]
  AIC_NULL   <- x$other$e_NULL[[1]][2,1]
  BIC_NULL   <- x$other$e_NULL[[1]][3,1]
  mod_NULL   <- x$other$e_NULL[[2]]
  if( x$other$y_type =="continuous" ){
    ml_name    <- ifelse( x$other$method=="reml", "(r)loglik: ", "loglik: " )
    cat( paste('\nNULL model: ', mod_NULL, sep="") )
    cat(paste("\n   ",ml_name,format(loglik_NULL,digits=7),sep=""))
    cat(paste(" ( AIC: ",
              format(AIC_NULL,digits=7), ",  BIC: ", format(BIC_NULL,digits=7)," )\n",sep=""))

  } else if( x$other$y_type =="count" ){
    ml_name0 <- ifelse( x$other$method=="reml", "(r)loglik", "loglik" )
    ml_name  <- paste("\n   Gaussian ",ml_name0, " approximating the model: ",sep="")
    cat( paste('\nNULL model: ', mod_NULL, sep="") )
    cat(paste(ml_name,format(loglik_NULL,digits=7), "\n",sep=""))
    cat(paste("   ( AIC: ",
              format(AIC_NULL,digits=7), ",  BIC: ", format(BIC_NULL,digits=7)," )\n",sep=""))
  }

  if( x$other$method=="reml"){
    cat('\nNote: AIC and BIC are based on the restricted/marginal likelihood.')
    cat('\n      Use method="ml" for comparison of models with different fixed effects (x)\n')
  }
  invisible(x)
}




