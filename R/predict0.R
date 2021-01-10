predict0   <- function( mod, meig0, x0 = NULL, xgroup0 = NULL ){

  {

    af <-function(par,y) par[1]+par[2]*y
    sc <-function(par,y) (y-par[1])/par[2]
    sa <-function(par,y) sinh(par[1]*asinh(y)-par[2])
    bc <-function(par,y, jackup){

      yy     <- y + ifelse(jackup, abs( par[2] ), 0 )
      if(par[1]==0){
        res <- log(yy)
      } else {
        res <- ( yy^par[1] - 1 ) / par[1]
      }
    }
    sal<-function(par,y, noconst=FALSE){
      if(noconst==FALSE){
        par[3]+ par[4]*sinh(par[1]*asinh(y)-par[2])
      } else {
        sinh(par[1]*asinh(y)-par[2])
      }
    }

    d_af<-function(par,y) par[2]
    d_sc <-function(par,y) 1/par[2]
    d_sa<-function(par,y) par[1]*cosh(par[1]*asinh(y)-par[2])/sqrt(1+y^2)
    d_bc<-function(par,y,jackup) (y + ifelse(jackup, abs( par[2] ), 0 ) )^(par[1]-1)#abs(y)^(par-1)
    d_sal<-function(par,y,noconst=FALSE){
      if(noconst==FALSE){
        d_af(par=par[3:4],sa(par=par[1:2],y))*d_sa(par=par[1:2],y)
      } else {
        d_af(par=c(0, 1),sa(par=par[1:2],y))*d_sa(par=par[1:2],y)
      }
    }

    i_af<-function(par,y) (y-par[1])/par[2]
    i_sc <-function(par,y) par[1]+par[2]*y
    i_sa<-function(par,y) sinh(1/par[1]*(asinh(y)+par[2]))
    i_bc<-function(par,y,jackup){
      if(par[1]==0){
        exp(y) - ifelse(jackup, abs( par[2] ), 0 )
      } else {
        (par[1]*y + 1)^(1/par[1]) - ifelse(jackup, abs( par[2] ), 0 )
      }
    }

    i_sal<-function(par,y,noconst=FALSE){
      if(noconst==FALSE){
        i_sa(par=par[1:2],i_af(par=par[3:4],y))
      } else {
        i_sa(par=par[1:2],i_af(par=c(0, 1),y))
      }
    }

    d_sal_k   <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL,y_ms=NULL,z_ms,jackup){
      d       <-1
      if( !is.null( bc_par[1] ) ){
        d     <- d_bc(bc_par,y,jackup=jackup)*d
        y     <-   bc(bc_par,y,jackup=jackup)
      }

      if(k>0){
        d       <- d_sc(par=y_ms,y)*d
        y       <-   sc(par=y_ms,y)

        for(kk in 1:k){
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          d     <-d_sal(par[[kk]],y,noconst=nc_dum)*d
          y     <-  sal(par[[kk]],y,noconst=nc_dum)
        }
      }
      d       <-d_sc(par=z_ms,y)*d
      return(d)
    }
    i_sal_k   <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL,y_ms=NULL,z_ms, jackup){

      y       <-i_sc(par=z_ms,y)
      if(k>0){
        for(kk in k:1) {
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          y <- i_sal(par[[kk]],y,noconst=nc_dum)
        }
      }

      if( !is.null( y_ms ) ){
        y     <-i_sc(par=y_ms,y)
      }

      if( !is.null(bc_par[1]) ){
        y <- i_bc(par=bc_par, y, jackup=jackup)
      }

      return(y)
    }
    sal_k     <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL, jackup){

      if( !is.null(bc_par[1]) ){
        y   <- bc(par=bc_par,y,jackup=jackup)
      }

      if( k > 0){
        y_ms<- c( mean(y), sd(y) )
        y   <- sc(par=y_ms,y)

        for(kk in 1:k) {
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          y <- sal(par[[kk]],y,noconst=nc_dum)
        }
      } else {
        y_ms<- NULL
      }

      z_ms  <- c( mean(y), sd(y) )
      y     <- sc(par=z_ms,y)

      return(list(y=y,y_ms=y_ms, z_ms=z_ms))
    }

    ######## Negative log-likelihood (SAL distribution)
    NLL_sal<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){
      n   <-length(y)
      par2<-list(NULL)
      for(kk in 1:k){
        if(noconst_last & (kk==k)){
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+2)]

        } else {
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par[(4*(kk-1)+4)]<- abs( par[(4*(kk-1)+4)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+4)]
        }
      }

      if(tr_nonneg){
        np_b  <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp<- sum( log( d_abs ) )
      nll <- ee/2 - comp
    }

    NLL_sal2<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){
      n   <-length(y)
      par2<-list(NULL)
      for(kk in 1:k){
        if(noconst_last & (kk==k)){
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+2)]
        } else {
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par[(4*(kk-1)+4)]<- abs( par[(4*(kk-1)+4)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+4)]
        }
      }

      if(tr_nonneg){
        np_b    <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp<- sum( log( d_abs ) )
      loglik<- ee/2 - comp
      return(list(z=z, b=b, loglik=loglik,comp=comp,y_ms=y_ms,z_ms=z_ms))
    }

    ######## Negative log-likelihood (BC distribution)
    NLL_bc <-function(par,y,M,Minv,m0,jackup){
      par[2]<-abs(par[2])
      z0  <- bc(par=par,y=y,jackup=jackup)
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]#z0#
      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,y ) )#abs( d_bc( par=par,y ) )#
      comp<- sum( log( d_abs ) )
      nll <- ee/2 - comp
    }
    NLL_bc2<-function(par,y, M, Minv,m0,jackup){
      par[2]<-abs(par[2])
      z0  <- bc(par=par,y=y,jackup=jackup)
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]#z0#
      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,y ) )#abs( d_bc( par=par,y ) )#
      comp<- sum( log( d_abs ) )
      loglik<- ee/2 - comp
      return(list(z=z, b=b, loglik=loglik,comp=comp,z_ms=z_ms))#
    }

    lik_resf_vc_tr <- function( par, k, tr_nonneg=FALSE,noconst_last=TRUE, X,
                                evSqrt, y0, n, nx, emet,M0,Minv,term1,null_dum3=NULL,jackup ){
      tr_par       <- par
      if(k > 0){
        tr_par2<-list(NULL)
        for(kk in 1:k){
          if(noconst_last & (kk==k)){
            tr_par2[[kk]]        <- tr_par[(4*(kk-1)+1):(4*(kk-1)+2)]
            tr_par2[[kk]][1]     <- abs(tr_par2[[kk]][1])
          } else {
            tr_par2[[kk]]        <- tr_par[(4*(kk-1)+1):(4*(kk-1)+4)]
            tr_par2[[kk]][c(1,4)]<- abs(tr_par2[[kk]][c(1,4)])
          }
        }

        if(tr_nonneg){
          np_b  <-length(tr_par)
          bc_par<-tr_par[(np_b-1):np_b]
          #if(bc_par[1] < -0.5) bc_par[1] <- -0.5###################### added 2020/12/14
          bc_par[2]<-abs(bc_par[2])
        } else {
          bc_par<-NULL
        }
        z0      <- sal_k(par=tr_par2,y=y0,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
        z       <- z0$y
        z_ms    <- z0$z_ms
        y_ms    <- z0$y_ms
        d_abs   <- abs( d_sal_k( par=tr_par2 ,y=y0, k=k,noconst_last=noconst_last,
                                 bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup ) )
        comp    <- (-2)*sum( log( d_abs ))

        if(max(abs(z_ms)) > 10^10){######### deleted
          comp  <- 10^100
        }

      } else if(tr_nonneg){
        z0      <- bc(par=tr_par,y=y0,jackup=jackup)
        z_ms    <- c(mean(z0),sd(z0))
        z       <- (z0-z_ms[1])/z_ms[2]
        d_abs   <- abs( d_bc(par=tr_par,y0,jackup=jackup)*d_sc(par=z_ms,y0) )
        comp    <- (-2)*sum( log( d_abs ) )
      } else {
        z       <- y0
        comp    <- 0
      }

      m	        <- crossprod( X[,null_dum3], z )
      zz        <- sum( z^2 )

      ########################## reml
      m[-(1:nx)]		<- m[ -( 1:nx ) ] * evSqrt[null_dum3[-( 1:nx )]]
      b		<- Minv %*% m
      sse		<- zz - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
      dd		<- abs(sse) + sum( b[ -( 1:nx ) ] ^ 2 )
      if( emet == "reml" ){
        term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
      } else if( emet == "ml" ){
        term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
      }
      loglik		<- term1 + term2 + comp

      if(is.nan(loglik)) loglik <-10^100
      if(is.na(loglik))  loglik <-10^100

      return( loglik[ 1 ] )
    }

    lm_cw<-function(y, M, Minv, m0, k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){

      if(k != 0){
        par00<-rep(c(1,0,0,1),k)
        #lower<-rep(c(1e-10,-3,-Inf,1e-10),k)###### check!!!
        #upper<-rep(c(Inf, 3, Inf, Inf),k)

        if(noconst_last){
          par00<-par00[1:(4*k-2)]
          #lower<-lower[1:(4*k-2)]
          #upper<-upper[1:(4*k-2)]
        }

        if(tr_nonneg){
          par00<-c(par00,1, 0.001)
          #lower<-c(lower, -1, 0)
          #upper<-c(upper,  5, 10)
        }
        res    <-optim(par=par00,fn=NLL_sal, y=y,M=M,Minv=Minv,m0=m0,k=k,
                       noconst_last=noconst_last,tr_nonneg=tr_nonneg,jackup=jackup)
        #method="L-BFGS-B",lower=lower,upper=upper)
        est0   <-res$par
        est    <-list(NULL)
        for(kk in 1:k){
          if(noconst_last & (kk==k)){
            est0[(4*(kk-1)+1)]<- abs( est0[(4*(kk-1)+1)] )
            est[[kk]]    <- est0[(4*(kk-1)+1):(4*(kk-1)+2)]
          } else {
            est0[(4*(kk-1)+1)]<- abs( est0[(4*(kk-1)+1)] )
            est0[(4*(kk-1)+4)]<- abs( est0[(4*(kk-1)+4)] )
            est[[kk]]    <- est0[(4*(kk-1)+1):(4*(kk-1)+4)]
          }
        }

        if(tr_nonneg){
          np_be      <- length(est0)
          bc_par     <- est0[(np_be-1):np_be]
          if(jackup) bc_par[2]  <- abs(bc_par[2])
          #bc_par     <- est0[length(est0)]
          est[[kk+1]]<- bc_par
        } else {
          bc_par<-NULL
        }
        res2  <-NLL_sal2(par=unlist(est), y=y, M=M, Minv=Minv, m0=m0,k=k,
                         noconst_last=noconst_last, tr_nonneg=tr_nonneg,jackup=jackup)
        b     <-res2$b
        z     <-res2$z
        loglik<-res2$loglik
        comp  <-res2$comp
        y_ms  <-res2$y_ms
        z_ms  <-res2$z_ms
        #if(tr_nonneg){
        #  est <- est[1:(length(est)-2)]
        #}

      } else if(tr_nonneg){
        #res   <-optimize(f=NLL_bc, interval=c(-5,5),y=y,M=M,Minv=Minv,m0=m0 )#c(-5,5)
        #bc_par<-res$minimum
        #est   <-NULL

        res   <-optim(par=c(1,0.001),fn=NLL_bc,y=y,M=M,Minv=Minv,m0=m0,jackup=jackup )# interval=c(-5,5)
        bc_par<-res$par
        bc_par[2]<-abs(bc_par[2])
        est   <-NULL

        res2  <-NLL_bc2(par=bc_par, y=y, M=M, Minv=Minv, m0=m0,jackup=jackup )
        b     <-res2$b
        z     <-res2$z
        loglik<-res2$loglik
        comp  <-res2$comp
        y_ms  <-NULL
        z_ms  <-res2$z_ms# this part is different from Kaihatsuchushi

      } else {
        est    <- NULL
        res2   <- NULL
        b      <- NULL
        z      <- y
        loglik <- NULL
        comp   <- 0
        bc_par <- NULL
        y_ms  <-NULL
        z_ms  <-NULL
      }

      return(list(vpar=est,bc_par=bc_par,b=b,z=z,loglik=loglik,comp=comp,
                  y_ms=y_ms,z_ms=z_ms))
    }
  }

  n<-length(mod$other$coords[,1])
  if( mod$other$model == "esf" ){

    if( is.null( dim( mod$r ) ) ){
      sf_pred	<- meig0$sf[ mod$other$sf_id] %*% c( mod$r[,1] )
    } else {
      sf_pred	<- meig0$sf[ ,mod$other$sf_id ] %*% c( mod$r[,1] )
    }

    x0        <- as.matrix(x0)
    b_g0      <- NULL
    if( is.null( mod$other$x_id ) ){
      xb_pred	<- c( as.matrix( mod$b[ 1 ] ) )
      pred	  <- xb_pred + sf_pred
    } else {
      if( is.null( x0 ) ){
        message( " Note: Only spatial component (sf) is interpolated because x0 is missing")
        pred	<- sf_pred
      } else {
        if( is.numeric( x0 ) == FALSE ){
          mode( x0 ) <- "numeric"
        }

        if( length( c(sf_pred)) ==1 ){
          xb_pred <- c(1, x0[ ,mod$other$x_id ]) %*% mod$b[, 1 ]
        } else {
          xb_pred	<- as.matrix( cbind( 1, x0[ ,mod$other$x_id ] ) ) %*% mod$b[, 1 ]
        }
        pred	<- xb_pred + sf_pred
        pred	<- as.data.frame( cbind( pred, xb_pred, sf_pred ) )
        names( pred )<- c( "pred", "xb", "sf_residual" )
      }
    }
    c_vc  <- cse_vc <- ct_vc <- cp_vc <- NULL


  } else {

    if( is.null( dim( mod$r ) ) ){
      sf_pred	<- meig0$sf %*% c( mod$r )
    } else {
      sf_pred	<- meig0$sf %*% c( mod$r[, 1 ] )
    }
    n0        <- length( c(sf_pred) )

   if( !is.null( x0 )){
    if(dim( as.matrix( x0 ) )[2] != length( mod$other$res$other$xf_id )){
      stop("x and x0 must have the same number of columns")
    }

    X1 <- as.matrix( mod$other$res$other$xconst )[,mod$other$res$other$xf_id]
    X0 <- as.matrix( as.matrix( x0 )[,mod$other$res$other$xf_id] )
   } else {
     X0<-NULL
   }

   XX_0	  <- list( NULL )
   XX	    <- NULL
   nvc <- mod$other$res$other$nvc_xconst
   if( ( is.logical( nvc[ 1 ] ) ==FALSE )&( !is.null( x0 ) ) ){

     X1_nvc<- as.matrix( X1 )[ , nvc ]
     if( n0 == 1 ){
      X0_nvc<- X0[ nvc ]
    } else {
      X0_nvc<- as.matrix( X0[ , nvc ] )
    }

    xxfname	  <- names( as.data.frame( X1 ) )[ nvc ]
    nnxf      <- length( xxfname )
    X1        <- as.matrix( X1 )
    X1_nvc    <- as.matrix( X1_nvc )
    np_xx     <-apply( X1_nvc,2,function( x ) length( unique( x )))
    np_xx     <-ifelse( np_xx < mod$other$res$other$nvc_num/0.7, round( np_xx * 0.7 ) ,mod$other$res$other$nvc_num )

    np_xx_max <-round( n/nnxf ) - 2
    np_xx[ np_xx > np_xx_max ] <-np_xx_max
    np_xx[ np_xx < 2 ] <- 2

    #XB_c    <-NULL
    B_c     <-list(NULL)
    for( ii in 1:dim( X1_nvc )[ 2 ] ){
      if( np_xx[ ii ] <= 2 ){
        B_c[[ii+1]]     <- 0
        np_xx[[ii]]<- 0
      } else {

        test<-TRUE
        iiii<-0
        while(test){
          kkk      <- np_xx[ ii ]-iiii
          knots    <-seq(min( X1_nvc[ ,ii ] ),max( X1_nvc[ ,ii ] ),len=kkk+2)[2:(kkk+1)]
          testt<- try(XX1_00<- ns( X1_nvc[ , ii], knots = knots ), silent=TRUE)
          test <- class(testt)[1] == "try-error"
          iiii <- iiii+1
        }

        XX1_0 <- cbind( X1[,ii]    , XX1_00)

        if( n0 == 1 ){
          XX0_0 <- c( X0_nvc[ii], predict( XX1_00, newx= X0_nvc[ii]) )
        } else {
          XX0_0 <- cbind( X0_nvc[,ii], predict( XX1_00, newx= X0_nvc[,ii]) )
        }

        if( !is.na( mod$other$res$other$sel_basis_c[[ ii ]][ 1 ] ) ){
          XX1_0<- XX1_0[,mod$other$res$other$sel_basis_c[[ii]]]
          if( n0 == 1){
            XX0_0<- XX0_0[mod$other$res$other$sel_basis_c[[ii]]]
          } else {
            XX0_0<- XX0_0[,mod$other$res$other$sel_basis_c[[ii]]]
          }
        }

        if( n0 == 1){
          for( j in 1:length(XX0_0)){
            XX0_0[j]<- ( XX0_0[j] - mean(XX1_0[,j]))/sd(XX1_0[,j])
          }
        } else {
          for( j in 1:dim(XX0_0)[2]){
            XX0_0[,j]<- ( XX0_0[,j] - mean(XX1_0[,j]))/sd(XX1_0[,j])
          }
        }

        B_c[[ii]]    <- XX0_0
        #XB_c         <- cbind( XB_c , X0const_nvc[, ii ] * B_c[[ii]] )
      }
    }
  }

  if( is.null( mod$other$res$c_vc )|( is.null( x0 ) ) ){
    c_vc  <- cse_vc <- ct_vc <- cp_vc <- NULL
  } else {
    xc_vc	<- 0
    c_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    cse_vc<- matrix(0, nrow = n0, ncol = nnxf )
    ct_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    cp_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    for( i in 1:nnxf ){
      evSqrts_c<- mod$other$res$other$evSqrts_c[[ i ]]
      if(length( evSqrts_c ) == 1) evSqrts_c <- NULL

      if( length( mod$other$res$other$evSqrts_c[[ i ]] ) <= 1 ){#!= ne
        c_vc[ , i ]	  <- mod$other$res$other$b_c[[ i ]][ 1 ]
        cse_vc[ , i ]	<- sqrt( mod$other$res$other$b_covs_c[[ i ]] )
        ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
        cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$res$other$df )
      } else {
        if( n0 == 1 ){
          c_vc[ , i ]	<- mod$other$res$other$b_c[[ i ]][1] + c(B_c[[ i ]]) %*% c( mod$other$res$other$b_c[[ i ]][ -1 ])
          sf2		<- t( B_c[[ i ]] * mod$other$res$other$evSqrts_c[[ i ]] )
        } else {
          c_vc[ , i ]	<- mod$other$res$other$b_c[[ i ]][1] + B_c[[ i ]] %*% c( mod$other$res$other$b_c[[ i ]][ -1 ] )
          sf2		<- t( t( B_c[[ i ]] ) * mod$other$res$other$evSqrts_c[[ i ]] )
        }

        if( n0 == 1){
          x_sf		      <-t(as.matrix( c( 1, c( sf2 ) ) ))
          cse_vc[ , i ]	<- sqrt( x_sf %*% mod$other$res$other$b_covs_c[[ i ]] %*% t( x_sf ) )
        } else {
          x_sf		      <- as.matrix( cbind( 1, sf2 ) )
          cse_vc[ , i ]	<- sqrt( colSums( t( x_sf ) * ( mod$other$res$other$b_covs_c[[ i ]] %*% t( x_sf ) ) ) )
        }
        ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
        cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$res$other$df )
      }
    }
  }

  if( is.null( mod$b_g )==FALSE ){
      if( is.null(xgroup0) ){
        message( " Note: Group effects are ignored because xgroup0 is missing")
        b_g0     <- NULL
      } else {
        ng       <- length(mod$b_g)
        xgroup0  <- data.frame( xgroup0 )
        b_g0     <- xgroup0
        for(ggid in 1:ng){
          b_g_sub        <- as.data.frame( mod$b_g[[ggid]] )
          xg00_0	       <- factor( xgroup0[ , ggid ] )
          xg00_nam0      <- names( xgroup0 )[ ggid ]
          xg00_nam       <-ifelse( xg00_nam0 == "xgroup0", "xgroup", xg00_nam0 )
          xgroup0_nam    <- paste( xg00_nam, "_", xg00_0, sep = "" )
          xgroup0_dat    <- data.frame(id=1:length(xgroup0_nam),xg_nam=xgroup0_nam)
          xgroup0_dat2   <- merge(xgroup0_dat, b_g_sub, by.x = "xg_nam", by.y = "row.names", all.x = TRUE )
          b_g0[,ggid]    <- xgroup0_dat2[order(xgroup0_dat2$id), "Estimate" ]
        }
        xg_names         <- names( xgroup0 )
        if( length( xg_names ) == 1 ){
          xg_names       <- ifelse( xg_names == "xgroup0", "xgroup", paste("xgroup_",xg_names,sep="" ) )
        } else {
          xg_names       <- paste( "xgroup_", xg_names, sep = "" )
        }
        names(b_g0)      <- xg_names

      }
  } else {
      b_g0 <- NULL
  }


    if( is.null( mod$other$x_id ) ){
    	xb_pred	<- c( mod$b$Estimate[1] )#changed
    	pred	  <- xb_pred + sf_pred
    } else {
     if( is.null( X0 ) ){
    		message( " Note: Trend term (xb) is ignored because x0 is missing")
    		pred	<- data.frame(pred=NA, sf_residual=sf_pred)
   	 } else {
  		 if( is.numeric( X0 ) == FALSE ){
    			mode( X0 ) <- "numeric"
    	 }

   	   if( is.null( c_vc )){
   	     if( n0 ==1 ){
   	       xb_pred<- c( 1, X0[ mod$other$x_id ]) %*% mod$b[, 1 ]
   	     } else {
   	       xb_pred<- as.matrix( cbind( 1, X0[ ,mod$other$x_id ] ) ) %*% mod$b[, 1 ]
   	     }

   	   } else {
   	     if( n0 == 1 ){
   	       xb_pred<- sum( X0[ mod$other$x_id ] * c_vc ) + mod$b$Estimate[1]
   	     } else {
   	       xb_pred<- rowSums( X0[ ,mod$other$x_id ] * c_vc ) + mod$b$Estimate[1]
   	     }
   	   }

   	   pred	  <- xb_pred + sf_pred
   	   pred	  <- data.frame( pred = pred, xb = xb_pred, sf_residual = sf_pred )
    	}
    }

    if( is.null( b_g0 ) == FALSE ){
      pred      <- data.frame( pred, group = b_g0 )
      na_b_g0   <- is.na( b_g0 )
      if( sum( na_b_g0 ) > 0 ){
        b_g0[ na_b_g0 ] <- 0
        message( " Note: NAs are given to the groups that are not in xgroup")
      }
      pred[ ,1 ]<- pred[ ,1 ] + rowSums( b_g0 )
    }

   y0         <- mod$other$y
   tr_num     <- mod$other$tr_num
   tr_nonneg  <- mod$other$tr_nonneg
   tr_par     <- mod$tr_par
   tr_bpar    <- mod$tr_bpar$Estimate
   y_added    <- mod$other$y_added
   jackup     <- mod$other$jackup

   noconst_last<-TRUE
   if( tr_num > 0 ){######## transfer this part to prediction functions
     z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,jackup=jackup)
     z_ms     <- z0$z_ms
     y_ms     <- z0$y_ms
     pred2    <- i_sal_k(par=tr_par,y=pred[,1],k=tr_num,noconst_last=noconst_last,
                         bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup) - y_added
     if( tr_nonneg ) pred2[ pred2 < 0 ] <- 0
     pred	    <- data.frame( "pred" = pred2, pred )
     names(pred)[2]<-"pred_trans"

   } else if( tr_nonneg ){
     z0      <- bc(par=tr_bpar,y=y0,jackup=jackup)
     z_ms    <- c(mean(z0),sd(z0))
     y_ms     <- NULL
     pred0b   <- pred[,1]*z_ms[2] + z_ms[1]
     pred2     <- i_bc(par=tr_bpar,y=pred0b,jackup=jackup) - y_added
     pred2[pred2 < 0 ] <- 0
     pred	    <- data.frame( "pred" = pred2, pred )
     names(pred)[2]<-"pred_trans"

   }

  }
  return( list( pred = pred, c_vc = c_vc, cse_vc =cse_vc, ct_vc = ct_vc, cp_vc = cp_vc ) )
}

