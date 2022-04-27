predict0   <- function( mod, meig0, x0 = NULL, xgroup0 = NULL, offset0 = NULL,
                        weight0 = NULL, compute_quantile = FALSE ){

  #if( (class( mod ) !="resf")&(class( mod ) !="esf") ){
  if( (inherits(mod, "resf")==FALSE)&(inherits( mod, "esf")==FALSE) ){
      stop("Error: Input model must be an output from resf or esf function")
  }

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
    NLL_sal<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,y_nonneg=FALSE,jackup){
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

      if(y_nonneg){
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

    NLL_sal2<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,y_nonneg=FALSE,jackup){
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

      if(y_nonneg){
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

  }

  n        <-length(mod$other$coords[,1])
  nx       <-mod$other$nx
  eevSqrt  <-mod$other$eevSqrt
  #null_dum3<-mod$other$null_dum3
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
            test <- inherits( testt, "try-error" )#class(testt)[1] == "try-error"
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
      pred	  <- data.frame( pred = pred, xb = xb_pred, sf_residual = sf_pred )
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
      pred      <- data.frame( pred, b_g0 )
      na_b_g0   <- is.na( b_g0 )
      if( sum( na_b_g0 ) > 0 ){
        b_g0[ na_b_g0 ] <- 0
        message( " Note: b_g = 0 is assumed for samples who does not belong to any groups in xgroup")
      }
      pred[ ,1 ]<- pred[ ,1 ] + rowSums( b_g0 )
    }

    y0         <- mod$other$y
    tr_num     <- mod$other$tr_num
    y_nonneg   <- mod$other$y_nonneg
    y_type     <- mod$other$y_type
    tr_par     <- mod$tr_par
    tr_bpar    <- mod$tr_bpar$Estimate
    y_added    <- mod$other$y_added
    jackup     <- mod$other$jackup
    pred0      <- pred[,1]

    pred_name  <-names(pred)
    noconst_last<-TRUE
    if( tr_num > 0 ){######## transfer this part to prediction functions
      z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,jackup=jackup)
      z_ms     <- z0$z_ms
      y_ms     <- z0$y_ms
      pred2    <- i_sal_k(par=tr_par,y=pred0,k=tr_num,noconst_last=noconst_last,
                          bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup) - y_added
      if( y_nonneg ) pred2[ pred2 < 0 ] <- 0
      if( y_type=="count" ){
        pred2 <- exp( pred2 )
        if( !is.null( mod$other$offset ) ){
          if( is.null( offset0 ) ) stop( "offset0 is missing" )
          pred2 <- pred2 * offset0
        }
      }
      pred	    <- data.frame( pred= pred2, pred_transG=pred0, pred[,-1] )
      names(pred)[-c(1:2)]<-pred_name[ -1 ]

    } else if( y_nonneg ){
      z0        <- bc(par=tr_bpar,y=y0,jackup=jackup)
      pred2     <- i_bc(par=tr_bpar,y=pred0,jackup=jackup) - y_added
      pred2[ is.nan( pred2 ) ] <- 0
      pred2[pred2 < 0 ]        <- 0
      pred	    <- data.frame( pred = pred2, pred_transG=pred0, pred[,-1] )
      names(pred)[-c(1:2)]<-pred_name[ -1 ]

    } else if( y_type=="count" ){
      pred2     <- exp( pred0 )
      if( !is.null( mod$other$offset ) ){
        if( is.null( offset0 ) ) stop( "offset0 is missing" )
        pred2  <- pred2 * offset0
      }
      pred	    <- data.frame( pred = pred2, pred_transG=pred0, pred[,-1] )
      names(pred)[-c(1:2)] <- pred_name[ -1 ]
    }

    res      <- pred
    pq_dat   <- NULL
    if( compute_quantile ==TRUE ){
      if( mod$other$y_type == "count" ){
        message("Note: 'compute_quantile' is currently not supported for count data")
      } else {
        if( mod$other$is_weight ==TRUE ){
          if( is.null( weight0 ) ) stop( "Specify weight0 to compute quantile" )
        } else {
          weight0     <- NULL
        }

        if( !is.null( mod$other$x_id )&is.null( X0 ) ){
          stop( "x0 is required to compute quantile" )
        }

        B_covs<-mod$other$B_covs
        sig   <-mod$other$sig
        XX	   <- as.matrix( cbind( 1, X0, meig0$sf ) )

        #########Group
        if( !is.null( mod$b_g ) ){
          for( ggid in 1:ng ){
            skip_id     <-which(is.na(mod$b_g[[ggid]][,2]))
            xg_levels   <- mod$other$xg_levels[[ggid]][ -skip_id]
            Xg          <- matrix( 0, nrow =n0, ncol=length( xg_levels ) )
            for( ggid2 in 1:length( xg_levels ) ) Xg[ ,ggid2 ][ xgroup0[, ggid ] == xg_levels[ ggid2 ] ]<-1
            XX  <-cbind(XX, Xg)
          }
        }

        #########NVC
        if( !is.null( mod$other$res$other$evSqrts_c[[1]] ) ){
          for( i in 1:nnxf ){
            #if(i ==1){
            # evSqrts_n<- NULL
            #} else {
            evSqrts_n<- mod$other$res$other$evSqrts_c[[ i ]]
            if( length( evSqrts_n ) == 1 ) evSqrts_n <- NULL
            #}

            if( !is.null( evSqrts_n ) ){
              XX<- cbind( XX, X0[ , i ] * B_c[[ i ]] )
            }
          }
        }

        if( is.null( weight0 ) ){
          weight0  <- 1
        } else {
          weight0  <- weight0*mod$other$w_scale
        }

        X3      <- XX
        X3[,-( 1:nx )]<- t(t(XX[,-( 1:nx )])* eevSqrt[eevSqrt > 0])
        pred0_se<- sqrt( colSums( t( sqrt(weight0)*X3 ) * ( B_covs %*% t( sqrt(weight0)*X3 ) ) ) + sig )
        pred0_se<- pred0_se/sqrt( weight0 )
        if( sum(names( res ) %in% "pred_transG") == 0 ){
          res_name<-names( res )
          res     <- data.frame( res[,1], pred_se = pred0_se, res[,-1] )
          names( res ) <- c( res_name[1], "pred_se", res_name[ -1 ] )
        } else {
          res_name<-names( res )
          res     <- data.frame( res[,1:2], pred_se = pred0_se, res[,-(1:2)])
          names( res ) <- c( res_name[1:2], "pred_transG_se", res_name[-(1:2)] )
        }

        pquant    <- c(0.01, 0.025, 0.05, seq(0.1,0.9,0.1), 0.95, 0.975, 0.99)
        pq_dat0   <- NULL
        for(pq in pquant){
          pq_dat0<-cbind(pq_dat0,qnorm(pq,pred0,pred0_se))
        }
        pq_dat0       <- as.data.frame(pq_dat0)
        names(pq_dat0)<- paste("q",pquant,sep="")

        if( tr_num > 0 ){######## transfer this part to prediction functions
          tr_bpar0 <-tr_bpar
          z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
                            bc_par=tr_bpar0,jackup=jackup)
          z_ms     <- z0$z_ms
          y_ms     <- z0$y_ms
          pq_dat   <- pq_dat0
          for(pq in 1:ncol( pq_dat0 ) ){
            ptest<-try(pq_pred<- i_sal_k( par=tr_par,y=pq_dat0[,pq],k=tr_num,noconst_last=noconst_last,
                                          bc_par=tr_bpar0,y_ms=y_ms,z_ms=z_ms,jackup=jackup ) - y_added )
            if( !inherits(ptest, "try-error") ){#class(ptest)!="try-error"
              pq_dat[,pq]       <-pq_pred
            } else {
              pq_dat[,pq]       <-NA
            }
          }

        } else if( y_nonneg ==TRUE ){
          y        <- bc(par=tr_bpar,y=y0,jackup=jackup)
          pred     <- i_bc(par=tr_bpar,y=pred0,jackup=jackup) - y_added
          pred[ is.nan( pred ) ] <- 0
          pred[pred < 0 ]        <- 0

          pq_dat   <- pq_dat0
          for(pq in 1:ncol(pq_dat0)){
            ptest<-try(pq_pred<- i_bc(par=tr_bpar,y=pq_dat0[,pq],jackup=jackup) - y_added)
            if( !inherits(ptest, "try-error") ){#class(ptest)!="try-error"
              pq_pred[is.nan(pq_pred)&(pq_dat0[,pq] < 0)]<-0
              pq_pred[pq_pred < 0 ]        <- 0
              pq_dat[,pq]       <-pq_pred
            } else {
              pq_dat[,pq]       <-NA
            }
          }

        } else {
          pq_dat   <- pq_dat0
        }
      }
    }

  }

  result <- list( pred = res, pred_quantile=pq_dat,
                  c_vc = c_vc, cse_vc =cse_vc, ct_vc = ct_vc, cp_vc = cp_vc )
  return( result )
}

