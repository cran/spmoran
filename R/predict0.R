predict0   <- function( mod, meig0, x0 = NULL, xconst0 = NULL, xgroup0 = NULL,
                        offset0 = NULL, weight0 = NULL, compute_se=FALSE,
                        compute_quantile = FALSE ){

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

  if( !is.null( mod$other$interact_sel_all ) ){

    nx            <- 1
    if( inherits( mod, "resf_vc")|inherits( mod, "resf_vc_internal" ) ){
      nx          <- ncol(mod$b_vc)
    }

    evz_num       <- x_sel         <- NULL
    if( !is.null(meig0$ev_z[[1]]) ){
      evz_num     <- c(evz_num, rep(1,nx))
      x_sel       <- c(x_sel, (1:nx) )
    }
    if( length(meig0$ev_z) == 2  ){
      evz_num     <- c(evz_num, rep(2,nx))
      x_sel       <- c(x_sel, (1:nx) )
    }

    ev_z_list<-sf_z_list  <- evSqrts_list<-list(NULL)
    ev_z_list[1:2]<-sf_z_list[1:2]  <- evSqrts_list[1:2] <-NA
    if( 1 %in% evz_num[mod$other$interact_sel_all] ){#!is.null( mod$other$interact_sel_all ) &
      ev_sel       <- mod$other$int_ev_sel_list[[1]] # ev_z_all >= ev_min_val
      ev_z_id      <- rep(1:length(meig0$ev_z[[ 1 ]]), each=length(meig0$ev))[ev_sel]
      ev_id        <- rep(1:length(meig0$ev), length(meig0$ev_z[[ 1 ]]))[ev_sel]
      ev_z_all      <- meig0$ev_z[[ 1 ]] %x% meig0$ev

      ev_z_list[[1]]<- ev_z_all[ev_sel]
      sf_z_list[[1]]<- meig0$sf[ , ev_id ]*meig0$sf_z[[ 1 ]][ , ev_z_id ]
      evSqrts_list[[1]]<-mod$other$evSqrts_t_int
    }

    if( 2 %in% evz_num[mod$other$interact_sel_all] ){#!is.null(mod$other$ev_sel_list[[2]]) &
      ev_sel       <- mod$other$int_ev_sel_list[[2]]#ev_z_all >= ev_min_val
      ev_z_id      <- rep(1:length(meig0$ev_z[[ 2 ]]), each=length(meig0$ev))[ev_sel]
      ev_id        <- rep(1:length(meig0$ev), length(meig0$ev_z[[ 2 ]]))[ev_sel]
      ev_z_all     <- meig0$ev_z[[ 2 ]] %x% meig0$ev

      ev_z_list[[2]]<- ev_z_all[ev_sel]
      sf_z_list[[2]]<- meig0$sf[ , ev_id ]*meig0$sf_z[[ 2 ]][ , ev_z_id ]
      evSqrts_list[[2]]<-mod$other$evSqrts_tc_int
    }
  }

  if( inherits( mod, "esf") ){
    n        <-length(mod$other$coords[,1])
    nx       <-mod$other$nx

    sf      	<- meig0$sf
    ev      	<- meig0$ev
    sf_id     <- rep("s",length(ev))
    if(!is.null(meig0$ev_z)){
      sf      <-cbind(sf, meig0$sf_z[[1]])
      ev      <-c(ev, meig0$ev_z[[1]])
      sf_id   <-c(sf_id, rep("z1", length(meig0$ev_z[[1]])))

      if(length(meig0$ev_z)==2){
        sf    <-cbind(sf, meig0$sf_z[[2]])
        ev    <-c(ev, meig0$ev_z[[2]])
        sf_id <-c(sf_id, rep("z2", length(meig0$ev_z[[2]])))
      }
    }
    sf        <- as.matrix( sf )

    if( is.null( dim( mod$r ) ) ){
      sf_pred	<- rep(0, dim( sf )[1] )#dim(x0)[1]
    } else {
      if( length(mod$r[,1]) == 1){
        sf_pred	<- sf[ ,mod$other$sf_id ] * c( mod$r[,1] )
      }else{
        sf_pred	<- sf[ ,mod$other$sf_id ] %*% c( mod$r[,1] )
      }
    }

    b_g0      <- NULL
    if( is.null( mod$other$x_id ) ){
      xb_pred	<- c( as.matrix( mod$b[ 1 ] ) )
      pred	  <- xb_pred + sf_pred
    } else {
      if( is.null( x0 ) ){
        message( " Note: Only spatial component (sf) is interpolated because x0 is missing")
        pred	<- sf_pred
      } else {
        x0          <- as.matrix(x0)
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
    b_vc <- bse_vc <- bt_vc <- bp_vc <- NULL
    res      <- pred
    pq_dat   <- NULL

  }

  if( inherits(mod, "resf") ){
    n        <-length(mod$other$coords[,1])
    nx       <-mod$other$nx
    meig0$sf<-meig0$sf[,1:mod$other$nev0]
    meig0$ev<-meig0$ev[ 1:mod$other$nev0]

    sf_coef    <- mod$other$b_s[[ 1 ]]
    if( length(sf_coef)==1 ){
      sf_pred	<- rep(0, dim( meig0$sf )[1] )#dim(x0)[1]
    #  if( length(mod$r[,1]) == 1){
    #    sf_pred	<- meig0$sf * c( mod$r[,1] )
    #  }else{
    #    sf_pred	<- meig0$sf %*% c( mod$r[, 1 ] )
    #  }
    #}
    } else {
      i        <-1
      evSqrts  <-mod$other$evSqrts[[ i ]]
      evSqrts_n<- NA
      basis2     <- NULL
      evSqrts2   <- NULL
      if( !is.na( evSqrts[1] ) ){
        sel        <- mod$other$omit_list[[i]]
        basis2     <- cbind( basis2, meig0$sf[,sel] )
        evSqrts2   <- c( evSqrts2, evSqrts )#[sel]
      }
      #if( !is.na( evSqrts_n[1] ) ){
      #  sel        <- mod$other$omit_list[[ mod$other$id_nsv[i] ]]
      #  basis2     <- cbind( basis2, B_n[[i]][,sel] )
      #  evSqrts2   <- c( evSqrts2, evSqrts_n )#[sel]
      #}
      if( !is.na( mod$other$evSqrts_t[[i]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntv[i] ]]
        basis2     <- cbind( basis2, meig0$sf_z[[1]][,sel] )
        evSqrts2   <- c( evSqrts2, mod$other$evSqrts_t[[i]] )#[sel]
      }
      if( !is.na( mod$other$evSqrts_tc[[i]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntcv[i] ]]
        basis2     <- cbind( basis2, meig0$sf_z[[2]][,sel] )
        evSqrts2   <- c( evSqrts2, mod$other$evSqrts_tc[[i]] )#[sel]
      }

      if( !is.na(mod$other$evSqrts_t_int[[i]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntv_interact[i] ]]
        basis2     <- cbind( basis2, sf_z_list[[1]][,sel])
        evSqrts2   <- c( evSqrts2, ev_z_list[[1]][sel])#mod$other$evSqrts_t_int[[1]]
      }

      if( !is.na(mod$other$evSqrts_tc_int[[i]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntcv_interact[i] ]]
        basis2     <- cbind( basis2, sf_z_list[[2]][,sel])
        evSqrts2   <- c( evSqrts2, ev_z_list[[2]][sel])#mod$other$evSqrts_t_int[[1]]
      }
      sf2		     <- t( t( basis2 ) * evSqrts2 )
      sf_pred    <- basis2 %*% sf_coef[-1]# sf_coef[ 1 ] +
    }
    #ne       <-length(mod$r[,1])
    #meig0$sf <-meig0$sf[,1:ne]

    n0        <- length( c(sf_pred) )

    if( !is.null( x0 )){
      xx0     <- as.matrix( x0 )
      if( is.numeric( xx0 ) == FALSE ){
        mode( xx0 ) <- "numeric"
      }
      xx0     <- xx0[,mod$other$x_id0]#xf_id0
      X0      <- as.matrix(cbind( 1, xx0)[,mod$other$x_id] )
      X1      <- NULL
      if(!is.null(mod$other$res$other$xconst)){
        X1      <- as.matrix( mod$other$res$other$xconst )[,mod$other$x_id]
      }
    } else {
      X0      <- NULL
    }

    XX_0	    <- list( NULL )
    XX	      <- NULL
    nvc       <- mod$other$res$other$nvc_xconst
    if( (sum(mod$other$res$other$nvc_xconst)>0)&( is.logical( nvc[ 1 ] ) ==FALSE )&( !is.null( x0 ) ) ){

      X1_nvc  <- as.matrix( X1 )[ , nvc ]
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
            xb_pred<- c( 1, X0) %*% mod$b[, 1 ]#[ mod$other$x_id ]
          } else {
            xb_pred<- as.matrix( cbind( 1, X0 ) ) %*% mod$b[, 1 ]#[ ,mod$other$x_id ]
          }

        } else {
          if( n0 == 1 ){
            xb_pred<- sum( X0 * c_vc ) + mod$b$Estimate[1]#[ mod$other$x_id ]
          } else {
            xb_pred<- rowSums( X0* c_vc ) + mod$b$Estimate[1]#[ ,mod$other$x_id ]
          }
        }

        pred	  <- xb_pred + sf_pred
        pred	  <- data.frame( pred = pred, xb = xb_pred, sf_residual = sf_pred )
      }
    }

    if( is.null( b_g0 ) == FALSE ){
      na_b_g0   <- is.na( b_g0 )
      if( sum( na_b_g0 ) > 0 ){
        b_g0[ na_b_g0 ] <- 0
        message( " Note: b_g = 0 is assumed for samples who does not belong to any groups in xgroup")
      }
      pred      <- data.frame( pred, b_g0 )
      pred[ ,1 ]<- pred[ ,1 ] + rowSums( b_g0 )
    }

    y0         <- mod$other$y0
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
    if( !is.na(pred0[1]) ){
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
        z_ms      <- c(mean(z0),sd(z0))
        pred00    <- pred0*z_ms[2] + z_ms[1]
        pred2     <- i_bc(par=tr_bpar,y=pred00,jackup=jackup) - y_added
        pred2[ is.nan( pred2 ) ] <- 0
        pred2[pred2 < 0 ]        <- 0
        pred_test      <- data.frame(pred=pred2, pred_transG=pred00,pred[,-1])
        pred_test$xb   <- pred_test$xb*z_ms[2]+z_ms[1]
        if(ncol(pred_test)>=3){
          pred_test[,3:ncol(pred_test)] <- pred_test[,3:ncol(pred_test)]*z_ms[2]
        }
        names(pred_test)[-c(1:2)]<-pred_name[ -1 ]
        pred      <- pred_test

      } else if( y_type=="count" ){
        pred2     <- exp( pred0 )
        if( !is.null( mod$other$offset ) ){
          if( is.null( offset0 ) ) stop( "offset0 is missing" )
          pred2  <- pred2 * offset0
        }
        pred	    <- data.frame( pred = pred2, pred_transG=pred0, pred[,-1] )
        names(pred)[-c(1:2)] <- pred_name[ -1 ]
      }
    }

    ################## pred_se
    if( mod$other$y_type == "count" & compute_se ){
      compute_se      <- FALSE
    }
    if( compute_se == FALSE & compute_quantile == TRUE ){
      compute_se      <- TRUE
    }

    if( mod$other$y_type != "count" & compute_se ){
      if( mod$other$is_weight ==TRUE ){
        if( is.null( weight0 ) ) stop( "Error: weight0 is needed to compute SE/quantile" )
      } else {
        weight0     <- NULL
      }

      if( !is.null( mod$other$x_id )&is.null( X0 ) ){
        stop( "Error: x0 is needed to compute SE" )
      }

      evSqrts    <-mod$other$evSqrts[[ 1 ]]
      if(length( evSqrts ) == 1) evSqrts <- NULL

      sel        <- mod$other$omit_list[[ 1 ]]
      XX	   <- as.matrix( cbind( 1, X0, meig0$sf[,sel] ) )

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
      if( !is.na( mod$other$res$other$evSqrts_c[[1]][1] ) ){
        for( i in 1:nnxf ){
          evSqrts_n<- mod$other$res$other$evSqrts_c[[ i ]]
          if( length( evSqrts_n ) == 1 ) evSqrts_n <- NULL ############# check !!!!!!

          if( !is.null( evSqrts_n ) ){
            XX<- cbind( XX, X0[ , i ] * B_c[[ i ]] )
          }
        }
      }

      ######### T
      if( !is.na( mod$other$evSqrts_t[[1]][1] ) ){
        sel         <- mod$other$omit_list[[ mod$other$id_ntv[1] ]]
        XX<- cbind( XX, meig0$sf_z[[1]][,sel] )
        #evSqrts2   <- c( evSqrts2, mod$other$evSqrts_t[[i]] )#[sel]
      }
      ######### Tc
      if( !is.na( mod$other$evSqrts_tc[[1]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntcv[1] ]]
        XX<- cbind( XX, meig0$sf_z[[2]][,sel] )
        #evSqrts2   <- c( evSqrts2, mod$other$evSqrts_tc[[i]] )#[sel]
      }

      ######### T_int
      if( !is.na(mod$other$evSqrts_t_int[[1]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntv_interact[1] ]]
        XX<- cbind( XX, sf_z_list[[1]][,sel] )
        #evSqrts2   <- c( evSqrts2, ev_z_list[[1]] )#[,sel]#mod$other$evSqrts_t_int[[1]]
      }

      ######### Tc_int
      if( !is.na(mod$other$evSqrts_tc_int[[1]][1] ) ){
        sel        <- mod$other$omit_list[[ mod$other$id_ntcv_interact[1] ]]
        XX<- cbind( XX, sf_z_list[[2]][,sel] )
        #evSqrts2   <- c( evSqrts2, ev_z_list[[2]] )#[,sel]#mod$other$evSqrts_t_int[[1]]
      }


      if( is.null( weight0 ) ){
        weight0  <- 1
      } else {
        weight0  <- weight0*mod$other$w_scale
      }

      X3      <- XX
      X3[,-( 1:nx )]<- t(t(XX[,-( 1:nx )])* mod$other$eevSqrt[mod$other$eevSqrt > 0])
      pred0_se<- sqrt( colSums( t( sqrt(weight0)*X3 ) * ( mod$other$B_covs %*% t( sqrt(weight0)*X3 ) ) ) + mod$other$sig )
      pred0_se<- pred0_se/sqrt( weight0 )
      if( sum(names( pred ) %in% "pred_transG") == 0 ){
        pred_name<-names( pred )
        pred     <- data.frame( pred[,1], pred_se = pred0_se, pred[,-1] )
        names( pred ) <- c( pred_name[1], "pred_se", pred_name[ -1 ] )
      } else {
        pred_name<-names( pred )
        pred     <- data.frame( pred[,1:2], pred_se = pred0_se, pred[,-(1:2)])
        names( pred ) <- c( pred_name[1:2], "pred_transG_se", pred_name[-(1:2)] )
      }
    }

    ################## pred_quantiule
    res      <- pred
    pq_dat   <- NULL
    if( compute_quantile ==TRUE ){
      if( mod$other$y_type == "count" ){
        message("Note: 'compute_quantile' is currently not supported for count data")
      } else {

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
          z_ms     <- NULL
          y_ms     <- c(mean(y),sd(y))

          pred0    <- pred0*y_ms[2] + y_ms[1]
          pred     <- i_bc(par=tr_bpar,y=pred0,jackup=jackup) - y_added
          #pred[ is.nan( pred ) ] <- 0
          #pred[pred < 0 ]        <- 0

          pq_dat   <- pq_dat0
          for(pq in 1:ncol(pq_dat0)){
            ptest<-try(pq_pred<- i_bc(par=tr_bpar,y=pq_dat0[,pq]*y_ms[2] + y_ms[1],jackup=jackup) - y_added)
            if( !inherits(ptest, "try-error") ){#class(ptest)!="try-error"
              #pq_pred[is.nan(pq_pred)&(pq_dat0[,pq] < 0)]<-0
              #pq_pred[pq_pred < 0 ]        <- 0
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

    if( !is.null(c_vc) ){
      c_vc  <-as.data.frame(c_vc)  ;names(c_vc)  <-names(mod$c_vc)
      cse_vc<-as.data.frame(cse_vc);names(cse_vc)<-names(mod$c_vc)
      ct_vc <-as.data.frame(ct_vc) ;names(ct_vc) <-names(mod$c_vc)
      cp_vc <-as.data.frame(cp_vc) ;names(cp_vc) <-names(mod$c_vc)
    }
    b_vc <- bse_vc <- bt_vc <- bp_vc <- NULL
  }

  if( inherits( mod, "resf_vc")|inherits( mod, "resf_vc_internal" ) ){
    n        <- length(mod$other$coords[,1])
    nx       <- mod$other$nx
    meig0$sf<-meig0$sf[,1:mod$other$nev0]
    meig0$ev<-meig0$ev[ 1:mod$other$nev0]

    if( ( length(mod$other$evSqrts_n)==1 )&( is.null(mod$other$evSqrt_n[[1]]))&(length(mod$other$evSqrts) > 1 )){
      mod$other$evSqrts_n <-list(NULL)
      for(pp in 2:length(mod$other$evSqrts)){
        mod$other$evSqrts_n<-c(mod$other$evSqrts_n, list(NULL))
      }
    }

    if( !is.null( mod$other$x ) ){
      mod$other$x <- as.matrix(mod$other$x)
      if( is.numeric( mod$other$x ) == FALSE ){
        mode( mod$other$x ) <- "numeric"
      }

      x_id    <- apply(mod$other$x,2,sd)!=0
      X1_0    <- mod$other$x#[,apply(mod$other$x,2,sd)!=0]
      X1      <- as.matrix( X1_0[,x_id] )
    }

    n0        <- dim(meig0$sf)[1]
    if( !is.null( x0 ) ){
      if( n0 == 1 ){
         x0	<- t(as.matrix( c(x0) ))
      } else {
         x0	<- as.matrix( x0 )
      }

      if( is.numeric( x0 ) == FALSE ){
        mode( x0 ) <- "numeric"
      }

      if(dim( as.matrix( x0 ) )[2] != length( x_id )){
        stop("x and x0 must have the same number of columns")
      }
      if( n0 == 1 ){
        X0      <- t(as.matrix(x0[ x_id]))
      } else {
        X0      <- as.matrix( x0[,x_id] )
      }
    } else {
      X0        <- NULL
    }

    if( !is.null( mod$other$xconst ) ){
      X1const <- as.matrix( as.matrix( mod$other$xconst ) [,mod$other$xf_id] )
    }


    if( !is.null( xconst0 ) ){
      if( n0 == 1 ){
        xconst0<- t(as.matrix( c(xconst0) ))
      } else {
        xconst0<- as.matrix( xconst0 )
      }

      if( dim( as.matrix( xconst0 ) )[2] != length( mod$other$xf_id ) ){
        stop("xconst and xconst0 must have the same number of columns")
      }

      if( n0 == 1 ){
        X0const <- t( as.matrix( as.matrix( xconst0 )[ mod$other$xf_id] ))
      } else {
        X0const <- as.matrix( as.matrix( xconst0 )[,mod$other$xf_id] )
      }
      if( is.numeric( X0const ) == FALSE ){
        mode( X0const ) <- "numeric"
      }

    } else {
      X0const   <- NULL
    }

    XX1_0	      <- list(NULL)
    XX1	        <- NULL
    nvc_x       <- mod$other$nvc_x
    if( is.logical( nvc_x[1] ) == FALSE ){
      if(is.null(x0)) stop( "Error: x0 is missing" )

      X1_nvc  <- as.matrix( as.matrix(X1)[ , nvc_x ] )

      if( n0 == 1 ){
        X0_nvc  <- as.matrix(X0)[ nvc_x ]
      } else {
        X0_nvc  <- as.matrix( as.matrix(X0)[ , nvc_x ] )
      }
      xxname	<- names( as.data.frame( X1 ) )[ nvc_x ]
      nnsv    <- length( xxname )
      np_xx   <-apply( X1_nvc,2,function( x ) length( unique( x )))
      np_xx   <-ifelse( np_xx < mod$other$nvc_num/0.7, round( np_xx * 0.7 ) ,mod$other$nvc_num )

      np_xx_max <-round( n/nnsv ) - 2
      np_xx[ np_xx > np_xx_max ] <-np_xx_max
      np_xx[ np_xx < 2 ] <- 2

      #XB_n    <-NULL
      B_n     <-list(NULL)
      for( ii in 1:dim( X1_nvc )[ 2 ] ){
        if( np_xx[ ii ] <= 2 ){
          B_n[[ii+1]]<- 0
          np_xx[ ii ]<- 0
        } else {

          test<-TRUE
          iiii<-0
          while( test ){
            kkk       <- np_xx[ ii ]-iiii
            knots     <- seq(min( X1_nvc[ , ii ] ),max( X1_nvc[ , ii ] ),len=kkk+2)[2:(kkk+1)]
            testt     <- try(XX1_00    <- ns( X1_nvc[ ,ii ], knots = knots ), silent=TRUE)#replaced
            test      <- inherits( testt, "try-error")#class(testt)[1] == "try-error"
            iiii      <- iiii+1
          }

          XX1_0       <- cbind( X1_nvc[,ii], XX1_00)
          if( n0 == 1 ){
            XX0_0     <- cbind( X0_nvc[ ii], predict( XX1_00, newx = X0_nvc[ ii ] ) )
          } else {
            XX0_0     <- cbind( X0_nvc[,ii], predict( XX1_00, newx = X0_nvc[,ii ] ) )
          }

          #if( !is.na( mod$other$sel_basis_n[[ ii ]][ 1 ] ) ){
          #}

          if( n0 == 1 ){
            n_ids <-(1:length(c(XX0_0)))[(1:length(c(XX0_0))) %in% mod$other$sel_basis_n[[ii]]]
            for( j in length(c(XX0_0)) ){
              XX0_0[ j]<- ( XX0_0[ j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])### check
            }
          } else {
            n_ids <-(1:dim(XX0_0)[2])[(1:dim(XX0_0)[2]) %in% mod$other$sel_basis_n[[ii]]]
            for( j in  n_ids){
              XX0_0[,j]<- ( XX0_0[,j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])### check
            }
          }

          if( n0 == 1 ){
            B_n[[ii+1]]<- XX0_0[ mod$other$sel_basis_n[[ii]] ]
          } else {
            B_n[[ii+1]]<- XX0_0[ ,mod$other$sel_basis_n[[ii]] ]
          }

        }
      }
    }

    XXconst_0	  <- list( NULL )
    XXconst	    <- NULL
    nvc_xconst <- mod$other$nvc_xconst
    if( is.logical( nvc_xconst[ 1 ] ) ==FALSE ){
      if(is.null(X0const)) stop( "Error: xconst0 is missing" )

      X1const_nvc<- as.matrix( X1const[ , nvc_xconst ] )
      if( n0 == 1 ){
        X0const_nvc<- c( X0const[ nvc_xconst ] )
      } else {
        X0const_nvc<- as.matrix( X0const[ , nvc_xconst ] )
      }

      xxfname	  <- names( as.data.frame( X1const ) )[ nvc_xconst ]
      nnxf      <- length( xxfname )
      np_xxconst<-apply( X1const_nvc,2,function( x ) length( unique( x )))
      np_xxconst<-ifelse( np_xxconst < mod$other$nvc_num/0.7, round( np_xxconst * 0.7 ) ,mod$other$nvc_num )
      np_xxconst_max <-round( n/nnxf ) - 2
      np_xxconst[ np_xxconst > np_xxconst_max ] <-np_xxconst_max
      np_xxconst[ np_xxconst < 2 ] <- 2

      #XB_c    <-NULL
      B_c     <-list(NULL)
      for( ii in 1:dim( X1const_nvc )[ 2 ] ){
        if( np_xxconst[ ii ] <= 2 ){
          B_c[[ii+1]]     <- 0
          np_xxconst[[ii]]<- 0
        } else {

          test<-TRUE
          iiii<-0
          while(test){
            kkk      <- np_xxconst[ ii ]-iiii
            knots    <-seq(min( X1const_nvc[ ,ii ] ),max( X1const_nvc[ ,ii ] ),len=kkk+2)[2:(kkk+1)]
            testt<- try(XX1const_00<- ns( X1const_nvc[ , ii], knots = knots ), silent=TRUE)
            test <- inherits( testt, "try-error")#class(testt)[1] == "try-error"
            iiii <- iiii+1
          }

          XX1const_0   <- cbind( X1const[,ii] , XX1const_00)
          if( n0 == 1 ){
            XX0const_0 <- cbind( X0const_nvc[ ii], predict( XX1const_00, newx= X0const_nvc[ ii]) )
          } else {
            XX0const_0 <- cbind( X0const_nvc[,ii], predict( XX1const_00, newx= X0const_nvc[,ii]) )
          }

          if( !is.na( mod$other$sel_basis_c[[ ii ]][ 1 ] ) ){
            XX1const_0<- XX1const_0[,mod$other$sel_basis_c[[ii]]]
            if( n0 == 1){
              XX0const_0<- XX0const_0[ mod$other$sel_basis_c[[ii]]]
            } else {
              XX0const_0<- XX0const_0[,mod$other$sel_basis_c[[ii]]]
            }
          }

          if( n0 == 1 ){
            for( j in 1:length(XX0const_0)){
              XX0const_0[ j]<- ( XX0const_0[ j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])### check
            }
          } else {
            for( j in 1:dim(XX0const_0)[2]){
              XX0const_0[,j]<- ( XX0const_0[,j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])### check
            }
          }

          B_c[[ii]]    <- XX0const_0
          #XB_c         <- cbind( XB_c , X0const_nvc[, ii ] * B_c[[ii]] )
        }
      }
    }

    n	  <- length( mod$other$y )
    n0	<- length( meig0$sf[ , 1 ] )
    #ne	<- length( mod$other$b_s[[ 1 ]][ -1 ] )
    nsv	<- sum( x_id ) + 1#mod$other$x_id
    nxf	<- sum( mod$other$xf_id )
    meig0$sf<- as.matrix( meig0$sf)#[ , 1:ne ]

    #for( i in 1:length(mod$other$evSqrts_n)){
    #  if(length( mod$other$evSqrts_n[[i]] ) == 1) mod$other$evSqrts_n[[i]] <- NULL
    #}
    #for( i in 1:length(mod$other$evSqrts)){
    #  if(length( mod$other$evSqrts[[i]] ) == 1) mod$other$evSqrts[[i]] <- NULL
    #}

    b_vc	<- matrix(0, nrow = n0, ncol = nsv )
    bse_vc<- matrix(0, nrow = n0, ncol = nsv )
    bt_vc	<- matrix(0, nrow = n0, ncol = nsv )
    bp_vc	<- matrix(0, nrow = n0, ncol = nsv )
    for( i in 1:nsv ){
      evSqrts    <-mod$other$evSqrts[[ i ]]
      if( i == 1 ){
        evSqrts_n<- NA
      } else {
        evSqrts_n<- mod$other$evSqrts_n[[ i ]]
      }

      if( n0 == 1 ){
        b_vc_coef  <- mod$other$b_s[[ i ]]
        b_vc[ , i ]<- b_vc_coef[ 1 ]# + basis2 %*% b_vc_coef[-1]
        basis2     <- evSqrts2   <- NULL
        if( length(b_vc_coef)>1 ){
          if( !is.na( evSqrts[1] ) ){
            sel        <- mod$other$omit_list[[i]]
            basis2     <- c( basis2, meig0$sf[,sel] )
            evSqrts2   <- c( evSqrts2, evSqrts )#[sel]
          }
          if( !is.na( evSqrts_n[1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_nsv[i] ]]
            basis2     <- c( basis2, B_n[[i]][,sel] )
            evSqrts2   <- c( evSqrts2, evSqrts_n )#[sel]
          }
          if( !is.na( mod$other$evSqrts_t[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntv[i] ]]
            basis2     <- c( basis2, meig0$sf_z[[1]][,sel] )
            evSqrts2   <- c( evSqrts2, mod$other$evSqrts_t[[i]] )#[sel]
          }
          if( !is.na( mod$other$evSqrts_tc[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntcv[i] ]]
            basis2     <- c( basis2, meig0$sf_z[[2]][,sel] )
            evSqrts2   <- c( evSqrts2, mod$other$evSqrts_tc[[i]] )#[sel]
          }

          if( !is.na(mod$other$evSqrts_t_int[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntv_interact[i] ]]
            basis2     <- c( basis2, sf_z_list[[1]][,sel] )
            evSqrts2   <- c( evSqrts2, ev_z_list[[1]][sel] )#[,sel]#mod$other$evSqrts_t_int[[1]]
          }

          if( !is.na(mod$other$evSqrts_tc_int[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntcv_interact[i] ]]
            basis2     <- c( basis2, sf_z_list[[2]][,sel])
            evSqrts2   <- c( evSqrts2, ev_z_list[[2]][sel])#mod$other$evSqrts_t_int[[1]]
          }
          b_vc[ , i ]<- b_vc[ , i ] + sum( basis2 * b_vc_coef[-1])

          x_sf		        <- t( as.matrix( c( 1, basis2 * evSqrts2 ) ))#as.matrix( cbind( 1, sf2 ) )
          bse_vc[ , i ]	  <- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs[[ i ]] %*% t( x_sf ) ) ) )
        } else {
          bse_vc[ , i ]	  <- sqrt( mod$other$b_covs[[ i ]] )
        }

      } else {

        b_vc_coef  <- mod$other$b_s[[ i ]]
        b_vc[ , i ]<- b_vc_coef[ 1 ]# + basis2 %*% b_vc_coef[-1]
        basis2     <- evSqrts2   <- NULL
        if( length(b_vc_coef)>1){
          if( !is.na( evSqrts[1] ) ){
            sel        <- mod$other$omit_list[[i]]
            basis2     <- cbind( basis2, meig0$sf[,sel] )
            evSqrts2   <- c( evSqrts2, evSqrts )#[sel]
          }
          if( !is.na( evSqrts_n[1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_nsv[i] ]]
            basis2     <- cbind( basis2, B_n[[i]][,sel] )
            evSqrts2   <- c( evSqrts2, evSqrts_n )#[sel]
          }
          if( !is.na( mod$other$evSqrts_t[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntv[i] ]]
            basis2     <- cbind( basis2, meig0$sf_z[[1]][,sel] )
            evSqrts2   <- c( evSqrts2, mod$other$evSqrts_t[[i]] )#[sel]
          }
          if( !is.na( mod$other$evSqrts_tc[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntcv[i] ]]
            basis2     <- cbind( basis2, meig0$sf_z[[2]][,sel] )
            evSqrts2   <- c( evSqrts2, mod$other$evSqrts_tc[[i]] )#[sel]
          }

          if( !is.na(mod$other$evSqrts_t_int[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntv_interact[i] ]]
            basis2     <- cbind( basis2, sf_z_list[[1]][,sel] )
            evSqrts2   <- c( evSqrts2, ev_z_list[[1]][sel] )#[,sel]#mod$other$evSqrts_t_int[[1]]
          }

          if( !is.na(mod$other$evSqrts_tc_int[[i]][1] ) ){
            sel        <- mod$other$omit_list[[ mod$other$id_ntcv_interact[i] ]]
            basis2     <- cbind( basis2, sf_z_list[[2]][,sel])
            evSqrts2   <- c( evSqrts2, ev_z_list[[2]][sel])#mod$other$evSqrts_t_int[[1]]
          }
          b_vc[ , i ]<- b_vc[ , i ] + basis2 %*% b_vc_coef[-1]

          x_sf		        <- as.matrix( cbind( 1, t( t( basis2 ) * evSqrts2 ) ) )
          bse_vc[ , i ]	  <- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs[[ i ]] %*% t( x_sf ) ) ) )
        } else {
          bse_vc[ , i ]	  <- sqrt( mod$other$b_covs[[ i ]] )
        }
      }

      if( i ==1 ) b_vc[ , i ] <- b_vc[ , i ] + mod$other$Bias      #### added

      bt_vc[ , i ]	  <- b_vc[ , i ] / bse_vc[ , i ]
      bp_vc[ , i ]	  <- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    }

    if( is.null( mod$c_vc ) ){
      c_vc  <- cse_vc <- ct_vc <- cp_vc <- NULL
    } else {
      xc_vc	<- 0
      c_vc	<- matrix(0, nrow = n0, ncol = nnxf )
      cse_vc<- matrix(0, nrow = n0, ncol = nnxf )
      ct_vc	<- matrix(0, nrow = n0, ncol = nnxf )
      cp_vc	<- matrix(0, nrow = n0, ncol = nnxf )
      for( i in 1:nnxf ){
        evSqrts_c<- mod$other$evSqrts_c[[ i ]]
        if(length( evSqrts_c ) == 1) evSqrts_c <- NULL

        if( length( mod$other$evSqrts_c[[ i ]] ) <= 1 ){#!= ne
          c_vc[ , i ]	  <- mod$other$b_c[[ i ]][ 1 ]
          cse_vc[ , i ]	<- sqrt( mod$other$b_covs_c[[ i ]] )
          ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
          cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$df )
        } else {
          if( n0 == 1 ){
            c_vc[ , i ]	<- mod$other$b_c[[ i ]][1] + c(B_c[[ i ]]) %*% c( mod$other$b_c[[ i ]][ -1 ])
            sf2		<- t( B_c[[ i ]] * mod$other$evSqrts_c[[ i ]] )
          } else {
            c_vc[ , i ]	<- mod$other$b_c[[ i ]][1] + B_c[[ i ]] %*% c( mod$other$b_c[[ i ]][ -1 ] )
            sf2		<- t( t( B_c[[ i ]] ) * mod$other$evSqrts_c[[ i ]] )
          }
          x_sf		      <- as.matrix( cbind( 1, sf2 ) )
          cse_vc[ , i ]	<- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs_c[[ i ]] %*% t( x_sf ) ) ) )
          ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
          cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$df )
        }
      }
    }


    if( !is.null( mod$b_g ) ){
      if( is.null(xgroup0) ){
        message( "Note: Group effects are ignored because xgroup0 is missing")
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
    b_g    <- b_g0

    xb_vc     <- 0
    if( !is.null( x0 )&(sum(mod$other$x_id)>0) ){
      x_id    <- apply(mod$other$x,2,sd)!=0

      if( length( mod$other$x_id ) == 1 ){
        x0	  <- as.matrix( cbind( 1, x0[  x_id ] ))
      } else {
        if(n0==1){
          x0	<- t( as.matrix( c( 1, x0[, x_id ] )))
        } else {
          x0	<- as.matrix( cbind( 1, x0[, x_id ] ))
        }
      }
      if( is.numeric( x0 ) == FALSE ){
        mode( x0 )	<- "numeric"
      }

      if( n0 == 1 ){
        xb_vc   <- sum(x0[-1]*b_vc[-1])
      } else {
        xb_vc   <- rowSums(as.matrix(x0[,-1]*b_vc[,-1]))
      }
    }

    xb_const  <- 0
    if( !is.null( mod$other$xf_id ) ){
      if( !is.null( xconst0 ) ){
        if( length( mod$other$xf_id ) == 1 ){
          xconst0	<- as.matrix( xconst0[  mod$other$xf_id ] )
        } else {
          if(n0==1){
            xconst0	<- t( as.matrix( c( xconst0[, mod$other$xf_id ] )))
          } else {
            xconst0	<- as.matrix( xconst0[, mod$other$xf_id ] )
          }
        }
        if( is.numeric( xconst0 ) == FALSE ){
          mode( xconst0 )	<- "numeric"
        }
      }
    }

    if( !is.null( xconst0 ) ){
      if( !is.null( mod$c_vc ) ){
        if( n0 == 1 ){
          xb_const  <- sum( xconst0 * c_vc )
        } else {
          xb_const  <- rowSums( xconst0 * c_vc )
        }
      } else {
        if( !is.null( mod$c ) ){
          if( is.null( dim( mod$c ) ) ){
            xb_const<- xconst0 * mod$c[  1 ]
          } else {
            xb_const<- xconst0 %*% mod$c[, 1 ]
          }
        } else {
          xb_const  <- 0
        }
      }
    }

    if( n0 == 1 ){
      sf_pred <- b_vc[ 1]
    } else {
      sf_pred <- b_vc[,1]
    }

    x_problem      <- (length(mod$other$x_id)>0 ) & is.null( x0 )
    xconst_problem <- ( length(mod$other$xf_id) > 0 ) & is.null( xconst0 )
    #xconst_problem <- !is.null( mod$other$xf_id ) & is.null( xconst0 )
    if( ( x_problem == FALSE ) & ( xconst_problem == FALSE ) ){
      sf_pred	<- sf_pred - mean( mod$b_vc[,1] )
      xb	    <- xb_vc + xb_const + mean( mod$b_vc[,1] )
      pred0	<- xb + sf_pred
      if( is.null( b_g0 ) == FALSE ){
        na_b_g0   <- is.na( b_g0 )
        if( sum( na_b_g0 ) > 0 ){
          b_g0[ na_b_g0 ] <- 0
          message( "Note: NAs are given to the groups that are not in xgroup")
        }
        pred0   <- pred0 + rowSums( b_g0 )
      }

      y0         <- mod$other$y0
      tr_num     <- mod$other$tr_num
      y_nonneg   <- mod$other$y_nonneg
      y_type     <- mod$other$y_type
      tr_par     <- mod$tr_par
      tr_bpar    <- mod$tr_bpar$Estimate
      y_added    <- mod$other$y_added
      jackup     <- mod$other$jackup
      noconst_last<-TRUE
      if( tr_num > 0 ){
        z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,jackup=jackup)
        z_ms     <- z0$z_ms
        y_ms     <- z0$y_ms

        dif      <- 1/d_sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)
        pred     <- i_sal_k(par=tr_par,y=pred0,k=tr_num,noconst_last=noconst_last,
                            bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup) - y_added
        if( y_nonneg ) pred[ pred < 0 ] <- 0
        if( y_type=="count" ){
          pred <- exp( pred )
          if( !is.null( mod$other$offset ) ){
            if( is.null( offset0 ) ) stop( "offset0 is missing" )
            pred <- pred * offset0
          }
        }

      } else if( y_nonneg ){
        z0      <- bc(par=tr_bpar,y=y0,jackup=jackup)
        z_ms    <- c(mean(z0),sd(z0))
        y_ms    <- NULL

        dif      <- 1 / ( d_bc(par=tr_bpar,y=y0,jackup=jackup)*d_sc(par=z_ms,y) )#/z_ms[2]
        pred00   <- pred0*z_ms[2] + z_ms[1]
        pred     <- i_bc(par=tr_bpar,y=pred00,jackup=jackup) - y_added
        pred[is.nan(pred)]<-0
        pred[ pred < 0 ] <- 0

      } else {
        dif      <- 1
        z_ms     <- NULL
        y_ms     <- NULL
        pred     <- pred0
        if( y_type=="count" ){
          pred <- exp( pred )
          if( !is.null( mod$other$offset ) ){
            if( is.null( offset0 ) ) stop( "offset0 is missing" )
            pred <- pred * offset0
          }
        }
      }

      if( is.null( z_ms )&( y_nonneg==FALSE ) ){
        if( !is.null( b_g ) ){
          res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred, b_g )
        } else {
          res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred )
        }

      } else {
        if( !is.null( b_g ) ){
          res	<- data.frame( "pred" = pred, "pred_transG" = pred0, "xb" = xb, "sf_residual" = sf_pred, b_g )
        } else {
          res	<- data.frame( "pred" = pred, "pred_transG" = pred0, "xb" = xb, "sf_residual" = sf_pred )
        }
        if( y_nonneg & tr_num==0 ) res$pred_transG <- pred00
      }

    } else {
      res	<- NULL
      message( "Note: y is not predicted because x0 and/or xconst0 is missing" )
    }

    b_vc	  <- as.data.frame( b_vc )
    bse_vc	<- as.data.frame( bse_vc )
    bt_vc	  <- as.data.frame( bt_vc )
    bp_vc	  <- as.data.frame( bp_vc )
    names( b_vc )	 <- names( mod$b_vc )
    names( bse_vc )<- names( mod$b_vc )
    names( bt_vc ) <- names( mod$b_vc )
    names( bp_vc ) <- names( mod$b_vc )

    c_vc	  <- as.data.frame( c_vc )
    cse_vc	<- as.data.frame( cse_vc )
    ct_vc	  <- as.data.frame( ct_vc )
    cp_vc	  <- as.data.frame( cp_vc )

    if( dim(c_vc)[1] > 0){
      names( c_vc )	 <- names( mod$c_vc )
      names( cse_vc )<- names( mod$c_vc )
      names( ct_vc ) <- names( mod$c_vc )
      names( cp_vc ) <- names( mod$c_vc )
    } else {
      c_vc <- cse_vc <- ct_vc <- cp_vc <- NULL
    }

    ################## standard error
    if( mod$other$y_type == "count" & compute_se == TRUE ){
      compute_se      <- FALSE
    }
    if( compute_se == FALSE & compute_quantile == TRUE ){
      compute_se      <- TRUE
    }

    if( mod$other$y_type != "count" & compute_se & !is.null(res) ){
      if( mod$other$is_weight ==TRUE ){
        if( is.null( weight0 ) ) stop( "Error: weight0 is needed to compute SE/quantile" )
      } else {
        weight0     <- NULL
      }

      XX	<- as.matrix( cbind( rep(1,n0), X0const, X0 ) )

      #########SVC
      for( i in 1:nsv ){
        evSqrts    <-mod$other$evSqrts[[ i ]]
        if(length( evSqrts ) == 1) evSqrts <- NULL

        sel        <- mod$other$omit_list[[ i ]]
        if( !is.null( evSqrts ) ){
          if( i == 1 ) {
            if(n0==1){
              XX<- t(as.matrix( c(XX, meig0$sf[,sel] )))
            } else {
              XX<- cbind( XX, meig0$sf[,sel] )
            }

          } else {
            if(n0==1){
              XX<- t(as.matrix( c(XX, X0[ , i-1 ]*meig0$sf[,sel] )))
            } else {
              XX<- cbind( XX, X0[ , i-1 ] * meig0$sf[,sel] )
            }
          }
        }
      }

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

      #########NVConst
      if( !is.null( mod$c_vc ) ){
        for( i in 1:nnxf ){
          evSqrts_c<- mod$other$evSqrts_c[[ i ]]
          if(length( evSqrts_c ) == 1) evSqrts_c <- NULL

          if( length( mod$other$evSqrts_c[[ i ]] ) <= 1 ){
          } else {
            XX<- cbind( XX, X0const[ , i ] * B_c[[ i ]] )
          }
        }
      }

      #########NVC
      if( !is.null( mod$other$evSqrts_n[[1]] ) ){
        for( i in 1:nsv ){
          if(i ==1){
            evSqrts_n<- NULL
          } else {
            evSqrts_n<- mod$other$evSqrts_n[[ i ]]
            if(length( evSqrts_n ) == 1) evSqrts_n <- NULL
          }

          if( !is.null( evSqrts_n ) ){
            XX<- cbind( XX, X0[ , i-1 ] * B_n[[ i ]] )
          }
        }
      }

      ######### T
      for( i in 1:nsv ){
        if( !is.na( mod$other$evSqrts_t[[i]][1] ) ){
          sel        <- mod$other$omit_list[[ mod$other$id_ntv[i] ]]
          if( i == 1 ) {
            XX<- cbind( XX, meig0$sf_z[[1]][,sel] )
          } else {
            XX<- cbind( XX, X0[ , i-1 ] * meig0$sf_z[[1]][,sel] )
          }
          #evSqrts2   <- c( evSqrts2, mod$other$evSqrts_t[[i]] )#[sel]
        }
      }
      ######### Tc
      for( i in 1:nsv ){
        if( !is.na( mod$other$evSqrts_tc[[i]][1] ) ){
          sel        <- mod$other$omit_list[[ mod$other$id_ntcv[i] ]]
          if( i == 1 ) {
            XX<- cbind( XX, meig0$sf_z[[2]][,sel] )
          } else {
            XX<- cbind( XX, X0[ , i-1 ] * meig0$sf_z[[2]][,sel] )
          }
          #evSqrts2   <- c( evSqrts2, mod$other$evSqrts_tc[[i]] )#[sel]
        }
      }

      ######### T_int
      for( i in 1:nsv ){
        if( !is.na(mod$other$evSqrts_t_int[[i]][1] ) ){
          sel        <- mod$other$omit_list[[ mod$other$id_ntv_interact[i] ]]

          if( i == 1 ) {
            XX<- cbind( XX, sf_z_list[[1]][,sel] )
          } else {
            XX<- cbind( XX, X0[ , i-1 ] * sf_z_list[[1]][,sel] )
          }
          #evSqrts2   <- c( evSqrts2, ev_z_list[[1]] )#[,sel]#mod$other$evSqrts_t_int[[1]]
        }
      }

      ######### Tc_int
      for( i in 1:nsv ){
        if( !is.na(mod$other$evSqrts_tc_int[[i]][1] ) ){
          sel        <- mod$other$omit_list[[ mod$other$id_ntcv_interact[i] ]]

          if( i == 1 ) {
            XX<- cbind( XX, sf_z_list[[2]][,sel] )
          } else {
            XX<- cbind( XX, X0[ , i-1 ] * sf_z_list[[2]][,sel] )
          }
          #evSqrts2   <- c( evSqrts2, ev_z_list[[2]] )#[,sel]#mod$other$evSqrts_t_int[[1]]
        }
      }

      if( is.null( weight0 ) ){
        weight0  <- 1
      } else {
        weight0  <- c( weight0*mod$other$w_scale )
      }

      X3            <- XX
      X3[,-( 1:nx )]<- t(t(XX[,-( 1:nx )])* mod$other$eevSqrt[mod$other$eevSqrt > 0])

      B_covs  <- mod$other$B_covs
      sig     <- mod$other$sig
      pred0_se<- sqrt( colSums( t( sqrt(weight0)*X3 ) * ( B_covs %*% t( sqrt(weight0)*X3 ) ) ) + sig )
      pred0_se<- pred0_se/sqrt( weight0 )
      if( sum(names(res) %in% "pred_transG") == 0 ){
        res_name<-names( res )
        res   <- data.frame( res[,1], pred_se = pred0_se, res[,-1] )
        names( res ) <- c( res_name[1], "pred_se", res_name[-1] )
      } else {
        res_name<-names( res )
        res   <- data.frame( res[,1:2], pred_se = pred0_se, res[,-c(1:2)] )
        names( res ) <- c( res_name[1:2], "pred_transG_se", res_name[-c(1:2)] )
      }
    }

    ################## quantile
    pq_dat           <- NULL
    if( compute_quantile ==TRUE ){
      if( mod$other$y_type == "count" ){
        message("Note: 'compute_quantile' is currently not supported for count data")
      } else {
        if( is.null( res ) ){
          stop( "x0 and/or xconst0 are missing. They are required to compute quantile" )
        }

        if( mod$other$is_weight ==TRUE ){
          if( is.null( weight0 ) ) stop( "Specify weight0 to compute quantile" )
        } else {
          weight0     <- NULL
        }

        pquant  <- c(0.01, 0.025, 0.05, seq(0.1,0.9,0.1), 0.95, 0.975, 0.99)
        pq_dat0 <- NULL
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
            if( !inherits(ptest, "try-error")){#class(ptest)!="try-error"
              pq_dat[,pq]       <-pq_pred
            } else {
              pq_dat[,pq]       <-NA
            }
          }

        } else if( y_nonneg ==TRUE ){
          z0      <- bc(par=tr_bpar,y=y0,jackup=jackup)
          z_ms    <- c(mean(z0),sd(z0))
          y_ms    <- NULL

          #dif      <- 1 / ( d_bc(par=tr_bpar,y=y0,jackup=jackup)*d_sc(par=z_ms,y) )#/z_ms[2]
          pred0    <- pred0*z_ms[2] + z_ms[1]
          pred     <- i_bc(par=tr_bpar,y=pred0,jackup=jackup) - y_added
          pred[is.nan(pred)]<-0
          pred[ pred < 0 ] <- 0

          #y        <- bc(par=tr_bpar,y=y0,jackup=jackup)
          #pred     <- i_bc(par=tr_bpar,y=pred0,jackup=jackup) - y_added
          #pred[ is.nan( pred ) ]<- 0
          #pred[ pred < 0 ]      <- 0

          pq_dat   <- pq_dat0
          for(pq in 1:ncol(pq_dat0)){
            ptest<-try(pq_pred<- i_bc(par=tr_bpar,y=pq_dat0[,pq],jackup=jackup) - y_added)
            if( !inherits(ptest, "try-error") ){#class(ptest)!="try-error"
              pq_pred[is.nan(pq_pred)&(pq_dat0[,pq] < 0)]<-0
              pq_pred[pq_pred <0]<-0
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

  result	<- list( pred = res, pred_quantile=pq_dat,#XX=XX,
                  b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc,
                  c_vc = c_vc, cse_vc = cse_vc, ct_vc = ct_vc, cp_vc = cp_vc )
  return( result )
}
