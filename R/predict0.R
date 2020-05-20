predict0   <- function( mod, meig0, x0 = NULL, xgroup0 = NULL ){

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
    X0 <- as.matrix( x0 )[,mod$other$res$other$xf_id]
   } else {
     X0<-NULL
   }

   XX_0	  <- list( NULL )
   XX	    <- NULL
   nvc <- mod$other$res$other$nvc_xconst
   if( ( is.logical( nvc[ 1 ] ) ==FALSE )&( !is.null( x0 ) ) ){

     X1_nvc<- as.matrix( X1 )[ , nvc ]
     if( n0 == 1 ){
      X0_nvc<- as.matrix( X0 )[ nvc ]
    } else {
      X0_nvc<- as.matrix( X0 )[ , nvc ]
    }

    xxfname	  <- names( as.data.frame( X1 ) )[ nvc ]
    nnxf      <- length( xxfname )
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

   	     pred	  <- xb_pred + sf_pred
   	     pred	  <- as.data.frame( cbind( pred, xb_pred, sf_pred ) )
   	     names( pred )<- c( "pred", "xb", "sf_residual" )
   	   } else {
   	     if( n0 == 1 ){
   	       xb_pred<- sum( X0[ mod$other$x_id ] * c_vc ) + mod$b$Estimate[1]
   	     } else {
   	       xb_pred<- rowSums( X0[ ,mod$other$x_id ] * c_vc ) + mod$b$Estimate[1]
   	     }
   	     pred	  <- xb_pred + sf_pred
   	     pred	  <- as.data.frame( cbind( pred, xb_pred, sf_pred ) )
   	     names( pred )<- c( "pred", "xb", "sf_residual" )
   	   }
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
  }
  return( list( pred = pred, c_vc = c_vc, cse_vc =cse_vc, ct_vc = ct_vc, cp_vc = cp_vc ) )
}

