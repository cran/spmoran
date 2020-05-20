predict0_vc	<- function( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL ){

  n<-length(mod$other$coords[,1])
  if( ( length(mod$other$evSqrts_n)==1 )&( is.null(mod$other$evSqrt_n[[1]]))&(length(mod$other$evSqrts) > 1 )){
    mod$other$evSqrts_n <-list(NULL)
    for(pp in 2:length(mod$other$evSqrts)){
      mod$other$evSqrts_n<-c(mod$other$evSqrts_n, list(NULL))
    }
  }


  if( !is.null( mod$other$x ) ){
    X1      <- as.matrix( as.matrix( mod$other$x ) [,mod$other$x_id] )
  }

   n0        <- dim(meig0$sf)[1]
   if( !is.null( x0 ) ){
     if(dim( as.matrix( x0 ) )[2] != length( mod$other$x_id )){
       stop("x and x0 must have the same number of columns")
     }
     if( n0 == 1){
       X0      <- as.matrix( x0 )[ mod$other$x_id]
     } else {
       X0      <- as.matrix( as.matrix( x0 )[,mod$other$x_id] )
     }
   }

   if( !is.null( mod$other$xconst ) ){
     X1const <- as.matrix( as.matrix( mod$other$xconst ) [,mod$other$xf_id] )
   }


   if( !is.null( xconst0 ) ){
     if( dim( as.matrix( xconst0 ) )[2] != length( mod$other$xf_id ) ){
       stop("xconst and xconst0 must have the same number of columns")
     }

     if( n0 == 1 ){
       X0const <- as.matrix( xconst0 )[ mod$other$xf_id]
     } else {
       X0const <- as.matrix( as.matrix( xconst0 )[,mod$other$xf_id] )
     }
   }

  XX1_0	      <- list(NULL)
  XX1	        <- NULL
  nvc_x     <- mod$other$nvc_x
  if( is.logical( nvc_x[1] ) == FALSE ){
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
          test      <- class(testt)[1] == "try-error"
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
            XX0_0[ j]<- ( XX0_0[ j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])
          }
        } else {
          n_ids <-(1:dim(XX0_0)[2])[(1:dim(XX0_0)[2]) %in% mod$other$sel_basis_n[[ii]]]
          for( j in  n_ids){
            XX0_0[,j]<- ( XX0_0[,j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])
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
         test <- class(testt)[1] == "try-error"
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
           XX0const_0[ j]<- ( XX0const_0[ j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])
         }
       } else {
         for( j in 1:dim(XX0const_0)[2]){
           XX0const_0[,j]<- ( XX0const_0[,j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])
         }
       }

       B_c[[ii]]    <- XX0const_0
       #XB_c         <- cbind( XB_c , X0const_nvc[, ii ] * B_c[[ii]] )
     }
    }
  }

  n	  <- length( mod$pred )
  n0	<- length( meig0$sf[ , 1 ] )
  ne	<- length( mod$other$b_s[[ 1 ]][ -1 ] )
  nsv	<- sum( mod$other$x_id ) + 1
  nxf	<- sum( mod$other$xf_id )
  meig0$sf<- as.matrix( meig0$sf[ , 1:ne ])

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
    if(length( evSqrts ) == 1) evSqrts <- NULL

    if(i ==1){
      evSqrts_n<- NULL
    } else {
      evSqrts_n<- mod$other$evSqrts_n[[ i ]]
      if(length( evSqrts_n ) == 1) evSqrts_n <- NULL
    }

    evSqrt_i        <-c(evSqrts, evSqrts_n)
    if( length( evSqrt_i ) <= 1 ){
      b_vc[ , i ]	  <- mod$other$b_s[[ i ]][ 1 ]
      if( i ==1 ) b_vc[ , i ]<-b_vc[ , i ] + mod$other$Bias      #### added

      bse_vc[ , i ]	<- sqrt( mod$other$b_covs[[ i ]] )
      bt_vc[ , i ]	<- b_vc[ , i ] / bse_vc[ , i ]
      bp_vc[ , i ]	<- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    } else {
      if( n0 == 1 ){
        if( is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + c(meig0$sf) %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( meig0$sf * evSqrts )
        } else if( !is.null( evSqrts_n ) & is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + c(B_n[[i]]) %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( B_n[[i]] * evSqrts_n )
        } else if( !is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          basis2     <- c( meig0$sf, B_n[[i]] )
          evSvqr2    <- c( evSqrts, evSqrts_n )
          sf2		     <- t( basis2 * evSvqr2 )
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + basis2 %*% c(  mod$other$b_s[[ i ]][ -1 ] )
        }

      } else {
        if( is.null( evSqrts_n[[ i ]] ) & !is.null( evSqrts[[ i ]] ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + meig0$sf %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( meig0$sf ) * evSqrts[[ i ]] )
        } else if( !is.null( evSqrts_n[[ i ]] ) & is.null( evSqrts[[ i ]] ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + B_n[[i]] %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( B_n[[i]] ) * evSqrts_n[[ i ]] )
        } else if( !is.null( evSqrts_n[[ i ]] ) & !is.null( evSqrts[[ i ]] ) ){
          basis2     <- cbind( meig0$sf, B_n[[i]] )
          evSqrts2    <- c( evSqrts[[ i ]], evSqrts_n[[ i ]] )
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + basis2 %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( basis2 ) * evSqrts2 )
        }
      }

      if( i ==1 ) b_vc[ , i ] <- b_vc[ , i ] + mod$other$Bias      #### added

      x_sf		        <- as.matrix( cbind( 1, sf2 ) )
      bse_vc[ , i ]	  <- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs[[ i ]] %*% t( x_sf ) ) ) )
      bt_vc[ , i ]	  <- b_vc[ , i ] / bse_vc[ , i ]
      bp_vc[ , i ]	  <- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    }
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


  xb_vc	<- 0
  if( is.null( x0 ) == FALSE ){
    if( length( mod$other$x_id ) == 1 ){
      x0	<- as.matrix( cbind( 1, x0[  mod$other$x_id ] ))
    } else {
      x0	<- as.matrix( cbind( 1, x0[, mod$other$x_id ] ))
    }
    if( is.numeric( x0 ) == FALSE ){
      mode( x0 )	<- "numeric"
    }
  }

  if( is.null( mod$other$xf_id ) ){
    xb_const	<- 0
    if( !is.null( xconst0 ) ) message( "Note: xconst0 is ignored" )
  } else {
    if( is.null( xconst0 ) ){
      xb_const	<- 0
    } else {
      if( length( mod$other$xf_id ) == 1 ){
        xconst0	<- as.matrix( xconst0[  mod$other$xf_id ] )
      } else {
        xconst0	<- as.matrix( xconst0[, mod$other$xf_id ] )
      }
      if( is.numeric( xconst0 ) == FALSE ){
        mode( xconst0 )	<- "numeric"
      }
    }
  }

  if( n0 == 1 ){
    xb_vc   <- sum(x0[-1]*b_vc[-1])
    sf_pred <- b_vc[ 1]
  } else {
    xb_vc   <- rowSums(as.matrix(x0[,-1]*b_vc[,-1]))
    sf_pred <- b_vc[,1]
  }

  if( is.null(xconst0)){
    xb_const  <- 0
  } else {
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

  xconst_problem <- !is.null( mod$other$xf_id ) & is.null( xconst0 )
  if( ( !is.null( x0 ) ) & ( xconst_problem == FALSE ) ){
    	sf_pred	<- sf_pred - mean( mod$b_vc[,1] )
    	xb	    <- xb_vc + xb_const + mean( mod$b_vc[,1] )
    	pred	<- xb + sf_pred
    	if( is.null( b_g0 ) == FALSE ){
    		na_b_g0   <- is.na( b_g0 )
    		if( sum( na_b_g0 ) > 0 ){
    		  b_g0[ na_b_g0 ] <- 0
    		  message( "Note: NAs are given to the groups that are not in xgroup")
    		}
    		pred   <- pred + rowSums( b_g0 )
    	}

    	if( is.null(xgroup0) == FALSE ){
    		   res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred, "group" = b_g )
    	} else {
    		   res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred )
    	}
  } else if( is.null( x0 ) & ( xconst_problem == FALSE )){
    #sf_pred	<- sf_pred - mean( mod$b_vc[,1] )
    #xb	    <- xb_vc + xb_const + mean( mod$b_vc[,1] )
    #pred	<- xb + sf_pred
    #if( is.null( b_g0 ) == FALSE ){
    #  na_b_g0   <- is.na( b_g0 )
    #  if( sum( na_b_g0 ) > 0 ){
    #    b_g0[ na_b_g0 ] <- 0
    #    message( "Note: NAs are given to the groups that are not in xgroup")
    #  }
    #  pred   <- pred + rowSums( b_g0 )
    #}
    res	<- NULL
    message( "Note: y is not predicted because x0 and/or xconst0 is missing" )
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

  result	<- list( pred = res, b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc,
                  c_vc = c_vc, cse_vc = cse_vc, ct_vc = ct_vc, cp_vc = cp_vc )

 return( result )
}

