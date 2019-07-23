predict0_vc	<- function( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL ){
    if( mod$other$model == "resf" ){
    	res	<- predict0( mod = mod, meig0 = meig0, x0 = xconst0, xgroup0 = xgroup0 )
    	result	<- list( pred = res, b_vc = NULL, bse_vc = NULL, bt_vc = NULL, bp_vc = NULL, bg = NULL )

    } else if( mod$other$model == "resf_vc" ){
    	if( is.null( x0 ) == F ){
  		  if( is.numeric( x0 ) == F ){
    			x0		<- as.matrix( x0 )
    			mode( x0 )	<- "numeric"
    		}

    		if( length( mod$other$x_id ) == 1 ){
    			x0	<- as.matrix( cbind( 1, x0[  mod$other$x_id ] ))
    		} else {
    			x0	<- as.matrix( cbind( 1, x0[, mod$other$x_id ] ))
    		}

    		if( length( c( mod$vc ) ) != dim( as.matrix( x0 ) )[ 2 ] ){
    			stop( "x and x0 must have the same number of culumns" )
    		}

    		if( is.null( mod$other$xf_id ) == FALSE ){
    			if( is.null( xconst0 ) ){
    			  xb_const	<- 0
    			} else {
  			    if( is.numeric( xconst0 ) == F ){
    				  xconst0		<- as.matrix( xconst0 )
    				  mode( xconst0 )	<- "numeric"
    			  }

    			  if( length( mod$other$xf_id ) == 1 ){
    			  	xconst0	<- as.matrix( xconst0[  mod$other$xf_id ] )
    			  } else {
    			   	xconst0	<- as.matrix( xconst0[, mod$other$xf_id ] )
    			  }

    			  if( length( mod$other$xf_id ) != dim( as.matrix( xconst0 ) )[ 2 ] ){
    				  stop( "xconst and xconst0 must have the same number of culumns" )
    			  }

   			    if( is.null( dim( mod$b ) ) ){
    				  xb_const<- xconst0 * mod$b[  1 ]
    			  } else {
    				  xb_const<- xconst0 %*% mod$b[, 1 ]
    			  }
    		  }
    		} else {
    			xb_const	<- 0
    			if( is.null( xconst0 ) == FALSE ){
    				message( "Note: xconst0 is ignored" )
    			}
    		}
    	}

      if( is.null( mod$b_g )==FALSE ){
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


	    n	  <- length( mod$pred )
    	n0	<- length( meig0$sf[ , 1 ] )
    	ne	<- length( mod$other$b_s[[ 1 ]][ -1 ] )
    	nsv	<- sum( mod$other$x_id ) + 1
    	meig0$sf<- as.matrix( meig0$sf[ , 1:ne ])

    	xb_vc	<- 0
    	b_vc	<- matrix(0, nrow = n0, ncol = nsv )
    	bse_vc	<- matrix(0, nrow = n0, ncol = nsv )
    	bt_vc	<- matrix(0, nrow = n0, ncol = nsv )
    	bp_vc	<- matrix(0, nrow = n0, ncol = nsv )
    	for( i in 1:nsv ){
    		if( length( mod$other$evSqrts[[ i ]] ) != ne ){
        		b_vc[ , i ]	<- mod$other$b_s[[ i ]][ 1 ]
   			    bse_vc[ , i ]	<- sqrt( mod$other$b_covs[[ i ]] )
        		bt_vc[ , i ]	<- b_vc[ , i ] / bse_vc[ , i ]
        		bp_vc[ , i ]	<- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    		} else {
        		b_vc[ , i ]	<- mod$other$b_s[[ i ]][ 1 ] + meig0$sf %*% c( mod$other$b_s[[ i ]][ -1 ] )
        		sf2		<- t( t( meig0$sf ) * mod$other$evSqrts[[ i ]] )
        		x_sf		<- as.matrix( cbind( 1, sf2 ) )
   			    bse_vc[ , i ]	<- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs[[ i ]] %*% t( x_sf ) ) ) )
        		bt_vc[ , i ]	<- b_vc[ , i ] / bse_vc[ , i ]
        		bp_vc[ , i ]	<- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    		}

    		if( is.null( x0 ) == F ){
    			if( i == 1 ){
    				sf_pred	<- b_vc[, i ]
    			} else {
   				xb_vc	<- xb_vc + x0[ , i ]* b_vc[ , i ]
    			}
    		}
    	}

    	xconst_problem <- ( is.null( mod$other$xf_id )==FALSE ) & is.null( xconst0 )
    	if( ( is.null( x0 ) == FALSE ) & ( xconst_problem == FALSE ) ){
    		sf_mean	<- mean( sf_pred )
    		sf_pred	<- sf_pred - sf_mean
    		xb	<- xb_const + xb_vc + sf_mean
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
    		   res	<- data.frame( "pred" = pred, "xb" = xb, "sf" = sf_pred, b_g )
    		} else {
    		   res	<- data.frame( "pred" = pred, "xb" = xb, "sf" = sf_pred )
    		}
    	} else {
    		res	<- NULL
    		message( "Note: y is not predicted because x0 and/or xconst0 is missing" )
    	}

    	b_vc	<- as.data.frame( b_vc )
    	bse_vc	<- as.data.frame( bse_vc )
    	bt_vc	<- as.data.frame( bt_vc )
    	bp_vc	<- as.data.frame( bp_vc )
    	names( b_vc )	<- names( mod$b_vc )
    	names( bse_vc )<- names( mod$b_vc )
    	names( bt_vc )	<- names( mod$b_vc )
    	names( bp_vc )	<- names( mod$b_vc )
    	result	<- list( pred = res, b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc )
    }
    return( result )
}

