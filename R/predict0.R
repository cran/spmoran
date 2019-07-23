predict0<- function( mod, meig0, x0 = NULL, xgroup0 = NULL ){
    if( mod$other$model == "esf" ){
    	meig0$sf	<- as.matrix( meig0$sf[, mod$other$sf_id ] )
    }

    if( is.null( dim( mod$r ) ) ){
    	sf_pred	<- meig0$sf * mod$r[  1 ]
    } else {
    	sf_pred	<- meig0$sf %*% mod$r[, 1 ]
    }

    if( ( is.null( mod$b_g )==FALSE )&( mod$other$model != "esf" ) ){
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
    	xb_pred	<- c( as.matrix( mod$b[ 1 ] ) )
    	pred	  <- xb_pred + sf_pred
    } else {
    	if( is.null( x0 ) ){
    		message( " Note: Only spatial component (sf) is interpolated because x0 is missing")
    		pred	<- sf_pred
   	 } else {
   		 x02	  <- as.matrix( x0 )
  		 if( is.numeric( x02 ) == F ){
    			mode( x02 ) <- "numeric"
    	 }

   	 	xb_pred	<- as.matrix( cbind( 1, x02[ ,mod$other$x_id ] ) ) %*% mod$b[, 1 ]
    		pred	<- xb_pred + sf_pred
    		pred	<- as.data.frame( cbind( pred, xb_pred, sf_pred ) )
    		names( pred )<- c( "pred", "xb", "sf" )
    	}
    }

    if( is.null( b_g0 ) == FALSE ){
      pred      <- data.frame( pred, b_g0 )
      na_b_g0   <- is.na( b_g0 )
      if( sum( na_b_g0 ) > 0 ){
        b_g0[ na_b_g0 ] <- 0
        message( " Note: NAs are given to the groups that are not in xgroup")
      }
      pred[ ,1 ]<- pred[ ,1 ] + rowSums( b_g0 )
    }
    return( pred )
}

