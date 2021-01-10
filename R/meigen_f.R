meigen_f	<- function( coords, model = "exp", enum = 200, s_id = NULL ){
    n		<- nn<- dim( coords )[ 1 ]
    if( enum < 1 | enum != floor( enum ) ){
    	stop( " enum must be a positive integer" )
    }

    if( !is.null( s_id )[ 1 ] ){
      coords0 <-coords
      coords_x<-tapply(coords[,1],s_id,mean)
      coords_y<-tapply(coords[,2],s_id,mean)
      coords  <-as.matrix(cbind(coords_x, coords_y))
      s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
      s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
      s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
      s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
      nn      <- length( coords_x )
    }

    if( enum >= nn - 1 ){
      if( !is.null( s_id )[ 1 ] ) coords  <- coords0
    	result	<- meigen(coords = coords, model = model, s_id = s_id )
    	result$other$coordk	<- NULL
    	result$other$sfk	<- NULL
    	result$other$evk	<- NULL

    } else {
    	coordk	<- kmeans( coords, centers = enum + 1 )$centers
    	D	<- rdist( coordk )
    	h	<- max( spantree( D )$dist )
    	if( model == "exp" ){
    		C	<- exp( -D / h )
    	} else if( model == "gau" ){
    		C	<- exp( -( D / h) ^ 2 )
    	} else if( model == "sph" ){
    		C	<- ifelse( D < h , 1 - 1.5 * (D / h ) + 0.5 * ( D / h ) ^ 3, 0 )
    	} else {
    		stop( "model is not specified appropriately" )
    	}

    	Cmean	<- apply( C, 1, mean )
    	MCM	<- t( C - Cmean ) - Cmean + mean( Cmean )
    	eigenC	<- eigen( MCM )
    	sf0	<- eigenC$vectors[ , eigenC$values > 1e-08 ]
    	ev0	<- eigenC$values [   eigenC$values > 1e-08 ]
    	ev_full	<- ev0 * ( n / enum ) - 1
    	ev_ap	<- ev_full[   ev_full > 1e-08 ]

    	sf_ap0	<- sf0 %*% diag( 1 / ev0 )
    	#sf_ap0	<- t( exp( -rdist( coordk, coords ) / h ) - Cmean ) %*% sf_ap0

    	if( model == "exp" ){
    	  ccc0	<- exp( -rdist( coordk, coords ) / h )
    	} else if( model == "gau" ){
    	  C	    <- exp( -( D / h) ^ 2 )
    	  ccc0	<- exp( -( rdist( coordk, coords ) / h )^2 )
    	} else if( model == "sph" ){
    	  ddd0  <- rdist( coordk, coords )
    	  ccc0  <- ifelse( ddd0 < h , 1 - 1.5 * ( ddd0 / h ) + 0.5 * ( ddd0 / h ) ^ 3, 0 )
    	}
    	sf_ap0  <- t( ccc0 - Cmean ) %*% sf_ap0
    	sf_ap0	<- t( t( sf_ap0 ) - colMeans( sf_ap0 ) )
    	sf_ap	<- as.matrix( sf_ap0[ , ev_full > 1e-08 ] )

    	if( !is.null( s_id )[ 1 ]){
    	  sf_ap		<- sf_ap[ s_id_dat2$s_id_num,]
    	  #Cmean   <- Cmean[ s_id_dat2$s_id_num ]
    	  coords  <- coords0
    	}

   	  sf0	<- as.matrix( sf0   [ , ev_full > 1e-08 ] )
    	ev0	<- ev0   [   ev_full > 1e-08 ] -1
    	other	<- list( coordk = coordk, sfk = sf0, evk = ev0, Cmean = Cmean, h = h, model = model, fast = 1,
    	               coords = coords, s_id = s_id)

    	mes	<- paste( " ", length( ev_ap ), "/", nn, " eigen-pairs are approximated", sep = "" )
    	message( mes )
    	result	<- list( sf  =sf_ap, ev = ev_ap, ev_full = ev_full, other = other )
    }
    return( result )
}
