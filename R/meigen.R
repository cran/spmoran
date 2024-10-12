meigen	<- function( coords = NULL,model = "exp", enum = NULL, s_id = NULL,
                     threshold = 0, cmat = NULL, coords_z=NULL, enum_z=NULL,
                     interact=TRUE, interact_max_dim = 600 ){#

  if( threshold > 1 | threshold < 0 ){
    stop( "threshold must be between 0 and 1" )
  } else if ( threshold ==1) {
    threshold <- threshold - 1e-07
  }

  #if( is.null( coords_z ) ) interact <- FALSE

  if( is.null( cmat ) ){
    if( !is.null( s_id )[ 1 ] ){
      coords0 <-coords
      coords_x<-tapply(coords[,1],s_id,mean)
      coords_y<-tapply(coords[,2],s_id,mean)
      coords  <-as.matrix(cbind(coords_x, coords_y))
      s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
      s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
      s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
      s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    } else {
      coords_uni<- unique(coords)
      n         <- dim(coords)[1]
      nn        <- dim(coords_uni)[1]
      if( n > nn ){
        coords0 <- coords
        s_id    <- get.knnx(coords_uni,coords,1)$nn.index
        coords_x<-tapply(coords[,1],s_id,min)# min just because it seems the fastest

        coords  <- coords_uni
        s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
        s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
        s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
        s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
      }
    }

    D	<- rdist( coords )
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
  } else {
    coords0 <-NULL
    if( !is.null( s_id )[ 1 ] ){
      s_id_uni<-unique(s_id)
      s_id2   <-data.frame( s_id=s_id_uni, s_id_num = 1:length(s_id_uni) )
      s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
      s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
      s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    }

    if( isSymmetric( unname( cmat ) )==FALSE ) {
      C 	<- 0.5 * ( cmat + t( cmat ) )
      message( " Note: cmat is symmetrized by ( cmat + t( cmat ) ) / 2" )
    } else {
      C	  <- as.matrix( cmat )
    }
    model	<- "other"
    coords<- NULL
    h	    <- NULL
  }

  diag( C )<- 0
  Cmean	   <- apply( C, 1, mean )
  MCM		   <- t( C - Cmean ) - Cmean + mean( Cmean )
  eigenC	 <- eigen( MCM )
  eigenC$values	<- Re( eigenC$values )
  eigenC$vectors<- Re( eigenC$vectors )
  sel		   <- ( eigenC$values >= threshold + 1e-07 )# / max( eigenC$values )
  if( is.null( enum ) == FALSE ){
    if( sum(sel) > enum ){
      sel[ -c(1:enum) ] <- FALSE
    }
  }
  sf		   <- as.matrix( eigenC$vectors[ , sel ] )
  ev		   <- eigenC$values [ sel ]# / max( eigenC$values [ sel ] )
  if( !is.null( s_id )[ 1 ]){
    sfk    <- sf
    Cmeank <- Cmean#/ max( eigenC$values [ sel ] ) #######
    coordk <- coords

    sf		 <- sf[ s_id_dat2$s_id_num,]
    Cmean  <- Cmean[ s_id_dat2$s_id_num ]#/ max( eigenC$values [ sel ] ) #######
    coords <- coords0
  } else {
    sfk    <-NULL
    Cmeank <-NULL
    coordk <-NULL
  }
  mes		<- paste0( " ", length( ev ), " spatial eigen-pairs" )
  message( mes )

  z_use  <- NULL
  if( !is.null(coords_z) ){
    coords_z  <- as.matrix(coords_z)
    z_uni     <- apply(coords_z,2,function(x) length(unique(x)))
    coords_z  <- as.matrix(coords_z[,z_uni>=4])
    if(min(z_uni) < 4 ){
      z_om<-paste(which( z_uni < 4 ), collapse = ",")
      message( paste0("   Note: coords_z[,c(",z_om,")] is ignored because its unique value is too small (<=3)" ))
      if(dim(coords_z)[2]==0) coords_z<-NULL
      z_use  <- z_uni >= 4
      #if(sum(z_use)==0) z_use<-NULL
    }
  }

  if( !is.null(coords_z) ){
    ev_z <- sf_z <- ev0_z <- sf0_z <- h_z <- Cmean_z <- id_z <- coordk_z<-list( NULL )
    nz        <- dim(coords_z)[2]

    for(j in 1:nz){
      coords_z2   <- sort(unique(coords_z[,j]))
      nn_z        <- length( coords_z2 )
      coords_z2_id<- 1:nn_z

      if( nn_z < 2000 ){
        Dz       <- rdist( coords_z2 )
        hz	     <- max( spantree( Dz )$dist )
        if( model == "exp"|!is.null(cmat) ){
          Cz	<- exp( -Dz / hz )
        } else if( model == "gau" ){
          Cz	<- exp( -( Dz / hz) ^ 2 )
        } else if( model == "sph" ){
          Cz	<- ifelse( Dz < hz , 1 - 1.5 * ( Dz / hz ) + 0.5 * ( Dz / hz ) ^ 3, 0 )
        } else {
          stop( "model is not specified appropriately" )
        }
        diag( Cz )<- 0
        Cmeanz	  <- apply( Cz, 1, mean )
        MCMz		  <- t( Cz - Cmeanz ) - Cmeanz + mean( Cmeanz )
        eigenCz	  <- eigen( MCMz )
        selz		  <- ( eigenCz$values ) >= threshold + 1e-07# / max( eigenCt$values )
        selz_max  <- min( round((nn_z + 0.1)/2), 6, ncol(Cz)-1 )
        if( sum(selz) <  selz_max ){
          selz_sel<-abs(eigenCz$values) > threshold + 1e-07
          eigenCz$values[selz_sel]<-(eigenCz$values[selz_sel]+1)/(eigenCz$values[selz_sel][selz_max+1]+1)-1
          selz		  <- ( eigenCz$values ) >= threshold + 1e-07# / max( eigenCt$values )
        }

        if( !is.null( enum_z ) ){
          if( sum(selz) > enum_z ) selz[ which(selz)[-(1:enum_z)] ] <- FALSE
        }

        nn_z        <- get.knnx(coords_z2, coords_z[,j], 1)$nn.index
        sf_z[[j]]   <- as.matrix( Re( eigenCz$vectors[nn_z,selz]) )
        ev_z[[j]]   <- Re( eigenCz$values[selz] ) #################/max(eigenCz$values[selz])
        id_z[[j]]   <- coords_z2_id[nn_z]
        coordk_z[[j]]<- coords_z2#coords_z[,j]
        h_z[[j]]    <- hz
        Cmean_z[[j]]<-Cmeanz#/max(eigenCz$values[selz]) ###################

      } else {

        suppressWarnings(coordk_z2	<- kmeans( coords_z2, centers = enum_z + 1 )$centers)
        Dz	     <- rdist( coordk_z2 )
        hz	     <- max( spantree( Dz )$dist )
        if( model == "exp" ){
          Cz	<- exp( -Dz / hz )
        } else if( model == "gau" ){
          Cz	<- exp( -( Dz / hz) ^ 2 )
        } else if( model == "sph" ){
          Cz	<- ifelse( Dz < hz , 1 - 1.5 * ( Dz / hz ) + 0.5 * ( Dz / hz ) ^ 3, 0 )
        } else {
          stop( "model is not specified appropriately" )
        }
        #diag( Cz )<- 0: this is for approximation
        Cmeanz	  <- apply( Cz, 1, mean )
        MCMz		  <- t( Cz - Cmeanz ) - Cmeanz + mean( Cmeanz )
        eigenCz	  <- eigen( MCMz )
        sf0	      <- eigenCz$vectors[ , eigenCz$values > 1e-08 ]
        ev0	      <- eigenCz$values [   eigenCz$values > 1e-08 ]
        ev_full	  <- ev0 * ( nn_z / enum_z ) - 1
        selz		  <- ev_full >= 1e-08
        ev_ap_z	  <- ev_full[   selz ]
        sf_ap0	  <- sf0 %*% diag( 1 / ev0 )

        if( model == "exp" ){
          ccc0	<- exp( -rdist( coordk_z2, coords_z2 ) / hz )
        } else if( model == "gau" ){
          C	    <- exp( -( D / hz) ^ 2 )
          ccc0	<- exp( -( rdist( coordk_z2, coords_z2 ) / hz )^2 )
        } else if( model == "sph" ){
          ddd0  <- rdist( coordk_z2, coords_z2 )
          ccc0  <- ifelse( ddd0 < hz , 1 - 1.5 * ( ddd0 / hz ) + 0.5 * ( ddd0 / hz ) ^ 3, 0 )
        }
        sf_ap0  <- t( ccc0 - Cmean ) %*% sf_ap0
        sf_ap0	<- t( t( sf_ap0 ) - colMeans( sf_ap0 ) )
        sf_ap_z	  <- as.matrix( sf_ap0[ , selz ] )

        #sfsf        <- rowSums( sf_ap_z^2 )
        #sfsf[sfsf>1]<-1

        nn_z        <- get.knnx(coords_z2, coords_z[,j], 1)$nn.index
        sf_z[[j]]   <- sf_ap_z[nn_z,]
        ev_z[[j]]   <- ev_ap_z#/max(ev_ap_z) ###################
        sf0_z[[j]]  <- as.matrix( sf0[ , selz ] )
        ev0_z[[j]]  <- ev0[  selz ] -1#/max(ev_ap_z) ###################
        id_z[[j]]   <- coords_z2_id[nn_z]
        coordk_z[[j]]<- coordk_z2
        h_z[[j]]    <- hz
        Cmean_z[[j]]<-Cmeanz#/max(ev_ap_z)###################
      }
      mes		<- paste0( " ", length( ev_z[[j]] ), " eigen-pairs for coords_z[,", j,"]" )
      message( mes )
    }
  } else {
    ev_z <- sf_z <- ev0_z <- sf0_z <- h_z <- Cmean_z <- id_z <-coordk_z<- NULL
    interact<-FALSE
  }

  other	<- list( coords = coords, Cmean = Cmean,
                 h = h, model = model, fast = 0,s_id = s_id,
                 sfk = sfk, evk = ev, Cmeank = Cmeank, coordk = coordk, coordk_z=coordk_z,
                 h_z = h_z, coords_z=coords_z, Cmean_z=Cmean_z, sfk_z=sf0_z, evk_z=ev0_z, z_use=z_use,
                 interact=interact, interact_max_dim=interact_max_dim, id_z=id_z )
  #other$coordk	<- NULL
  #other$sfk	<- NULL
  #other$evk	<- NULL

  result<-list( sf = sf, ev = ev, sf_z=sf_z, ev_z=ev_z, other = other )
  return( result )
}

