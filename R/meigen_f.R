meigen_f<- function( coords, model = "exp", enum = 200, s_id = NULL, threshold = 0,
                     coords_z = NULL, enum_z = 200, interact = TRUE,
                     interact_max_dim = 600){

  n		<- nn<- dim( coords )[ 1 ]
  if( enum < 1 | enum != floor( enum ) ){
    stop( " enum must be a positive integer" )
  }

  #if( is.null( coords_z ) ) interact <- FALSE

  coords0   <-coords
  if( !is.null( s_id )[ 1 ] ){
    coords_x<-tapply(coords0[,1],s_id,mean)
    coords_y<-tapply(coords0[,2],s_id,mean)
    coords  <-as.matrix(cbind(coords_x, coords_y))
    s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
    s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
    s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
    s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    nn      <- length( coords_x )
  } else {
    coords_uni<- unique(coords0)
    nn        <- dim(coords_uni)[1]
    if( n > nn ){
      s_id    <- get.knnx(coords_uni,coords0,1)$nn.index
      coords_x<-tapply(coords0[,1],s_id,min)# min just because it seems the fastest

      coords  <- coords_uni
      s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
      s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
      s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
      s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    }
  }

  if( enum >= nn - 1| nn < 2000  ){
    result	<- meigen(coords = coords0, model = model, s_id = s_id,
                      coords_z = coords_z, enum_z = enum_z, interact = interact )
    #result$other$coordk	<- NULL
    #result$other$sfk	<- NULL
    #result$other$evk	<- NULL

  } else {

    suppressWarnings(coordk	<- kmeans( coords, centers = enum + 1 )$centers)

    D	<- rdist( coordk )
    h	<- max( spantree( D )$dist )#if( is.null(h) ) h	<- max( spantree( D )$dist )
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
    ev_full	<- ev0 * ( nn / enum ) - 1
    ev_ap	  <- ev_full[   ev_full > 1e-08 ]
    ev_ap   <- ev_ap#/max(ev_full) ###############################

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
    sf_ap	  <- as.matrix( sf_ap0[ , ev_full > 1e-08 ] )

    if( !is.null( s_id )[ 1 ]){
      sf_ap		  <- sf_ap[ s_id_dat2$s_id_num,]
      #Cmean    <- Cmean[ s_id_dat2$s_id_num ]
    }
    coords      <- coords0
    sfsf        <- rowSums( sf_ap^2 )
    sfsf[sfsf>1]<-1

    sf0	        <- as.matrix( sf0   [ , ev_full > 1e-08 ] )
    ev0	        <- ev0   [   ev_full > 1e-08 ] -1
    ev0         <- ev0#/max(ev_full) ###############################
    Cmean       <- Cmean#/max( eigenCz$values[selz] ) ###############################
    #other	<- list( coordk = coordk, sfk = sf0, evk = ev0, Cmean = Cmean, h = h, model = model, fast = 1,
    #               coords = coords, s_id = s_id, sfsf=sfsf)

    mes		 <- paste0( " ", length( ev_ap ), " spatial eigen-pairs" )
    message( mes )
    #result	<- list( sf  =sf_ap, ev = ev_ap, ev_full = ev_full )

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
      }
    }

    if( !is.null(coords_z) ){
      ev_z      <- sf_z  <-ev0_z <- sf0_z    <- h_z <-Cmean_z<- id_z<-coordk_z<-list( NULL )
      nz        <- dim(coords_z)[2]
      for(j in 1:nz){
        coords_z2  <- sort(unique(coords_z[,j]))
        nn_z       <- length( coords_z2 )
        coords_z2_id<- 1:nn_z

        if( nn_z < 2000 ){
          Dz       <- rdist( coords_z2 )
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
          diag( Cz )<- 0
          Cmeanz	  <- apply( Cz, 1, mean )
          MCMz		  <- t( Cz - Cmeanz ) - Cmeanz + mean( Cmeanz )
          eigenCz	  <- eigen( MCMz )
          selz		  <- ( eigenCz$values ) >= threshold + 1e-07
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
          sf_z[[j]]   <- as.matrix( Re( eigenCz$vectors[nn_z,selz] ))
          ev_z[[j]]   <- Re( eigenCz$values[selz])#/max(eigenCz$values[selz] )) ###############################
          id_z[[j]]   <- coords_z2_id[nn_z]
          h_z[[j]]    <- hz
          coordk_z[[j]]<-coords_z2
          Cmean_z[[j]]<-Cmeanz#/max( eigenCz$values[selz] ) ###############################

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
          #diag( Cz )<- 0
          Cmeanz	  <- apply( Cz, 1, mean )
          MCMz		  <- t( Cz - Cmeanz ) - Cmeanz + mean( Cmeanz )
          eigenCz	  <- eigen( MCMz )
          sf0_z0	  <- eigenCz$vectors[ , eigenCz$values > 1e-08 ]
          ev0	      <- eigenCz$values [   eigenCz$values > 1e-08 ]
          ev_full	  <- ev0 * ( nn_z / enum_z ) - 1
          selz		  <- ev_full >= 1e-08
          selz_max  <- min( round(nn_z/2), 6 )
          if( sum(selz) <  selz_max ){
            ev_full <- (eigenCz$values+1)/(eigenCz$values[selz_max+1]+1)-1
            selz		<- ev_full >= 1e-08# / max( eigenCt$values )
          }
          ev_ap_z	  <- ev_full[   selz ]
          sf_ap0	  <- sf0_z0 %*% diag( 1 / ev0 )

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
          ev_z[[j]]   <- ev_ap_z#/max(ev_ap_z) ###############################
          sf0_z[[j]]  <- as.matrix( sf0_z0[ , selz ] )
          ev0_z[[j]]  <- ev0[  selz ] -1#/max(ev_ap_z) ###############################
          id_z[[j]]   <- coords_z2_id[nn_z]
          h_z[[j]]    <- hz
          coordk_z[[j]]<-coordk_z2
          Cmean_z[[j]] <-Cmeanz /max(ev_ap_z) ##################3
        }
        mes		<- paste0( " ", length( ev_z[[j]] ), " eigen-pairs for coords_z[,", j,"]" )
        message( mes )
      }
    } else {
      ev_z<- sf_z  <-ev0_z<-sf0_z    <- h_z <-Cmean_z<- id_z<-coordk_z<-NULL
      interact     <-FALSE
    }

    other	<- list( coords = coords, Cmean = Cmean,sfsf=sfsf,
                   h = h, model = model, fast = 1,s_id = s_id,
                   sfk = sf0, evk=ev0, coordk = coordk,#Cmeank = Cmeank,
                   h_z = h_z, coords_z=coords_z, coordk_z=coordk_z,
                   Cmean_z=Cmean_z, sfk_z=sf0_z, evk_z=ev0_z,z_use=z_use,
                   interact=interact, interact_max_dim=interact_max_dim, id_z =id_z )

    result<-list( sf = sf_ap, ev = ev_ap, sf_z=sf_z, ev_z=ev_z, other = other )
  }

  return( result )
}
