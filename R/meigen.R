meigen		<- function( coords = NULL, model = "exp", threshold = 0, enum = NULL, cmat = NULL, s_id = NULL ){
  if( threshold > 1 | threshold < 0 ){
    stop( "threshold must be a value between 0 and 1" )
  } else if ( threshold ==1) {
    threshold <- threshold - 1e-07
  }

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
  sel		   <- ( eigenC$values / max( eigenC$values ) >= threshold + 1e-07 )
  if( is.null( enum ) == FALSE ){
    if( sum(sel) > enum ){
      sel[ -c(1:enum) ] <- FALSE
    }
  }
  sf		   <- as.matrix( eigenC$vectors[ , sel ] )
  if( !is.null( s_id )[ 1 ]){
    sfk    <- sf
    Cmeank <- Cmean
    coordk <- coords

    sf		 <- sf[ s_id_dat2$s_id_num,]
    Cmean  <- Cmean[ s_id_dat2$s_id_num ]
    coords <- coords0
  } else {
    sfk    <-NULL
    Cmeank <-NULL
    coordk<-NULL
  }

  ev		<- eigenC$values [ sel ]
  other	<- list( coords = coords, Cmean = Cmean, h = h, model = model, fast = 0, s_id = s_id,
                 sfk = sfk, Cmeank = Cmeank, coordk = coordk )
  mes		<- paste( " ", length( ev ), "/", length( eigenC$values ), " eigen-pairs are extracted", sep = "" )
  message( mes )
  return( list( sf = sf, ev = ev, ev_full = eigenC$values, other = other ) )
}

