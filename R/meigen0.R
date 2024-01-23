
meigen0	<- function( meig, coords0, s_id0 = NULL ){

  if( meig$other$model == "other" ){
    	stop( "meigen0 is not supported for user-specified spatial proximity matrix" )
  }

  if( !is.null( meig$other$s_id )[ 1 ] ){
    if( is.null( s_id0 ) ) {
      stop( "s_id0 must be provided" )
    } else {
      coords0_0 <-coords0
      coords0_x<-tapply(coords0[,1],s_id0,mean)
      coords0_y<-tapply(coords0[,2],s_id0,mean)
      coords0  <-as.matrix(cbind(coords0_x, coords0_y))
      s_id2   <-data.frame( s_id=rownames(coords0_x), s_id_num = 1:length(coords0_x) )
      s_id_dat<-data.frame( s_id = s_id0, id = 1:length(s_id0))
      s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
      s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    }
  }

  if( meig$other$fast == 0 ){
    if( !is.null( meig$other$s_id )[ 1 ] ){
      sfk   <-meig$other$sfk
      evk   <-meig$ev
      coordk<-meig$other$coordk
      Cmean <-meig$other$Cmeank
    } else {
      sfk	  <- meig$sf
      evk	  <- meig$ev
      coordk<- meig$other$coords
      Cmean <- meig$other$Cmean
    }
  } else {
    	sfk	<- meig$other$sfk
    	evk	<- meig$other$evk
    	coordk	<- meig$other$coordk
    	Cmean <- meig$other$Cmean
  }
  h		<- meig$other$h
  nm	<- dim( coords0 )[ 1 ]
  sfk	<- sfk %*% diag( 1 / (evk+1) )
  Dk	<- rdist( coordk, coords0 )
  if( meig$other$model == "exp" ){
    	sfk	<- t( exp( -Dk / h ) - Cmean ) %*% sfk
  } else if( meig$other$model == "gau" ){
    	sfk	<- t( exp( - ( Dk / h ) ^ 2 ) - Cmean ) %*% sfk
  } else if( meig$other$model == "sph" ){
    	sfk	<- t( ifelse( Dk < h , 1 - 1.5 * ( Dk / h ) + 0.5 * ( Dk / h ) ^ 3, 0 ) - Cmean ) %*% sfk
  }
  if( !is.null( meig$other$s_id )[ 1 ] ){
      sfk		<- sfk[ s_id_dat2$s_id_num,]
  }

  other = list( coords0 = coords0 )
  return( list( sf = sfk, ev = meig$ev, ev_full = meig$ev_full, other = other )  )
}
