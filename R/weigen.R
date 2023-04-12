
weigen <-function( x = NULL, type = "knn", k = 4,
                   threshold = 0.25, enum = NULL ){

    nb2mat2 <- function( nb, n = n ){

      listw     <- nb2listw( nb, style = "B" )$neighbours
      listw2    <- cbind( 0, unlist( listw ) )
      listw2_len<- sapply( listw, length )
      id_id	    <- 1
      for( i in 1:n ){
        id_end	<- id_id + listw2_len[ i ] - 1
        listw2[ id_id:id_end, 1 ] <- i
        id_id	<- id_end + 1
      }

      C   <- sparseMatrix(i={},j={},dims=c(n,n),x=0)
      C[ listw2 ]<-1
      return(C)
    }

    if( type %in% c( "knn", "tri" ) == FALSE ){
      stop( "type must be knn or tri" )
    }

    if( threshold > 1 | threshold < 0 ){
        stop( "threshold must be a value between 0 and 1" )
    } else if ( threshold ==1) {
        threshold <- threshold - 1e-07
    }

    if( "vector" %in% is( x ) ){
      x       <- as.matrix( x )
      coords  <- x
      polys   <- NULL
      cmat    <- NULL
      if( dim( x )[ 1 ] == dim( x )[ 2 ] ){
        coords<- NULL
        polys <- NULL
        cmat  <- x
      }
    } else if( "sf" %in% is( x ) ){#SpatialPolygons
      coords  <- NULL
      polys   <- x
      cmat    <- NULL
    }

    if( is.null( polys ) ){
      if( is.null( cmat ) ){
        n		  <- length( coords[ , 1 ] )
        if( type == "knn" ){
          message( "--------------knn-based W-------------" )
          knn	<- knearneigh( as.matrix( coords ), k = k )
          nb	<- knn2nb(knn)
          cmat<- nb2mat2( nb = nb, n = n )

        } else if( type == "tri" ){
          message( "---- Delaunay triangulation-based W ----" )
          nb	<- tri2nb( coords )
          cmat<- nb2mat2( nb = nb, n = n )
        }
      } else {
        n       <- dim( cmat )[ 1 ]
        if( isSymmetric( cmat ) == FALSE ){
          cmat  <- 0.5 * ( cmat + t( cmat ) )
        }
        diag( cmat ) <- 0
        message( "---------- User specified W -----------" )
      }

    } else {
          message( "----------- Adjacency-based W ----------" )
        n     <- dim(polys)[ 1 ]#length( polys )
        nb	  <- poly2nb( polys )
        cmat  <- nb2mat2( nb = nb, n = n )
        if( is.null( coords ) == FALSE ){
          message( " note: coords is ignored because polys is provided" )
        }
    }

    enum	        <- min( n - 1, ifelse( is.null( enum ), 200, enum ) )
    eigenC     	  <- eigs_sym( cmat, enum, which = "LA" )
    eigenC$values <- Re( eigenC$values )
    eigenC$vectors<- Re( eigenC$vectors )
    sel           <- ( eigenC$values / max( eigenC$values ) >= ( threshold + 1e-07) )
    if( sum( sel ) < enum ){
      eigenC$values	<- eigenC$values [ sel ]
      eigenC$vectors<- eigenC$vectors[ , sel ]
    } else {
      eigenC$values	<- eigenC$values [ 1:enum ]
      eigenC$vectors<- eigenC$vectors[ , 1:enum ]
    }

    other       <- list( wdum =1 )
    mes         <- paste( " ", length( eigenC$values ), "/", n, " eigen-pairs are extracted", sep = "" )
    message( mes )
    return( list( sf = eigenC$vectors, ev = eigenC$values, other = other ) )
}
