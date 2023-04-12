
besf_vc		<- function( y, x, xconst = NULL, coords, s_id = NULL,
                      x_nvc = FALSE, xconst_nvc = FALSE,
                      x_sel = TRUE, x_nvc_sel = TRUE, xconst_nvc_sel = TRUE, nvc_num=5,
                      method = "reml", penalty = "bic", maxiter = 30,
                      covmodel="exp",enum = 200, bsize = 4000, ncores=NULL ){

  n     	    <- length( y )
  nvc_x       <- x_nvc
  nvc_xconst  <- xconst_nvc
  xgroup      <- NULL
  ng	        <- 0

  if( is.logical( x_sel ) ){
    allsvc     <- !x_sel
    allsvc_flag<-0
  } else {
    allsvc    <- x_sel
    allsvc_flag<-1
  }


  if( is.logical( x_nvc_sel ) ){
    allnvc_x    <- !x_nvc_sel
    allnvc_x_flag<-0
  } else {
    allnvc_x    <- x_nvc_sel
    allnvc_x_flag<-1
  }


  if( is.logical( xconst_nvc_sel ) ){
    allnvc_xconst<- !xconst_nvc_sel
    allnvc_xconst_flag<-0
  } else {
    allnvc_xconst<- xconst_nvc_sel
    allnvc_xconst_flag<-1
  }

  #if( is.logical( sel_xgroup ) ){
  #  allg<- !sel_xgroup
  #} else {
  #  allg<- sel_xgroup
  #}

  if( x_nvc[1]      == FALSE ) x_nvc_sel = FALSE
  if( xconst_nvc[1] == FALSE ) xconst_nvc_sel = FALSE

  sel_basis_vif <-function(test0, vif_max=15){
    tres3    <-TRUE
    rem_ord  <-NULL
    while(tres3){
      suppressWarnings(cor_test0<-cor(test0))
      testt2<-try(vif <- diag(solve(cor_test0)), silent=TRUE)
      #testt2<-try(vif <- diag(solve(cor(test0))), silent=TRUE)
      if( inherits(testt2, "try-error") ){#class(testt2)=="try-error"
        suppressWarnings(corr_0<-cor(test0))
        #corr_0   <- cor(test0)
        diag(corr_0)<- 0
        cor_max  <- apply(corr_0,2,max)
        rem_ord0 <- which(cor_max == max(cor_max))[1]
        test0    <- test0[, -rem_ord0]
        rem_ord  <- c( rem_ord, rem_ord0 )

      } else {
        if( max(vif)>vif_max ){
          rem_ord0<- which( vif == max( vif ))[ 1 ]
          test0   <- test0[ ,-rem_ord0 ]
        } else {
          tres3<-FALSE
        }
      }
    }
    return(test0)#list(selvar=test0, rem_ord=rem_ord)
  }

  Mdet_f0	  	<- function( M, M0, id, par0_sel, emet ){
    if( emet == "ml" ){
      M	  <- M [ id != 0, id != 0 ]
      M0	<- M0[ id != 0, id != 0 ]
      id	<- id[ id != 0 ]
    }
    if(sum( id != par0_sel ) == 0 ){
      term2	  <- NULL
      term3_0	<- M[ id == par0_sel, id == par0_sel ]
    } else {
      M0_sub	<- as.matrix( M0[ id != par0_sel, id != par0_sel ] )
      term2		<- determinant( M0_sub )$modulus
      Msub_00	<- M[ id == par0_sel, id == par0_sel ]
      Msub_01	<- M[ id == par0_sel, id != par0_sel ]
      term3_0	<- Msub_00 - Msub_01 %*% solve( M0_sub, tol = 1e-30 ) %*% t( Msub_01 )#, tol = 1e-30
    }
    return(list(term2 = term2, term3_0 = term3_0))
  }

  Mdet_f	  	<- function( evSqrt, id, term2, term3_0, par0_sel ){
    term1		<- sum( log( evSqrt ) ) * 2
    diag( term3_0 ) <- diag( term3_0 ) + 1/evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
    term3	  <- determinant( term3_0 )$modulus
    if( is.null( term2 )  ){
      Mdet    <- term1 + term3
    } else {
      Mdet	  <- term1 + term2 + term3
    }
    return(Mdet)
  }

  lik_resf_vc		<- function( par0, par0_est, par0_id, par0_sel, ev, M, M0inv, M0inv_01, M0inv_00,
                            m, yy, b_01, b_02, n, nx, nsv, nnxf, nnsv, ng, emet, term2, term3_0, null_dum2, id ){
    par		<- par0 ^ 2
    par_est		<- par0_est ^ 2
    par[ par0_id == par0_sel ]  <- par_est
    evSqrt	<- NULL
    for( i in ( 1:nsv )[ null_dum2[ 1:nsv ] == 0 ] ){
      evv  	<- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
      evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
    }

    if( ng != 0 ){
      for( j in (1:ng)[ null_dum2[ (( nsv + 1 ):( nsv + ng )) ] == 0 ] ){#null_dum2[ -(1:nsv) ] == 0
        xgg	<- rep( 1, sum( id == nsv + j ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + j ] * xgg )
      }
    }

    if( nnxf != 0 ){
      for( i2 in (1:nnxf)[ null_dum2[ (nsv+ng+1):(nsv+ng+nnxf) ] == 0 ] ){
        xgg	<- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + ng + i2 ] * xgg )
      }
    }

    if( nnsv != 0 ){
      for( i2 in (1:nnsv)[ null_dum2[ (nsv+ng+nnxf+1):(nsv+ng+nnxf+nnsv) ] == 0 ] ){
        xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + ng + nnxf + i2 ] * xgg )
      }
    }

    Mdet		<- Mdet_f( id = id, par0_sel=par0_sel, term2 = term2, term3_0 = term3_0, evSqrt = evSqrt )
    M2		<- M
    for( j in 1:max( id ) ){
      diag( M[ id == j, id == j ])<-diag( M[ id == j, id == j ] ) + 1/evSqrt[ id[ -c( 1:nx ) ] == j ] ^ 2
    }

    diag(M0inv_00)	<- diag( M0inv_00 ) + evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
    b_02_b		<- solve( M0inv_00, tol = 1e-30 ) %*% b_02
    b_02		<- M0inv_01 %*% b_02_b
    b		<- b_01 - b_02
    sse		<- yy - 2 * t( b ) %*% m + t( b ) %*% M2 %*% b
    dd		<- sse + sum( ( b[ -( 1:nx ) ] / evSqrt ) ^ 2 )
    if( emet == "reml" ){
      loglik	<- Mdet + ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
    } else if( emet == "ml" ){
      loglik	<- Mdet + n * ( 1 + log( 2 * pi * dd / n ) )
    }
    return( loglik )
  }

  lik_resf_vc0	<- function( par0, ev, M, m, yy, n, nx, nsv, ng, nnxf, nnsv, emet, null_dum4 ){
    par		<- par0 ^ 2
    evSqrt		<- NULL
    for( i in ( 1:nsv )[ null_dum4[ 1:nsv ] == 0 ] ){
      evv  	<- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
      evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
    }

    if( ng != 0 ){
      for( j in (1:ng)[ null_dum4[ (( nsv + 1 ):( nsv + ng )) ] == 0 ] ){#null_dum4[ -(1:nsv) ] == 0
        xgg	<- rep( 1, sum( id == nsv + j ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + j ] * xgg )
      }
    }
    if( nnxf != 0 ){
      for( i2 in (1:nnxf)[ null_dum4[ (nsv+ng+1):(nsv+ng+nnxf) ] == 0 ] ){
        xgg	<- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + ng + i2 ] * xgg )
      }
    }

    if( nnsv != 0 ){
      for( i2 in (1:nnsv)[ null_dum4[ (nsv+ng+nnxf+1):(nsv+ng+nnxf+nnsv) ] == 0 ] ){
        xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
        evSqrt	<- c( evSqrt, par[ 2 * nsv + ng + nnxf + i2 ] * xgg )
      }
    }

    M[ -( 1:nx ), -( 1:nx ) ]	<- t(M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    M[ -( 1:nx ),    1:nx   ]	<-   M[ -( 1:nx ),    1:nx   ] * evSqrt
    M[    1:nx  , -( 1:nx ) ]	<- t(M[ -( 1:nx ),    1:nx   ] )
    M0		<- M
    diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1

    m[-(1:nx)]		<- m[ -( 1:nx ) ] * evSqrt
    test			<-try(Minv	<- solve( M, tol = 1e-30 ))
    #if(class(test)[1]=="try-error"){
    if( inherits(test,"try-error") ){
      loglik  	<- Inf
    } else {
      b		<- Minv %*% m
      sse		<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
      dd		<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
      if( emet == "reml" ){
        term1	<- determinant( M )$modulus
        term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
      } else if( emet == "ml" ){
        term1	<- determinant( as.matrix( M[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
        term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
      }
      loglik		<- term1 + term2
    }
    return( loglik[ 1 ] )
  }

  if( method == "reml" ){
    lik_nam	<- "rlogLik"
  } else if( method == "ml" ){
    lik_nam	<- "logLik"
  }


  if(is.null(xconst)){
    nvc_xconst <- FALSE
  }

  nx0   <- 0
  X1	  <- x
  if( is.null( X1 ) == FALSE ){
    X1	<- as.matrix( X1 )
    if( is.numeric( X1 ) == FALSE ){
      mode( X1 ) <- "numeric"
    }
    x_id	<- apply( X1, 2, sd ) != 0
    nx0	  <- sum( x_id )
    if( nx0 == 0 ){
      X1	<- NULL
      xname	<- NULL
      x_id	<- NULL
    } else {
      X1	<- as.matrix( X1[ , x_id ] )
      xname	<- names( as.data.frame( X1 ) )
    }

    if( allsvc_flag == 1 ){
      allsvc_id        <- rep(FALSE,length(x_id))
      allsvc_id[allsvc]<- TRUE
      allsvc_id        <- allsvc_id[ x_id ]
    }

    if( allnvc_x_flag == 1 ){
      allnvc_x_id        <- rep(FALSE,length(x_id))
      allnvc_x_id[allnvc_x]<- TRUE
      allnvc_x_id        <- allnvc_x_id[ x_id ]
    }

  } else {
    xname	<- NULL
    x_id	<- NULL
  }

  Xconst	<- xconst
  if( is.null( xconst ) == FALSE ){
    Xconst	<- as.matrix( Xconst )
    if( is.numeric( Xconst ) == FALSE ){
      mode( Xconst ) <- "numeric"
    }
    xf_id	<- apply( Xconst, 2, sd ) != 0
    nxf	<- sum( xf_id )
    if( nxf == 0 ){
      Xconst	<- NULL
      xfname	<- NULL
      xf_id	  <- NULL
    } else {
      Xconst	<- as.matrix( Xconst[ , xf_id ] )
      xfname	<- names( as.data.frame( Xconst ) )
    }

    if( allnvc_xconst_flag == 1 ){
      allnvc_xconst_id               <- rep(FALSE,length(xf_id))
      allnvc_xconst_id[allnvc_xconst]<- TRUE
      allnvc_xconst_id  <- allnvc_xconst_id[ xf_id ]
    }
  } else {
    xfname	<- NULL
    xf_id	  <- NULL
    nxf	    <- 0
  }

  #XX1_0	      <- list(NULL)
  #XX1	        <- NULL
  #nvc_x_use   <- apply( X1, 2, function( x ) length( unique( x ))) >= 3
  if( is.logical( nvc_x[ 1 ] )&( nvc_x[ 1 ] == TRUE) ) nvc_x     <- 1:nx0

  pmod_x   <- list(NULL)
  nsamp    <- NULL
  if( is.logical( nvc_x[ 1 ] ) == FALSE ){
    X1_nvc  <- as.matrix(X1[ , nvc_x ])
    xxname	<- names( as.data.frame( X1 ) )[ nvc_x ]
    nnsv    <- length( xxname )
    np_xx     <-apply( X1_nvc,2,function( x ) length( unique( x )))
    np_xx     <-ifelse( np_xx < nvc_num/0.7, round( np_xx * 0.7 ) ,nvc_num )#chganged

    for( ii in 1:dim( X1_nvc )[ 2 ] ){
      test    <- TRUE
      iiii    <- 0

      while(test){
        kkk      <- np_xx[ ii ]-iiii
        knots    <-seq(min( X1_nvc[ , ii ] ),max( X1_nvc[ , ii ] ),len=kkk+2)[2:(kkk+1)]
        #if( n > 2000000 ){############################# to be implemented #################
        #  if( !is.null(nsamp) ) nsamp   <- sample( 2000000 )
        #  testt<-try(XX1_00_b<- ns( X1_nvc[ nsamp,ii ], knots=knots ), silent=TRUE)
        #} else {
        testt<-try(XX1_00<- ns( X1_nvc[ ,ii ], knots = knots ), silent=TRUE)
        #}
        test <- inherits(testt, "try-error")#class(testt)[1] == "try-error"
        iiii <- iiii+1
      }

      #if( n > 2000000 ){
      #  testt2<-try(XX1_00     <- sel_basis_vif( cbind( X1_nvc[nsamp,ii], XX1_00_b), vif_max=15),silent=TRUE)
      #} else{
      testt2<-try(XX1_00     <- sel_basis_vif( cbind( X1_nvc[,ii], XX1_00), vif_max=15),silent=TRUE )
      #}

      if( !inherits( testt2, "try-error" ) ){#class(testt2)[1] != "try-error"
        np_xx[ ii ]   <- dim(XX1_00)[2]
      } else {
        np_xx[ ii ]   <- 0
      }
      if( np_xx[ ii ] > 2 ){
        pmod_x[[ii]]<- scale( XX1_00 )
      } else {
        pmod_x[[ii]]<- 0
        np_xx[ ii ] <- 0
      }
    }

  } else {
    nnsv      <- 0
    np_xx     <- NULL
  }

  if( is.logical( nvc_xconst[ 1 ] )&( nvc_xconst[ 1 ] == TRUE) ) nvc_xconst<- 1:nxf

  pmod_xconst <- list(NULL)
  if( is.logical( nvc_xconst[ 1 ] ) ==FALSE ){
    Xconst_nvc<- as.matrix( Xconst[ , nvc_xconst ] )
    xxfname	  <- names( as.data.frame( Xconst ) )[ nvc_xconst ]
    nnxf      <- length( xxfname )
    np_xxconst       <-apply( Xconst_nvc,2,function( x ) length( unique( x )))
    np_xxconst       <-ifelse( np_xxconst < nvc_num/0.7, round( np_xxconst * 0.7 ) ,nvc_num )

    rem_ord_XXc<-list(NULL)
    for( ii in 1:dim( Xconst_nvc )[ 2 ] ){
      test<-TRUE
      iiii<-0

      while(test){
        kkk      <- np_xxconst[ ii ]-iiii
        knots    <-seq(min( Xconst_nvc[ ,ii ] ),max( Xconst_nvc[ ,ii ] ),len=kkk+2)[2:(kkk+1)]
        #if( n > 2000000 ){############################# to be implemented #################
        #  if( !is.null(nsamp) ) nsamp <- sample( 2000000 )
        #  testt<-try(XXconst_00_b<- ns( Xconst_nvc[ nsamp,ii ], knots = knots ), silent=TRUE)
        #} else {
        testt<-try(XXconst_00<- ns( Xconst_nvc[ ,ii ], knots = knots ), silent=TRUE)
        #}
        test <- inherits(testt, "try-error")#class(testt)[1] == "try-error"
        iiii <- iiii+1
      }

      #if( n > 2000000 ){
      #  testt2<-try(XXconst_00     <- sel_basis_vif( cbind( Xconst_nvc[nsamp,ii], XXconst_00_b ), vif_max=15),silent=TRUE)
      #} else{
      testt2<-try(XXconst_00     <- sel_basis_vif( cbind( Xconst_nvc[,ii], XXconst_00 ), vif_max=15),silent=TRUE)
      #}

      if( !inherits(testt2, "try-error") ){#class(testt2)[1] != "try-error"
        np_xxconst[ ii ]  <- dim(XXconst_00)[2]
      }else{
        np_xxconst[ ii ] <- 0
      }

      if( np_xxconst[ ii ] > 2 ){
        pmod_xconst[[ii]]<- scale( XXconst_00 )
      } else {
        pmod_xconst[[ii]]<- 0
        np_xxconst[ ii ] <- 0
      }
    }

  } else {
    nnxf         <- 0
    np_xxconst   <- NULL
  }

  if( ( nx0 > 0 ) & ( nxf > 0 ) ){
    dup      <-duplicated(cbind(X1,Xconst),MARGIN=2)
    if(sum(dup)>0){
      stop( "x and xconst cannot have the same column" )
    }
  }

  #X0	<- as.matrix( cbind( 1, Xconst, X1 ) )
  #Mo 	<- crossprod( X0 )
  #mo	<- crossprod( X0, y )
  #parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )

  nsv	<- ifelse( is.null( X1 ), 1, dim( X1 )[ 2 ] + 1 )

  if( !is.null( s_id )[ 1 ] ){
    coords_gx<-tapply(coords[,1],s_id,mean)
    coords_gy<-tapply(coords[,2],s_id,mean)
    coords_g <-as.matrix(cbind(coords_gx, coords_gy))
    #s_id2   <-data.frame( s_id=rownames(coords_x), s_id_num = 1:length(coords_x) )
    #s_id_dat<-data.frame( s_id = s_id, id = 1:length(s_id))
    #s_id_dat2<-merge(s_id_dat, s_id2, by="s_id", all.x=TRUE)
    #s_id_dat2<-s_id_dat2[order(s_id_dat2$id),c("s_id", "s_id_num")]
    enum    <- min( enum, length( coords_g ) - 1 )
    suppressWarnings(coordk	<- kmeans( coords, centers = enum + 1 )$centers)
  } else {
    suppressWarnings(coordk	<- kmeans( coords, centers = enum + 1 )$centers)
  }

  D	      <- rdist( coordk )
  h	      <- max( spantree( D )$dist )
  if( covmodel == "exp" ){
    C	    <- exp( -D / h )
  } else if( covmodel == "gau" ){
    C	    <- exp( -( D / h) ^ 2 )
  } else if( covmodel == "sph" ){
    C	    <- ifelse( D < h , 1 - 1.5 * (D / h ) + 0.5 * ( D / h ) ^ 3, 0 )
  } else {
    stop( "covmodel is not specified appropriately" )
  }

  Cmean	  <- apply( C, 1, mean )
  MCM	    <- t( C - Cmean ) - Cmean + mean( Cmean )
  eigenC	<- eigen( MCM )
  sf0	    <- eigenC$vectors[ , eigenC$values > 1e-08 ]
  ev0	    <- eigenC$values [   eigenC$values > 1e-08 ]
  sf_ap0	<- sf0 %*% diag( 1 / ev0 )
  ev_full	<- ev0 * ( n / enum ) - 1
  ev	    <- ev_full[   ev_full > 1e-08 ]
  nev		<- length( ev )

  if( is.null( Xconst )&is.null( X1 )){
    X0	<- as.matrix( rep(1,n) )
  } else {
    X0	<- as.matrix( cbind( 1, Xconst, X1 ) )
  }

  #X0	<- as.matrix( cbind( 1, Xconst, X1 ) )
  nx	<- dim( X0 )[ 2 ]
  Mo 	<- crossprod( X0 )
  mo	<- crossprod( X0, y )
  parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )

  nblock  <-max(1,round((n*(nxf + length(ev)*nsv)/bsize^2)))
  ids   <- unique(c(round(seq(1,n,len=nblock+1)),n))

  if(is.null(ncores)) {
    ncores <- makeCluster(detectCores(),setup_strategy = "sequential")
  } else {
    ncores <- makeCluster(ncores,setup_strategy = "sequential")
  }
  registerDoParallel(ncores)

  im     <- NULL
  sfsf   <- foreach(im = 1:(length(ids)-1), .packages=c("fields", "splines"),.combine="+") %dopar%  {
    if(im==1){
      ids_sub<-ids[im]:ids[im+1]
    } else {
      ids_sub<-(ids[im]+1):ids[im+1]
    }

    if( covmodel == "exp" ){
      ccc0	<- exp( -rdist( coordk, coords[ids_sub,] ) / h )
    } else if( covmodel == "gau" ){
      C	    <- exp( -( D / h) ^ 2 )
      ccc0	<- exp( -( rdist( coordk, coords[ids_sub,] ) / h )^2 )
    } else if( covmodel == "sph" ){
      ddd0  <- rdist( coordk, coords[ids_sub,] )
      ccc0  <- ifelse( ddd0 < h , 1 - 1.5 * ( ddd0 / h ) + 0.5 * ( ddd0 / h ) ^ 3, 0 )
    }
    sf  <- t( ccc0 - Cmean ) %*% sf_ap0[ , ev_full > 1e-08 ]

    X2   <- cbind(X0[ids_sub,  ], sf)
    cSum_sub0<- rep( 0, nx )

    cSum_sub1 <- colSums( sf )
    if( nsv >= 2 ){
      for( im2 in 1:( nsv - 1 ) ){
        X2       <- cbind( X2, X1[ids_sub, im2 ] * sf )
        cSum_sub1 <- c( cSum_sub1, colSums( sf ) )
      }
    }

    cSum_sub2 <- NULL
    if( is.logical( nvc_xconst[ 1 ] ) == FALSE ){
      for( ii in 1:length( np_xxconst )){
        if(np_xxconst[ii]>0){
          basis  <- as.matrix( pmod_xconst[[ ii ]][ids_sub, ] )
          X2     <- cbind( X2, Xconst_nvc[ids_sub, ii ]*basis )
          cSum_sub2 <- c( cSum_sub2, colSums( basis ) )
        }
      }
    }

    cSum_sub3 <- NULL
    if( is.logical( nvc_x[ 1 ] ) == FALSE ){
      for( ii in 1:length( np_xx ) ){
        if(np_xx[ii]>0){
          basis  <- as.matrix( pmod_x[[ii]][ids_sub,  ] )
          X2     <- cbind( X2, X1_nvc[ids_sub, ii ]*basis )
          cSum_sub3 <- c( cSum_sub3, colSums( basis ) )
        }
      }
    }

    cSum <- c(cSum_sub0, cSum_sub1, cSum_sub2, cSum_sub3)
    sfsf1<-cbind(cSum, crossprod( X2, y[ids_sub] ),crossprod( X2 ))
    #sfsf1<-cbind(rep(0,length(sfsf0[,1])),sfsf0)

    #clsf<-colSums(sf)
    #sfsf1[1:length(clsf),1]<-clsf

    sfsf1
  }

  sfmean <-sfsf[,1]/n
  m0     <-c(sfsf[,2])
  EX     <-sfsf[, -(1:2)][ ,1:nx ]

  id	   <- c( rep( 0, nx ), rep( 1:nsv, each = nev ) )
  sf_id  <- c(rep(1,nx),rep(1, nev))
  sf_id2 <- 1
  sf_id3 <- 1
  if( nsv >=2 ){
    sf_id <-c( sf_id, nxf + rep(2:nsv,each=nev) )
    sf_id2 <- c(sf_id2, nxf + 2:nsv)
    sf_id3 <- c(sf_id3, nxf + 2:nsv)
  }

  ii_count <- max( id )
  if( is.logical( nvc_xconst[ 1 ] ) == FALSE ){
    for( ii in 1:length( np_xxconst ) ){
      ii_count <- ii_count + 1
      if(np_xxconst[ii]>0){
        id	   <- c( id, rep( ii_count, np_xxconst[ii] ))
        sf_id  <- c(sf_id, rep( 1+ii, np_xxconst[ii] ) )
        sf_id2 <- c( sf_id2, 1+ii )
      } else {
        sf_id2 <- c( sf_id2, NA )
      }
      sf_id3   <- c( sf_id3, 1+ii )
    }
  }

  if( is.logical( nvc_x[ 1 ] ) == FALSE ){
    for( ii in 1:length( np_xx ) ){
      ii_count <- ii_count + 1
      if(np_xx[ii]>0){
        id	   <- c( id, rep( ii_count, np_xx[ii] ))
        sf_id  <-c(sf_id, rep( 1+nxf+ii, np_xx[ii] ) )
        sf_id2 <- c( sf_id2, 1+nxf+ii )
      } else {
        sf_id2 <- c( sf_id2, NA )
      }
      sf_id3   <- c( sf_id3, 1+nxf+ii )
    }
  }

  XX2    <-Mo[sf_id,sf_id]*outer(sfmean,sfmean)
  EX2    <-t(t(as.matrix(EX)[,sf_id])*sfmean)
  M      <-sfsf[, -(1:2)] - EX2 - t(EX2) + XX2
  m      <-m0 - sfmean*mo[sf_id]
  yy     <- sum( y ^ 2 )

  if( penalty == "aic" ){
    pen	<- 2
  } else if( penalty == "bic" ){
    pen	<- log( n )
  }
  par0	<- rep( 1, 2 * nsv )
  par0[1] <- parVmax_sq / 3
  par0_est<- c( 1, 1 )
  par0_id	<- rep( 1:nsv, 2 )
  null_dum<- rep( 0, nsv )

  #if( is.null( xgroup ) == FALSE ){
  #  par0	<- c( par0, rep( parVmax_sq/3, ng ) )## changed
  #  par0_id	<- c( par0_id , max(par0_id) + ( 1 : ng ) )
  #  null_dum<- c( null_dum, rep( 0, ng ) )
  #}

  #pp_id    <-max( id )
  if( is.logical( nvc_xconst[ 1 ] ) == FALSE ){#nvc_xconst[ 1 ] != FALSE
    par0	  <- c( par0, rep( 1, nnxf ) )
    par0_id	<- c( par0_id , max(par0_id) + ( 1:nnxf ) )
    null_dum0<- rep( 0, nnxf ) # added
    null_dum0[np_xxconst==0]<-1# added
    null_dum<- c( null_dum, null_dum0 )
  }
  if( is.logical( nvc_x[ 1 ] ) == FALSE ){#nvc_x[ 1 ] != FALSE
    par0	<- c( par0, rep( 1, nnsv ) )
    par0_id	<- c( par0_id , max(par0_id) + ( 1:nnsv ) )
    null_dum0<- rep( 0, nnsv ) # added
    null_dum0[np_xx==0]<-1     # added
    null_dum<- c( null_dum, null_dum0 )
  }

  id_exist<-unique(id)[unique(id)>0]
  par0[(par0_id %in% id_exist) ==FALSE]<-0
  Par    <-par0

  obj	   <- sd( y )
  LL	   <- NULL
  res_old<- Inf
  n_omit <- 0
  iter	 <- 1
  gmess  <- 1
  warn   <- FALSE
  while( (( obj > sd( y ) * 0.0001 ) & ( iter <= maxiter ))|( iter <= 4 ) ){

    if( !is.null(x) ) print( paste( "-------  Iteration ", iter, "  -------", sep = "" ) )

    LL0	<- LL
    n_omit	<- 0
    for( par0_sel in (1:( nsv + nnxf+ ng + nnsv ))[id_exist]){
      evSqrt	<- NULL
      par0_sq	<- par0 ^ 2
      for( i in 1:nsv ){
        evv  	  <- ev ^ par0_sq[ nsv + i ] * sum( ev ) / sum( ev ^ par0_sq[ nsv + i ] )
        evSqrt	<- c( evSqrt, par0_sq[ i ] * sqrt( evv ) )
      }

      M0	<- M
      for( j in 1:nsv ){
        if( j != par0_sel ){
          id_sub	<- ( id == j )
          diag( M0[ id_sub, id_sub ] )<-
            diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == j ] ^ 2
        }
      }

      if( ng != 0 ){
        for( i2 in 1:ng ){
          xgg	<- rep( 1, sum( id == ( nsv + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ 2 * nsv + i2 ] * xgg )
        }
        for( j2 in 1:ng ){
          if( j2 != ( par0_sel - nsv ) ){
            id_sub	<- ( id == ( j2 + nsv ) )
            diag( M0[ id_sub, id_sub ] )<-
              diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
          } else {############################### added
            id_sub	<- ( id == ( j2 + nsv ) )
            g_add   <- (1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2)/100
            diag( M0[ id_sub, id_sub ] )<-
              diag( M0[ id_sub, id_sub ] ) + g_add
          }       ############################### added
        }
      }

      if( nnxf != 0 ){
        for( i2 in 1:nnxf ){
          xgg	<- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng ) + i2 ] * xgg )
        }

        for( j2 in 1:nnxf ){
          if( j2 != ( par0_sel - nsv - ng ) ){
            id_sub	<- ( id == ( j2 + nsv + ng ) )
            diag( M0[ id_sub, id_sub ] )<-
              diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng) ] ^ 2
          }
        }
      }

      if( nnsv != 0 ){
        for( i2 in 1:nnsv ){
          xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng + nnxf ) + i2 ] * xgg )
        }

        for( j2 in 1:nnsv ){
          if( j2 != ( par0_sel - nsv - ng - nnxf ) ){
            id_sub	<- ( id == ( j2 + nsv + ng + nnxf) )
            diag( M0[ id_sub, id_sub ] )<-
            diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng + nnxf) ] ^ 2
          }
        }
      }

      null_dum2<- null_dum
      null_dum2[ par0_sel ] <- 0
      null_dum3<- c( rep( 0, nx ), null_dum2[ id ]) == 0
      MM	<- M [ null_dum3, null_dum3 ]
      MM0	<- M0[ null_dum3, null_dum3 ]
      mm	<- m [ null_dum3 ]
      idd	<- id[ null_dum3 ]
      id_omit1<- c( diff( id ), 1)
      #if( ng != 0 )	id_omit1[ id > nsv ] <- 0
      #if( ng != 0 )	id_omit1[( nsv < id ) & ( id <= nsv + ng )] <- 0
      id_omit2_0<- which( id_omit1 != 0 )                           ### changed
      id_omit2<-id_omit2_0[2:(nsv + 1)]                             ### changed
      id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0### changed

      sstop <-FALSE
      if( ( ng != 0 )&( nsv < par0_sel )&( par0_sel <= nsv + ng ) ){
        if( rcond( MM0 ) < 1e-30 ){
          loglik	<- ( - 1 / 2 ) * res_old
          par0[ par0_id == par0_sel ] 	<- 0
          null_dum[ par0_sel ]		<- 1
          loglik 	<- ( -1 / 2 ) * c( res_ref )
          res_old	<- res_ref
          score   <- res_ref+ pen * (nx+sum(par0!=0)+1)

          #if(gmess == 1 ){
          message( paste( "Note:", "group effect", par0_sel - nsv, "is omitted to stablize the estimates", sep = " " ) )
          #}
          gmess   <- gmess + 1
          sstop   <-TRUE
        }
      }

      if( sstop==FALSE ){

        er_dum	<- TRUE
        n_omit0	<- 0
        while( er_dum == TRUE ){
          try1	<- try( M0inv	<- solve( MM0,  tol = 1e-30 ), silent = TRUE )
          try2	<- try( Mdet0	<- Mdet_f0( M = MM, M0 = MM0, id = idd,
                                         par0_sel = par0_sel, emet = method ), silent = TRUE)
          #er_dum	<- ( class(try1)[1] =="try-error" ) | ( class(try2)[1] =="try-error" )
          er_dum	<- inherits( try1, "try-error" )|inherits( try2, "try-error" )

          if( er_dum == TRUE ){
            M	  <- M [ id_omit1 == 0, id_omit1 == 0 ]
            M0	<- M0[ id_omit1 == 0, id_omit1 == 0 ]
            m	  <- m [ id_omit1 == 0 ]
            id	<- id[ id_omit1 == 0 ]
            #X	<- X [, id_omit1 == 0 ] ### difference
            null_dum3<- null_dum3[ id_omit1 == 0 ]

            MM	<- M [ null_dum3, null_dum3 ]
            MM0	<- M0[ null_dum3, null_dum3 ]
            mm	<- m [ null_dum3 ]
            idd	<- id[ null_dum3 ]

            id_omit1<- c( diff( id ), 1)
            #if( ng != 0 )	id_omit1[ id > nsv ] <- 0
            #if( ng != 0 )	id_omit1[( nsv < id ) & ( id <= nsv + ng )] <- 0 ## deleted

            id_omit2_0<- which( id_omit1 != 0 )                           ### changed
            id_omit2<-id_omit2_0[2:(nsv + 1)]                             ### changed
            id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0### changed
            n_omit0	<- n_omit0 + 1
            nev     <-nev-1   ## added
            warn    <- TRUE
          }

          if( n_omit0 >= 1000){
            stop( " Singular fit. Simplify the model" )
          }
        }

        if( n_omit0 > 0 ){
          message( paste( "Note: ", n_omit0, " eigenvectors are omitted to stablize the estimates (dummy variable in x is a typical cause)", sep = "" ) )
        }
        n_omit	<- c( n_omit, n_omit0 )

        ev	   <- ev     [   1:sum( id == 1 ) ]
        #meig$sf	<- meig$sf[ , 1:sum( id == 1 ) ] ##################3 changed
        term2	  <- Mdet0$term2
        term3_0	<- Mdet0$term3_0

        if( min( par0[ par0_id == par0_sel ] ) >= 1e-5 ){
          par00	<- par0[ par0_id == par0_sel ]
        } else {
          if(par0_sel<=nsv){
            par00_id<- max( which( c(apply(Par[,par0_id == par0_sel], 1, min)) != 0 ) )#fixed
          } else {
            par00_id<- max( which( Par[,par0_id == par0_sel] != 0 ) )## fixed
          }

          par00	<- Par[ par00_id, par0_id == par0_sel ]
        }

        if( n_omit0 > 0 ){
          res_old <- lik_resf_vc0( par0, ev = ev, M = MM, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                                   n = n, nx = nx, nsv = nsv, ng = ng, emet = method, null_dum4 = null_dum2 )
        }

        M0inv_01<- M0inv[ 		, idd == par0_sel ]
        M0inv_00<- M0inv[ idd ==par0_sel, idd == par0_sel ]
        b_01	<- M0inv %*% mm
        b_02	<- t( M0inv_01 ) %*% mm

        if( par0_sel == 1){
          llim    <- c( parVmax_sq / 1000, 1e-03)
          ulim    <- c( parVmax_sq, 4 )
          omethod <- "L-BFGS-B"##### error length of b[-(1:nx)] and evSqrt are inconsistent

          #M = MM; M0inv = M0inv; m = mm; emet = method; id = idd
          #lower = llim; upper = ulim; method = omethod
          res    	<- optim( fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                            M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                            term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                            n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                            par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2,
                            lower = llim, upper = ulim, method = omethod )
          res_int <- res

        } else if( par0_sel <= nsv ){
          res    	<- optim( fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                            M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                            term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                            n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                            par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2 )

        } else if( par0_sel <= nsv+ng ){#### changed
          res    	<- optimize( f = lik_resf_vc, lower = parVmax_sq/1000, upper = parVmax_sq, par00,
                               par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                               M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                               term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                               n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                               par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2 )
          res$value	<- c( res$objective )
          res$par	<- c( res$minimum )

        } else {
          res    	<- optimize( f = lik_resf_vc, lower = parVmax_sq/10^5, upper = 10^10, par00, par0 = par0, ev = ev,
                               M = MM, M0inv = M0inv,
                               M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                               term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                               n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                               par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2 )
          res$value	<- c( res$objective )
          res$par	<- c( res$minimum )
        }

        if( ( iter > 3 ) & ( res$value > res_old ) ){
          loglik	<- ( - 1 / 2 ) * res_old
        } else {
          if( par0_sel != 1 ){#### ) & ( par0_sel <= nsv ))|( par0_sel > nsv + ng )
            MM_ref 	<- MM [ idd != par0_sel, idd != par0_sel ]
            mm_ref	<- mm [ idd != par0_sel ]
            null_dum4<- null_dum2
            null_dum4[ par0_sel ] <- 1
            res_ref <- lik_resf_vc0( par0, ev = ev, M = MM_ref, m = mm_ref, yy = yy,nnxf=nnxf,nnsv=nnsv,
                                     n = n, nx = nx, nsv = nsv, ng = ng, emet = method, null_dum4 = null_dum4 )
          } else {
            res_ref	<- Inf
          }

          np_add       <- length(par00)
          flag         <- 0
          if( par0_sel <= nsv ){
            if( allsvc_flag == 0 ){   ####### added
              allvc  <- allsvc
            } else {
              flag   <-1
              allvc  <- c(TRUE, allsvc_id)[ par0_sel ] # including intercept
            }
          } else if( ( nsv < par0_sel )&( par0_sel <= nsv + ng )){
            allvc  <- TRUE#allg
          } else if( ( nsv +ng < par0_sel )&( par0_sel <= nsv + ng + nnxf )){
            if( allnvc_xconst_flag == 0 ){   ####### added
              allvc  <- allnvc_xconst
            } else {
              flag   <-1
              allvc  <- allnvc_xconst_id[ par0_sel - nsv - ng ]
            }
          } else {
            if( allnvc_x_flag == 0 ){   ####### added
              allvc  <- allnvc_x
            } else {
              flag   <-1
              allvc  <- allnvc_x_id[ par0_sel - nsv - ng - nnxf ]
            }
          }

          if( ((flag==0)&(res_ref < res$value + pen * np_add)&(allvc[1]==FALSE))|((flag==1)&(allvc[1]==FALSE))){### reject
            par0[ par0_id == par0_sel ] 	<- 0
            null_dum[ par0_sel ]		<- 1
            loglik 	<- ( -1 / 2 ) * c( res_ref )
            res_old	<- res_ref
            score   <- res_ref+ pen * (nx+sum(par0!=0)+1)
          } else {      ####### accept
            par0[ par0_id == par0_sel ]	<- res$par
            null_dum[ par0_sel ]		<- 0
            res_old	<- res$value
            loglik 	<- ( -1 / 2 ) * res$value
            score   <- res_old + pen * (nx+sum(par0!=0)+1)
          }
        }
      }

      LL0	<- c( LL0, loglik )
      if( !is.null(x) ) print( paste( par0_sel, "/", nsv + ng + nnxf + nnsv, sep = "" ) )
    }

    if( iter > 1 ){
      if( sum( n_omit ) == 0 ){
        obj 	<- abs( loglik - loglik0 )
      } else {
        obj	<- Inf
      }
    }

    Par	    <- rbind( Par, par0 )
    loglik0	<- loglik
    LL	    <- c( LL, loglik )
    if( !is.null(x) ) print( paste( toupper( penalty ), ": ", round( score, 3 ), sep = "" ) )
    iter	  <- iter + 1
  }

  par2   	<- par0 ^ 2
  evSqrt	<- NULL
  evSqrt2	<- list(NULL)
  for( i in 1:nsv ){
    evv  	<- ev ^ par2[ nsv + i ] * sum( ev ) / sum( ev ^ par2[ nsv + i ] )
    evSqrt<- c( evSqrt, par2[ i ] * sqrt( evv ) )
    evSqrt2[[ i ]]	<- par2[ i ] * sqrt( evv )
  }

  if( ng != 0 ){
    for( i2 in 1:ng ){
      xgg	<- rep( 1, sum( id == ( nsv + i2 ) ) )
      evSqrt	<- c( evSqrt, par2[ 2 * nsv + i2 ] * xgg )
      evSqrt2[[ nsv + i2 ]]	<- par2[ 2 * nsv + i2 ] * xgg
    }
  }

  if( nnxf != 0 ){
    for( i2 in 1:nnxf ){
        xgg	<- rep( 1, sum( id == ( nsv + ng + i2 ) ) )
        evSqrt	<- c( evSqrt, par2[ 2 * nsv + ng + i2 ] * xgg )
        evSqrt2[[ nsv + ng + i2 ]]	<- par2[ 2 * nsv + ng + i2 ] * xgg
    }
  }

  if( nnsv != 0 ){
    for( i2 in 1:nnsv ){
      xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
      evSqrt	<- c( evSqrt, par2[ 2 * nsv + ng + nnxf + i2 ] * xgg )
      evSqrt2[[ nsv + ng + nnxf + i2 ]]	<- par2[ 2 * nsv + ng + nnxf + i2 ] * xgg
    }
  }


  MM	<- M [ null_dum3, null_dum3 ]
  MM0	<- M0[ null_dum3, null_dum3 ]
  mm	<- m [ null_dum3 ]
  idd	<- id[ null_dum3 ]

  M[ -( 1:nx ), -( 1:nx ) ]	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
  M[ -( 1:nx ),    1:nx   ]	<-    M[ -( 1:nx ),    1:nx   ] * evSqrt
  M[    1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ),    1:nx   ] )
  diag( M[ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1

  MM		<- M[ null_dum3, null_dum3 ]
  MMinv		<- solve( MM, tol = 1e-30 )

  mm		<- m[ null_dum3 ]
  eevSqrt		<- evSqrt[ null_dum3[ -( 1:nx )] ]
  mm[ -( 1:nx ) ]	<- mm[ -( 1:nx ) ] * eevSqrt
  b		<- MMinv %*% mm
  b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * eevSqrt

  null_dum3_id<- which(null_dum3)
  sf_id4    <- c(nx,rep(nev,nsv),np_xxconst,np_xx) #fixed id_omit1
  sf_id5    <- 1:nsv
  eff_id    <- rep(1,nsv)
  if( nnxf != 0) {
    sf_id5 <- c( sf_id5, nsv + 1:nnxf )
    eff_id <- c( eff_id, rep(2, nnxf))
  }
  if( nnsv != 0){
    sf_id5 <- c( sf_id5, 2:nsv )#not nsv
    eff_id <- c( eff_id, rep(3, nnsv))#here is nnsv
  }
  sec_id   <- 1:( nsv+ ng + nnxf + nnsv )

  #b_s		    <- list(NULL)
  #b_covs		<- list(NULL)
  bvc_all   <- foreach(i = 1:(length(ids)-1), .packages=c("fields","splines"),.combine="rbind") %dopar%  {
    if(i==1){
      ids_sub<-ids[i]:ids[i+1]
    } else {
      ids_sub<-(ids[i]+1):ids[i+1]
    }
    n_sub       <- length(min(ids_sub):max(ids_sub) )

    b_vc_s0     <- matrix(0, nrow= length(ids_sub), ncol= nsv );j_s<-1
    b_vc_n0     <- matrix(0, nrow= length(ids_sub), ncol= nsv );j_n<-1
    bse_vc_s0   <- matrix(0, nrow= length(ids_sub), ncol= nsv )
    bse_vc_n0   <- matrix(0, nrow= length(ids_sub), ncol= nsv )
    b_vc0       <- matrix(0, nrow= length(ids_sub), ncol= nsv )
    bse_vc0     <- matrix(0, nrow= length(ids_sub), ncol= nsv )
    if(nnxf>0) {
      c_vc0     <- matrix(0, nrow= length(ids_sub), ncol= nnxf );j_c<-1
      cse_vc0   <- matrix(0, nrow= length(ids_sub), ncol= nnxf )
    } else {           ## added
      c_vc0     <-NULL
      cse_vc0   <-NULL
    }

    for( j in unique( sf_id5 ) ){#( nsv+ ng + nnxf + nnsv )
      bid0_NA     <- sf_id2[ sf_id5 == j ]#ifelse( j == 1, 1, nxf + j )
      bid0		    <- sf_id3[ sf_id5 == j ]#ifelse( j == 1, 1, nxf + j )

      null_dum_sub<- null_dum[ sf_id5 == j ]
      b_full_id   <- cumsum( sf_id4[ c(0, null_dum) == 0] )

      if( (j %in% sf_id5[null_dum==0])==FALSE ){### if vc is not selected
        if( j <= nsv ){

          b_vc0  [,j]<- b_vc_s0  [,j]<- b_vc_n0  [,j]<- b[bid0[1]] ##### added
          bse_vc0[,j]<- bse_vc_s0[,j]<- bse_vc_n0[,j]<- c(sqrt(MMinv[bid0[1], bid0[1]]))#### added

          #b_vc0[ ,j ]  <- b[ bid0[ 1 ] ]
          #bse_vc0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          j_s          <- j_s + 1
          j_n          <- j_n + 1
        } else {
          c_vc0[ ,j_c ]  <- b[ bid0[ 1 ] ]
          cse_vc0[ ,j_c ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          j_c          <- j_c + 1
        }
      } else {
        sf_id4_start<- b_full_id[ which(sf_id5[null_dum==0] == j) ] + 1
        sf_id4_end  <- b_full_id[ which(sf_id5[null_dum==0] == j) + 1 ]
        eff_id_sub  <- eff_id[null_dum==0][ sf_id5[null_dum==0] == j ]#(1,1,1,3,3)#(1=snvc, 3=nvc)
        sec_id_sub  <- sec_id[null_dum==0][ sf_id5[null_dum==0] == j ]#(1,2,3,4,5)
        bid         <- list( NULL )
        sfmean_sub  <- list(NULL)
        bid[[ 1 ]]  <- bid0[1]

        bid_len     <- length( sf_id4_start )
        for(jjj in 1:bid_len){
          sell              <- ifelse(eff_id_sub[jjj]==1, 2, 3)     #(1=snvc, 3=nvc)
          bid[[sell]]       <- sf_id4_start[ jjj ]:sf_id4_end[ jjj ]#selct column
          sfmean_sub[[sell]]<-sfmean[ bid[[sell]] ]
        }

        if( 1 %in% eff_id_sub ){
          if( (null_dum_sub[ 1 ] == 0)&( !is.na( bid0_NA[ 1 ] ) ) ){
            if( covmodel == "exp" ){
              ccc0	<- exp( -rdist( coordk, coords[ids_sub,] ) / h )
            } else if( covmodel == "gau" ){
              C	    <- exp( -( D / h) ^ 2 )
              ccc0	<- exp( -( rdist( coordk, coords[ids_sub,] ) / h )^2 )
            } else if( covmodel == "sph" ){
              ddd0  <- rdist( coordk, coords[ids_sub,] )
              ccc0  <- ifelse( ddd0 < h , 1 - 1.5 * ( ddd0 / h ) + 0.5 * ( ddd0 / h ) ^ 3, 0 )
            }
            basis_s <- t( ccc0 - Cmean ) %*% sf_ap0[ , 1:nev ]
            basis_s <- t( t( basis_s ) - sfmean_sub[[ 2 ]] )
            b_vc_s0[ ,j_s ]<- basis_s %*% b[ bid[[ 2 ]] ]
            MMinv_sub_s   <- MMinv[ unlist(bid[ 1:2 ]), unlist(bid[ 1:2 ])] ## Added below #########
            basis_sel_s   <- evSqrt2[which(sf_id5==j)][[1]]#[!is.null(basis_s)]
            basis2_s		  <- t( t( basis_s ) * unlist( basis_sel_s ))## added
            x_basis_s	    <- as.matrix( cbind( 1, basis2_s ) )###added
            bse_vc_s0[ ,j ]<- sqrt( colSums( t( x_basis_s ) * ( MMinv_sub_s %*% t( x_basis_s ))))### Added above ##

          } else {
            basis_s      <- NULL
            bse_vc_s0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          }
          j_s       <- j_s + 1
        } else if(1 %in% eff_id[ sf_id5==j ]){
          basis_s   <- NULL
          bse_vc_s0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          j_s       <- j_s + 1
        }

        if( 2 %in% eff_id_sub ){
          if( (null_dum_sub[ 1 ] == 0)&( !is.na( bid0_NA[ 1 ] ) ) ){
            jj     <- sec_id_sub[ 1 ] - nsv - ng
            basis_c<- as.matrix( pmod_xconst[[ jj ]][ids_sub,  ] )
            basis_c<- t( t(basis_c) - sfmean_sub[[ 3 ]])
            c_vc0[ ,j_c ]  <- b[ bid0[ 1 ] ] + basis_c %*% b[ bid[[ 3 ]] ]

            MMinv_sub    <-MMinv[ unlist(bid), unlist(bid)]
            basis2		   <- t( t( basis_c ) * unlist( evSqrt2[which(sf_id5==j)] ) )
            x_basis		   <- as.matrix( cbind( 1, basis2 ) )###added
            cse_vc0[ ,j_c ]<- sqrt( colSums( t( x_basis ) * ( MMinv_sub %*% t( x_basis ))))

          } else {
            basis_c       <- NULL
            c_vc0[ ,j_c ]  <- c(b[ bid0[ 1 ] ])#rep( b[ bid0[ 1 ] ],n_sub )
            cse_vc0[ ,j_c ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
                            #rep(sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ),n_sub)
          }
          j_c       <- j_c + 1
        }

        if( 3 %in% eff_id_sub ){
          if( (null_dum_sub[ 2 ] == 0)&( !is.na( bid0_NA[ 2 ] ) ) ){
            jj     <- sec_id_sub[ length(sec_id_sub) ] - nsv - ng - nnxf
            basis_n<- as.matrix( pmod_x[[ jj ]][ids_sub,  ] )
            #basis_n<- as.matrix( predict(pmod_x[[ jj ]],X1_nvc[ids_sub ,jj ]) )#[ ,-1 ]
            basis_n<- t( t(basis_n) - sfmean_sub[[ 3 ]])
            b_vc_n0[ ,j_n ]<- basis_n %*% b[ bid[[ 3 ]] ]
            MMinv_sub_n   <- MMinv[ unlist(bid[ c(1,3) ]), unlist(bid[ c(1,3) ])] ## Added below #########
            basis_sel_n   <- evSqrt2[which(sf_id5==j)][[2]]#[!is.null(basis_n)]
            basis2_n		  <- t( t( basis_n ) * unlist( basis_sel_n ))## added
            x_basis_n	    <- as.matrix( cbind( 1, basis2_n ) )###added
            bse_vc_n0[ ,j ]<- sqrt( colSums( t( x_basis_n ) * ( MMinv_sub_n %*% t( x_basis_n ))))### Added above ##

          } else {
            basis_n        <- NULL
            bse_vc_n0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          }
          j_n      <- j_n + 1
        } else if(1 %in% eff_id[sf_id5==j] ){#eff_id_sub
          basis_n  <- NULL
          bse_vc_n0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )
          j_n      <- j_n + 1
        }

        if( 1 %in% eff_id[sf_id5==j] ){#eff_id_sub
          b_vc0[ ,j ]  <- b[ bid0[ 1 ] ] + b_vc_s0[ ,j ] + b_vc_n0[ ,j ]
          b_vc_s0[ ,j ]<- b[ bid0[ 1 ] ] + b_vc_s0[ ,j ] ### added
          b_vc_n0[ ,j ]<- b[ bid0[ 1 ] ] + b_vc_n0[ ,j ] ### added

          MMinv_sub	   <- MMinv[ unlist(bid), unlist(bid) ]

          if( !is.null(basis_s)|!is.null(basis_n) ){
            basis2       <- cbind( basis_s, basis_n )
            basis_sel     <-evSqrt2[which(sf_id5==j)][c(!is.null(basis_s), !is.null(basis_n))]# fixed
            basis2		   <- t( t( basis2 ) * unlist( basis_sel ) )
            x_basis		   <- as.matrix( cbind( 1, basis2 ) )
            bse_vc0[ ,j ]<- sqrt( colSums( t( x_basis ) * ( MMinv_sub %*% t( x_basis ))))
          } else {
            bse_vc0[ ,j ]<- c( sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ) )#rep(sqrt( MMinv[ bid0[ 1 ], bid0[ 1 ] ] ),n_sub)
          }
        }
      }
    }
    cbind(b_vc0, bse_vc0, b_vc_s0, b_vc_n0, bse_vc_s0, bse_vc_n0, c_vc0, cse_vc0)
  }
  stopCluster(ncores)

  bvc_all   <-as.matrix(bvc_all)
  b_vc      <-as.matrix(bvc_all[,1:nsv])
  bse_vc    <-as.matrix(bvc_all[,(  nsv+1):(2*nsv)])
  b_vc_s0   <-as.matrix(bvc_all[,(2*nsv+1):(3*nsv)])
  b_vc_n0   <-as.matrix(bvc_all[,(3*nsv+1):(4*nsv)])
  bse_vc_s0 <-as.matrix(bvc_all[,(4*nsv+1):(5*nsv)])
  bse_vc_n0 <-as.matrix(bvc_all[,(5*nsv+1):(6*nsv)])
  if( nnxf !=0 ){
    c_vc      <-as.matrix(bvc_all[,(6*nsv+1):(6*nsv + nnxf)])
    cse_vc    <-as.matrix(bvc_all[,(6*nsv + nnxf + 1):(6*nsv + 2*nnxf)])
  } else {
    c_vc    <- NULL
    cse_vc  <- NULL
  }

  pred    <- b_vc[,1]
  if( nsv > 1  ) pred  <- pred + rowSums( X1 * b_vc[, -1 ] )
  if( nnxf !=0 ){
    pred    <- pred + rowSums( Xconst_nvc * c_vc )
  } else {
    c_vc    <- NULL
    if( nxf > 0){
      pred    <- pred + as.matrix( Xconst ) %*% c(b[ 2:( nxf + 1 ) ])
    }
  }

  resid		<- y - pred
  SSE		  <- sum( resid ^ 2 )
  SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
  sig		  <- SSE / ( n - nx )
  bse_vc  <-bse_vc * sqrt( sig )      #### added

  moran_vc  <- rep( 0, nsv )
  j0		    <- 1
  for(j in 1:nsv){
    bid0		<- ifelse( j == 1, 1, nxf + j )
    if( null_dum[ j ] == 0 ){
      bid0_vc	<- ( nx + nev * ( j0 - 1 ) + 1 ):( nx + nev * j0 )
      moran_vc[ j ] <- sum(b[ bid0_vc ]^2*ev)/(ev[1]*sum(b[ bid0_vc ]^2))
      j0		<- j0 + 1
    } else {
      moran_vc[ j ] <- NA
    }
  }

  parV		<- par2[   1:nsv  ]
  parR		<- par2[ (nsv + 1):( 2 * nsv ) ]
  if( nnsv !=0 ){
    parNsv    <- par2[ (2 * nsv + ng + nnxf + 1 ):( 2 * nsv + ng + nnxf + nnsv ) ]
  } else {
    parNsv    <- NULL
  }
  if( nnxf !=0 ){
    parNxf    <- par2[ (2 * nsv +ng + 1 ):( 2 * nsv + ng + nnxf) ]
  } else {
    parNxf    <- NULL
  }

  nsv2		<- sum( parV != 0 )
  nnxf2   <- sum( parNxf != 0 )
  nnsv2   <- sum( parNsv != 0 )
  np		  <- nx + nsv2 * 2 + ng + nnxf2 + nnsv2 + 1#
  AIC		  <- -2 * loglik + np * 2
  BIC		  <- -2 * loglik + np * log( n )
  r2_0		<- 1 - SSE / SSY
  r2		  <- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1 )
  b_cov		<- sig * MMinv
  bse		  <- sqrt( diag( b_cov ) )

  if( is.null( x ) ){
    #b_vc[,1]<- b_vc[,1]# + Bias
    mean_bvc<- mean( b_vc[,1] )
    b_vc[,1]<- b_vc[,1] - mean_bvc

    b_c		  <- b[ 1:( nxf + 1 ) ]
    b_c[1]  <- mean_bvc
    r       <- b[( nx + 1 ):( nx + nev) ]

    bse_c	  <- bse[ 1:( nxf + 1 ) ]
    bz_c	  <- b_c / bse_c
    bp_c	  <- pnorm(abs(bz_c),lower.tail=FALSE) * 2
    b_par		<- data.frame( Estimate = b_c, SE = bse_c, t_value = bz_c, p_value = bp_c )
    rownames( b_par )<- c("(Intercept)", xfname)

  } else {
    if( nxf != 0 ) {
      b_c	  <- b  [ 2:( nxf + 1 ) ]
      bse_c	<- bse[ 2:( nxf + 1 ) ]
      bz_c	<- as.matrix( b_c / bse_c )
      bp_c  <- pnorm(abs(bz_c),lower.tail=FALSE) * 2
      b_par		<- data.frame( Estimate = b_c, SE = bse_c, z_value = bz_c, p_value = bp_c )
      rownames( b_par )<- xfname
    } else {
      b_par		<- NULL
    }
    r       <- NULL
  }

  bz_vc		<- b_vc / bse_vc
  bp_vc   <- pnorm(abs(bz_vc),lower.tail=FALSE) * 2
  b_vc		<- data.frame( b_vc )
  bse_vc	<- data.frame( bse_vc )
  bz_vc		<- data.frame( bz_vc )
  bp_vc		<- data.frame( bp_vc )
  names( b_vc )	<- c( "(Intercept)", xname )
  names( bse_vc )	<- c( "(Intercept)", xname )
  names( bz_vc )	<- c( "(Intercept)", xname )
  names( bp_vc )	<- c( "(Intercept)", xname )

  if(is.null(b_vc_s0) ==FALSE ){
    bz_vc_s0<- as.matrix( b_vc_s0 / bse_vc_s0 )
    bp_vc_s0<- pnorm(abs( bz_vc_s0 ),lower.tail=FALSE) * 2

    b_vc_s0		      <- data.frame( b_vc_s0 )
    bse_vc_s0       <- data.frame( bse_vc_s0 )
    bz_vc_s0		    <- data.frame( bz_vc_s0 )
    bp_vc_s0		    <- data.frame( bp_vc_s0 )
    names( b_vc_s0 )	<- c( "(Intercept)", xname )
    names( bse_vc_s0 )<- c( "(Intercept)", xname )
    names( bz_vc_s0 )	<- c( "(Intercept)", xname )
    names( bp_vc_s0 )	<- c( "(Intercept)", xname )

    B_vc_s <-list(b_vc = b_vc_s0, bse_vc = bse_vc_s0, bz_vc = bz_vc_s0, bp_vc = bp_vc_s0)
  }
  if(is.null(b_vc_n0) ==FALSE ){
    bz_vc_n0<- as.matrix( b_vc_n0 / bse_vc_n0 )
    bp_vc_n0<- pnorm(abs( bz_vc_n0 ),lower.tail=FALSE) * 2

    b_vc_n0		      <- data.frame( b_vc_n0 )
    bse_vc_n0       <- data.frame( bse_vc_n0 )
    bz_vc_n0		    <- data.frame( bz_vc_n0 )
    bp_vc_n0		    <- data.frame( bp_vc_n0 )
    names( b_vc_n0 )	<- c( "(Intercept)", xname )
    names( bse_vc_n0 )<- c( "(Intercept)", xname )
    names( bz_vc_n0 )	<- c( "(Intercept)", xname )
    names( bp_vc_n0 )	<- c( "(Intercept)", xname )

    B_vc_n <-list(b_vc = b_vc_n0, bse_vc = bse_vc_n0, bz_vc = bz_vc_n0, bp_vc = bp_vc_n0)
  }

  if( nnxf != 0 ) {
    cse_vc  <- cse_vc * sqrt( sig )      #### added
    cz_vc		<- as.matrix( c_vc / cse_vc )
    cp_vc   <- pnorm(abs(cz_vc),lower.tail=FALSE)*2
    c_vc		<- data.frame( c_vc )
    cse_vc	<- data.frame( cse_vc )
    cz_vc		<- data.frame( cz_vc )
    cp_vc		<- data.frame( cp_vc )
    names( c_vc )	  <- xfname
    names( cse_vc )	<- xfname
    names( cz_vc )	<- xfname
    names( cp_vc )	<- xfname
  } else {
    c_vc    <- NULL
    cse_vc  <- NULL
    cz_vc   <- NULL
    cp_vc   <- NULL
  }

  parV		  <- parV * sqrt( sig )
  sf_par		<- data.frame( rbind( parV, moran_vc ) )
  names( sf_par )	<- c( "(Intercept)", xname )
  rownames( sf_par )<- c( "random_SE", "Moran.I/max(Moran.I)" )

  e_stat		<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
  rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", lik_nam, "AIC", "BIC")
  vc		<- data.frame(ifelse( sf_par[1,] ==0, 0, 1) )
  names( vc )	<- names( sf_par )
  rownames( vc )	<- "Spatial"


  if( is.null( parNsv ) ==FALSE ){
    vc  <- rbind( vc, c(0, ifelse( parNsv ==0, 0, 1 )))
    rownames( vc )[2]	<- "Non-spatial"

    parN		<- t( parNsv * sqrt( sig ) )
    bn_par	<- data.frame( parN )
    rownames( bn_par )	<- "random_SE"
    names( bn_par )	<- xxname
    s_n    <- bn_par

    totv   <-c(t(sf_par[1,]+c(0,s_n)))#as.vector(sf_par[1,]+c(0,s_n))
    s_rat  <-c(c(t(sf_par[1,]))/totv)#as.vector(sf_par[1,])/totv
    n_rat  <-c(c(0,t(s_n))/totv)#as.vector(c(0,s_n))/totv
    s_rat[totv==0]<-NA;n_rat[totv==0]<-NA
    sn_rat <-data.frame(rbind(s_rat,n_rat))
    rownames(sn_rat)<-c("Share (Spatial)","Share (Non-spatial)")
    names(sn_rat)[1]<-"(Intercept)"
    s       <- list(sf_par, s_n, sn_rat)
  } else {
    s       <- sf_par
  }

  s_nxconst<- NULL
  if( is.null( parNxf ) ==FALSE ){
    vc_n      <- as.data.frame( t(ifelse( parNxf ==0, 0, 1 )))
    names( vc_n )	  <- xxfname
    rownames( vc_n )<- "Non-spatial"
    vc        <- list( vc, vc_n )
    parNf		  <- t( parNxf * sqrt( sig ) )
    bnf_par		<- data.frame( parNf )
    rownames( bnf_par )	<- "random_SE"
    names( bnf_par )	  <- xxfname
    s_nxconst<- bnf_par
  }

  if( sum( n_omit ) >5 ) {
    if( e_stat[2,1]> -0.2 ){
      message( "Note: The model is nearly singular. Consider simplifying the model (dummy variable in x is a typical cause)" )
    } else {
      message( "Note: Singular fit. Simplify the model (dummy variable in x is a typical cause)")
    }
  } else if( e_stat[2,1]< -0.2 ){
    message("Note: Singular fit. Simplify the model (dummy variable in x is a typical cause)")
  }

  #evSqrt2 # evSqrts = evSqrts, evSqrts_n = evSqrts_n,### b_s = bb, b_covs = bb_cov,
  other		<- list( res_int =res_int, r = r, sf_alpha = parR, x_id = x_id, nxf = nxf, xf_id = xf_id,# df = df,
                  model = "resf_vc", nvc_x=nvc_x, nvc_xconst=nvc_xconst, nvc_num = nvc_num, method=method,
                  x = x, xconst = xconst, coords = coords )#, Bias=bias

  result    <- list( b_vc = b_vc, bse_vc = bse_vc, z_vc = bz_vc, p_vc = bp_vc,
                     B_vc_s = B_vc_s, B_vc_n = B_vc_n, c = b_par, c_vc = c_vc, cse_vc = cse_vc,
                     cz_vc = cz_vc, cp_vc = cp_vc, s = s, s_c=s_nxconst, vc = vc, e = e_stat,
                     pred = pred, resid = resid, other = other, call = match.call() )
  class(result)<-"besf_vc"
  return( result )
}

print.besf_vc <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  if( is( x$s )[1] == "list" ){
    cat("\n----Spatially and non-spatially varying coefficients on x (summary)----\n")
  } else {
    cat("\n----Spatially varying coefficients on x (summary)----\n")
  }
  cat("\nCoefficient estimates:\n")
  print( summary( x$b_vc ) )
  cat("\nStatistical significance:\n")
  p01<-apply(x$p_vc,2,function(x) sum(x<0.01))
  p05<-apply(x$p_vc,2,function(x) sum(x<0.05)) - p01
  p10<-apply(x$p_vc,2,function(x) sum(x<0.10)) - p01 - p05
  p90<-length(x$p_vc[,1]) - p01 - p05 - p10
  pv <-data.frame(rbind( p90, p10, p05, p01))
  names(pv)[1]  <- "Intercept"
  row.names(pv) <- c("Not significant", "Significant (10% level)",
                     "Significant ( 5% level)","Significant ( 1% level)")
  print(pv)
  if( is.null(x$s_c) & !is.null(x$c) ){
    cat("\n----Constant coefficients on xconst----------------------------\n")
    print( x$c )

  } else if( !is.null(x$c_vc ) ){
    cat("\n----Non-spatially varying coefficients on xconst (summary)----\n")
    cat("\nCoefficient estimates:\n")
    print( summary( x$c_vc ) )
    cat("\nStatistical significance:\n")
    cp01<-apply(x$cp_vc,2,function(x) sum(x<0.01))
    cp05<-apply(x$cp_vc,2,function(x) sum(x<0.05)) - cp01
    cp10<-apply(x$cp_vc,2,function(x) sum(x<0.10)) - cp01 - cp05
    cp90<-length(x$cp_vc[,1]) - cp01 - cp05 - cp10
    cpv <-data.frame(rbind( cp90, cp10, cp05, cp01))
    row.names(cpv) <- c("Not significant", "Significant (10% level)",
                        "Significant ( 5% level)","Significant ( 1% level)")
    print(cpv)
  }
  cat("\n----Variance parameters----------------------------------\n")
  cat("\nSpatial effects (coefficients on x):\n")
  if( is( x$s )[1] != "list" ){
    print( x$s )
  } else {
    print(x$s[[1]])
    cat("\nNon-spatial effects (coefficients on x):\n")
    print(x$s[[2]])
  }
  if( !is.null(x$s_c) ){
    cat("\nNon-spatial effects (coefficients on xconst):\n")
    print(x$s_c)
  }
  #if( !is.null(x$s_g) ){
  #  cat("\nGroup variation:\n")
  #  print(x$s_g)
  #}
  cat("\n----Error statistics-------------------------------------\n")
  print(x$e)
  if( x$other$method=="reml"){
    cat('\nNote: AIC and BIC are based on the restricted/marginal likelihood.')
    cat('\n      Use method="ml" for comparison of models with different fixed effects (x and xconst)\n')
  }
  invisible(x)
}
