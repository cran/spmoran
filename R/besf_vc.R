
besf_vc		<- function( y, x = NULL, xconst = NULL, coords, method = "reml", penalty = "bic", allsvc=FALSE, maxiter = 30, sizelimit = 2000,
                        covmodel="exp",enum = 200, bsize = 3000, cl=NULL ){

  Mdet_f0	  	<- function( M, M0, id, par0_sel, emet ){
    if( emet == "ml" ){
      M	<- M [ id != 0, id != 0 ]
      M0	<- M0[ id != 0, id != 0 ]
      id	<- id[ id != 0 ]
    }

    if(sum( id != par0_sel ) == 0 ){
      term2	<- NULL
      term3_0	<- NULL
    } else {
      M0_sub		<- M0[ id != par0_sel, id != par0_sel ]
      term2		<- determinant( M0_sub )$modulus

      Msub_00		<- M[ id == par0_sel, id == par0_sel ]
      Msub_01		<- M[ id == par0_sel, id != par0_sel ]
      term3_0		<- Msub_00 - Msub_01 %*% solve( M0_sub, tol = 1e-30 ) %*% t( Msub_01 )
    }
    return(list(term2 = term2, term3_0 = term3_0))
  }

  Mdet_f	  	<- function( evSqrt, id, term2, term3_0, par0_sel ){
    if( is.null( term2 )  ){
      Mdet		<- sum( log( evSqrt ) ) * 2
    } else {
      term1		<- sum( log( evSqrt ) ) * 2
      diag( term3_0 ) <- diag( term3_0 ) + 1/evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
      term3	  	<- determinant( term3_0 )$modulus
      Mdet	  	<- term1 + term2 + term3
    }
    return(Mdet)
  }

  lik_resf_vc		<- function( par0, par0_est, par0_id, par0_sel, ev, M, M0inv, M0inv_01, M0inv_00,
                            m, yy, b_01, b_02, n, nx, nsv, emet, term2, term3_0, null_dum2, id ){
    par		<- par0 ^ 2
    par_est		<- par0_est ^ 2
    par[ par0_id == par0_sel ]  <- par_est
    evSqrt	<- NULL
    for( i in ( 1:nsv )[ null_dum2[ 1:nsv ] == 0 ] ){
      evv  	<- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
      evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
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

  lik_resf_vc0	<- function( par0, ev, M, m, yy, n, nx, nsv, emet, null_dum4 ){
    par		<- par0 ^ 2
    evSqrt		<- NULL
    for( i in ( 1:nsv )[ null_dum4[ 1:nsv ] == 0 ] ){
      evv  	<- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
      evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
    }

    M[ -( 1:nx ), -( 1:nx ) ]	<- t(M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    M[ -( 1:nx ),    1:nx   ]	<-   M[ -( 1:nx ),    1:nx   ] * evSqrt
    M[    1:nx  , -( 1:nx ) ]	<- t(M[ -( 1:nx ),    1:nx   ] )
    M0		<- M
    diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1

    m[-(1:nx)]		<- m[ -( 1:nx ) ] * evSqrt
    test			<-try(Minv	<- solve( M, tol = 1e-30 ))
    if(class(test)[ 1 ] == "try-error"){
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


  if( is.null( x ) ){
  	result	<- besf( y = y, x = xconst, method = "reml", coords=coords,
                   covmodel = covmodel, enum = enum, bsize = bsize, cl = cl )
  	b_par	  <- result$b
  	sf_par	<- result$s
  	e_stat	<- result$e
  	b_vc	  <- NULL
  	bse_vc	<- NULL
  	bz_vc	  <- NULL
  	bp_vc	  <- NULL
  	pred	  <- result$pred
  	resid	  <- result$resid
  	vc	    <- NULL
  	r	      <- result$r
  	other	  <- result$other
  	other$nxf<-other$nx

  } else {
  if( method == "reml" ){
    lik_nam	<- "rlogLik"
  } else if( method == "ml" ){
    lik_nam	<- "logLik"
  }
  n     	<- length( y )
  X1	<- x
  if( is.null( X1 ) == FALSE ){
    X1	<- as.matrix( X1 )
    if( is.numeric( X1 ) == FALSE ){
      mode( X1 ) <- "numeric"
    }
    x_id	<- apply( X1, 2, sd ) != 0
    nx0	<- sum( x_id )
    if( nx0 == 0 ){
      X1	<- NULL
      xname	<- NULL
      x_id	<- NULL
    } else {
      X1	<- as.matrix( X1[ , x_id ] )
      xname	<- names( as.data.frame( X1 ) )
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
      xf_id	<- NULL
    } else {
      Xconst	<- as.matrix( Xconst[ , xf_id ] )
      xfname	<- names( as.data.frame( Xconst ) )
    }
  } else {
    xfname	<- NULL
    xf_id	<- NULL
    nxf	<- 0
  }

  if( ( nx0 > 0 ) && ( nxf > 0 ) ){
    for( dd  in 1:nx0  ){
      for( ddf in 1:nxf ){
        if( sum( X1[ , dd ] != Xconst[ , ddf ] ) == 0 ){
          stop( " x and xconst cannot have the same column" )
        }
      }
    }
  }

  nsv	<- ifelse( is.null( X1 ), 1, dim( X1 )[ 2 ] + 1 )
  coordk	<- kmeans( coords, centers = enum + 1 )$centers
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

  nev0	<- min( round( n / nsv ) - 4, length( ev ))
  if( is.null( sizelimit ) == FALSE ){
    if( nev0 * nsv + nxf > sizelimit ){
      nev1	<- round(( sizelimit - nxf ) / nsv )
      if( nev0 - nev1 > 0 ){
        if( nev1 < 10){
          nev1	<- 10
          message( paste( " Note: sizelimit was too small. It was replaced with ", nsv * 10 + nxf, sep="" ) )
        }
        nev0	<- nev1
      }
    }
  }


  id	<- c( rep( 0, nsv ), rep( 0, nxf ), rep( 1, length( ev ) ) )
  if( nsv >= 2 ){
    for( i in 1:( nsv - 1 ) ){
      id	<- c( id, rep( i+1, length( ev )))
    }
  }

  X0	<- as.matrix( cbind( 1, Xconst, X1 ) )
  nx	<- dim( X0 )[ 2 ]
  Mo 	<- crossprod( X0 )
  mo	<- crossprod( X0, y )
  parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )

  nblock  <-max(1,round((n*(nxf + length(ev)*nsv)/bsize^2)))
  ids   <- unique(c(round(seq(1,n,len=nblock+1)),n))

  if(is.null(cl)) {
    cl <- makeCluster(detectCores())
  } else {
    cl <- makeCluster(cl)
  }
  registerDoParallel(cl)

  im     <- NULL
  sfsf   <- foreach(im = 1:(length(ids)-1), .packages=c("fields"),.combine="+") %dopar%  {
    if(im==1){
      ids_sub<-ids[im]:ids[im+1]
    } else {
      ids_sub<-(ids[im]+1):ids[im+1]
    }
    sf	 <- t( exp( -rdist( coordk, coords[ids_sub,] ) / h ) - Cmean ) %*% sf_ap0[ , ev_full > 1e-08 ]
    X2   <- cbind(X0[ids_sub,  ], sf)
    if( nsv >= 2 ){
      for( im2 in 1:( nsv - 1 ) ){
        X2  <- cbind( X2, X1[ids_sub, im2 ] * sf )
      }
    }
    sfsf0<-cbind(crossprod( X2, y[ids_sub] ),crossprod( X2 ))
    sfsf1<-cbind(rep(0,length(sfsf0[,1])),sfsf0)

    clsf<-colSums(sf)
    sfsf1[1:length(clsf),1]<-clsf
    sfsf1
  }

  sfmean<-sfsf[1:sum(ev_full > 1e-08),1]/n
  m0     <-c(sfsf[,2])
  EX     <-sfsf[, -(1:2)][,1:nsv]
  sf_id  <-c(1,rep(1,nxf),2:nsv,rep(1:nsv,each=sum(ev_full > 1e-08)))
  sf_ide <-rep(1:sum(ev_full > 1e-08), nsv)
  sfmean2<-c(rep(0,nxf + nsv), sfmean[sf_ide])
  EX2    <-t(t(EX[,sf_id])*sfmean2)
  XX2    <-Mo[sf_id,sf_id]*outer(sfmean2,sfmean2)
  M      <-sfsf[, -(1:2)] - EX2 - t(EX2) + XX2
  m      <-m0 - sfmean2*mo[sf_id]
  yy     <- sum( y ^ 2 )
  if( penalty == "aic" ){
    pen	<- 2
  } else if( penalty == "bic" ){
    pen	<- log( n )
  }
  par0	<- rep( 1, 2 * nsv )
  Par	<- par0
  par0_est<- c( 1, 1 )
  par0_id	<- rep( 1:nsv, 2 )
  null_dum<- rep( 0, nsv )

  obj	<- sd( y )
  LL	<- NULL
  res_old	<- Inf
  n_omit	<- 0
  iter	<- 1
  while( ( obj > sd( y ) * 0.0001 ) & ( iter <= maxiter ) ){

    print( paste( "-------  Iteration ", iter, "  -------", sep = "" ) )
    LL0	<- LL
    n_omit	<- 0
    for( par0_sel in 1:nsv ){
      evSqrt	<- NULL
      par0_sq	<- par0 ^ 2
      for( i in 1:nsv ){
        evv  	<- ev ^ par0_sq[ nsv + i ] * sum( ev ) / sum( ev ^ par0_sq[ nsv + i ] )
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

      null_dum2<- null_dum
      null_dum2[ par0_sel ] <- 0
      null_dum3<- c( rep( 0, nx ), null_dum2[ id ]) == 0
      MM	<- M [ null_dum3, null_dum3 ]
      MM0	<- M0[ null_dum3, null_dum3 ]
      mm	<- m [ null_dum3 ]
      idd	<- id[ null_dum3 ]
      id_omit1<- c( diff( id ), 1)

      id_omit2<- which( id_omit1 != 0 )
      id_omit1[ id_omit2[ 1 ]] <- 0

      er_dum	<- TRUE
      n_omit0	<- 0
      while( er_dum == TRUE ){
        try1	<- try( M0inv	<- solve( MM0,  tol = 1e-30 ), silent = TRUE )
        try2	<- try( Mdet0	<- Mdet_f0( M = MM, M0 = MM0, id = idd,
                                       par0_sel = par0_sel, emet = method ), silent = TRUE)
        er_dum	<- ( class(try1)[1] =="try-error" ) | ( class(try2)[1] =="try-error" )
        if( er_dum == TRUE ){
          M	<- M [ id_omit1 == 0, id_omit1 == 0 ]
          M0	<- M0[ id_omit1 == 0, id_omit1 == 0 ]
          m	<- m [ id_omit1 == 0 ]
          id	<- id[ id_omit1 == 0 ]
          X	<- X [, id_omit1 == 0 ]
          null_dum3<- null_dum3[ id_omit1 == 0 ]

          MM	<- M [ null_dum3, null_dum3 ]
          MM0	<- M0[ null_dum3, null_dum3 ]
          mm	<- m [ null_dum3 ]
          idd	<- id[ null_dum3 ]

          id_omit1<- c( diff( id ), 1)

          id_omit2<- which( id_omit1 != 0 )
          id_omit1[ id_omit2[ 1 ] ] 	<- 0
          n_omit0	<- n_omit0 + 1
        }
      }

      if( n_omit0 > 0 ){
        message( paste( "Note: ", n_omit0, " eigenvectors are omitted to stablize the estimates", sep = "" ) )
      }
      n_omit	<- c( n_omit, n_omit0 )

      #ev	<- ev     [   1:sum( id == 1 ) ]
      #meig$sf	<- meig$sf[ , 1:sum( id == 1 ) ]
      term2	<- Mdet0$term2
      term3_0	<- Mdet0$term3_0

      if( min( par0[ par0_id == par0_sel ] ) >= 1e-5 ){
        par00	<- par0[ par0_id == par0_sel ]
      } else {
        par00_id<- max( which( apply(Par[,par0_id == par0_sel], 1, min) != 0 ) )
        par00	<- Par[ par00_id, par0_id == par0_sel ]
      }

      if( n_omit0 > 0 ){
        res_old <- lik_resf_vc0( par0, ev = ev, M = MM, m = mm, yy = yy,
                                 n = n, nx = nx, nsv = nsv, emet = method, null_dum4 = null_dum2 )
      }

      M0inv_01<- M0inv[ 		, idd == par0_sel ]
      M0inv_00<- M0inv[ idd ==par0_sel, idd == par0_sel ]
      b_01	<- M0inv %*% mm
      b_02	<- t( M0inv_01 ) %*% mm

      if( par0_sel == 1){
        llim    <- c( parVmax_sq / 1000, 1e-03)
        ulim    <- c( parVmax_sq, 4 )
        omethod <- "L-BFGS-B"
        res    	<- optim( fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                          M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                          term2 = term2, term3_0 = term3_0, m = mm, yy = yy,
                          n = n, nx = nx, nsv = nsv, emet = method, id = idd,
                          par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2,
                          lower = llim, upper = ulim, method = omethod )
      } else if( par0_sel <= nsv ){
        res    	<- optim( fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                          M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                          term2 = term2, term3_0 = term3_0, m = mm, yy = yy,
                          n = n, nx = nx, nsv = nsv, emet = method, id = idd,
                          par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2 )
      } else {
        res    	<- optimize( f = lik_resf_vc, lower = 0.001, upper = 100,par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                             M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                             term2 = term2, term3_0 = term3_0, m = mm, yy = yy,
                             n = n, nx = nx, nsv = nsv, emet = method, id = idd,
                             par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2 )
        res$value	<- c( res$objective )
        res$par	<- c( res$minimum )
      }

      if( ( iter > 3 ) & ( res$value > res_old ) ){
        loglike	<- ( - 1 / 2 ) * res_old
      } else {
        if( ( par0_sel != 1 ) & ( par0_sel <= nsv ) ){
          MM_ref 	<- MM [ idd != par0_sel, idd != par0_sel ]
          mm_ref	<- mm [ idd != par0_sel ]
          null_dum4<- null_dum2
          null_dum4[ par0_sel ] <- 1
          res_ref <- lik_resf_vc0( par0, ev = ev, M = MM_ref, m = mm_ref, yy = yy,
                                   n = n, nx = nx, nsv = nsv, emet = method, null_dum4 = null_dum4 )
        } else {
          res_ref	<- Inf
        }

        if( (res_ref < res$value + pen)&(allsvc==FALSE)){
          par0[ par0_id == par0_sel ] 	<- 0
          null_dum[ par0_sel ]		<- 1
          loglik 	<- ( -1 / 2 ) * c( res_ref )
          res_old	<- res_ref
        } else {
          par0[ par0_id == par0_sel ]	<- res$par
          null_dum[ par0_sel ]		<- 0
          res_old	<- res$value
          loglik 	<- ( -1 / 2 ) * res_old
        }
      }

      LL0	<- c( LL0, loglik )
      print( paste( par0_sel, "/", nsv, sep = "" ) )

    }

    if( iter > 1 ){
      if( sum( n_omit ) == 0 ){
        obj 	<- abs( loglik - loglik0 )
      } else {
        obj	<- Inf
      }
    }

    Par	<- rbind( Par, par0 )
    loglik0	<- loglik
    LL	<- c( LL, loglik )
    print( paste( lik_nam, ": ", round( loglik, 3 ), sep = "" ) )
    iter	<- iter + 1
  }

  par2   	<- par0 ^ 2
  evSqrt	<- NULL
  evSqrt2	<- list(NULL)
  for( i in 1:nsv ){
    evv  	<- ev ^ par2[ nsv + i ] * sum( ev ) / sum( ev ^ par2[ nsv + i ] )
    evSqrt	<- c( evSqrt, par2[ i ] * sqrt( evv ) )
    evSqrt2[[ i ]]	<- par2[ i ] * sqrt( evv )
  }

  MM	<- M [ null_dum3, null_dum3 ]
  MM0	<- M0[ null_dum3, null_dum3 ]
  mm	<- m [ null_dum3 ]
  idd	<- id[ null_dum3 ]

  M[ -( 1:nx ), -( 1:nx ) ]	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
  M[ -( 1:nx ),    1:nx   ]	<-    M[ -( 1:nx ),    1:nx   ] * evSqrt
  M[    1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ),    1:nx   ] )
  diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1

  MM		<- M [ null_dum3, null_dum3 ]
  MMinv		<- solve( MM, tol = 1e-30 )

  mm		<- m [ null_dum3 ]
  eevSqrt		<- evSqrt[ null_dum3[ -( 1:nx )] ]
  mm[ -( 1:nx ) ]	<- mm[ -( 1:nx ) ] * eevSqrt
  b		<- MMinv %*% mm
  b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * eevSqrt

  nxx <- nx
  nev		<- length( ev )
  b_s		<- list(NULL)
  b_covs		<- list(NULL)
  #evSqrts		<- list(NULL)
  bvc_all   <- foreach(i = 1:(length(ids)-1), .packages="fields",.combine="rbind") %dopar%  {
    j0		<- 1
    for(j in 1:nsv){
    bid0		<- ifelse( j == 1, 1, nxf + j )
    if( null_dum[ j ] == 0 ){
      bid0_vc	<- ( nx + nev * ( j0 - 1 ) + 1 ):( nx + nev * j0 )
      bid		  <- c( bid0, bid0_vc )

      if(i==1){
        ids_sub<-ids[i]:ids[i+1]
      } else {
        ids_sub<-(ids[i]+1):ids[i+1]
      }

      if(j==1){
        bvc_all0<- matrix(0, ncol=nsv, nrow= length(ids_sub))
        bse_all0<- matrix(0, ncol=nsv, nrow= length(ids_sub))
      }
      sf	 <- t( exp( -rdist( coordk, coords[ids_sub,] ) / h ) - Cmean ) %*% sf_ap0[ , ev_full > 1e-08 ]
      sf  <- t( t( sf ) - sfmean )

      b_vc	<- b[ bid0 ] + sf %*% b[ bid0_vc ]
      MMinv_sub	<- MMinv[ bid, bid ]
      sf2		<- t( t( sf ) * evSqrt2[[ j ]] )
      x_sf		<- as.matrix( cbind( 1, sf2 ) )
      bse_vc	<- sqrt( colSums( t( x_sf ) * ( MMinv_sub %*% t( x_sf ) ) ) )

      #b_s[[ j ]]	<- c(b[ bid0 ], b[ bid0_vc ] )
      #b_covs[[ j ]]	<- b_cov_sub
      #evSqrts[[ j ]]	<- evSqrt2[[ j ]]
      j0		<- j0 + 1
    } else {
      n_sub <-length(min(ids_sub):max(ids_sub))
      b_vc	  <- rep(b[ bid0 ],n_sub)
      bse_vc	<- rep(sqrt( MMinv[ bid0, bid0 ] ),n_sub)
      #b_s[[ j ]]	  <- rep(b[ bid0 ], n_sub)
      #b_covs[[ j ]]	<- rep(b_cov[ bid0, bid0 ], n_sub)
     # evSqrts[[ j ]]<- 0
    }
    bvc_all0[,j]<-b_vc
    bse_all0[,j]<-bse_vc
    }
    cbind(bvc_all0, bse_all0)
  }
  stopCluster(cl)

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

  b_vc    <- bvc_all[,(1:nsv)]
  if( nxf !=0  ){
    pred    <- b_vc[,1] + rowSums( X1 * b_vc[,-1] ) + as.matrix(Xconst) %*% c(b[ 2:( nxf + 1 ) ])
  } else {
    pred    <- b_vc[,1] + rowSums( X1 * b_vc[,-1] )
  }
  resid		<- y - pred
  SSE		<- sum( resid ^ 2 )
  SSY		<- sum( ( y - mean( y ) ) ^ 2 )
  sig		<- SSE / ( n - nxx )
  bse_vc  <-bvc_all[,-(1:nsv)]*sqrt( sig )
  #b_cov		<- sig * MMinv

  parV		<- par2[   1:nsv  ]
  parR		<- par2[ (nsv + 1):( 2 * nsv ) ]
  nsv2		<- sum( parV != 0 )

  np		<- nxx + nsv2 * 2 + 1
  AIC		<- -2 * loglik + np * 2
  BIC		<- -2 * loglik + np * log( n )
  r2_0		<- 1 - SSE / SSY
  r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1 )

  if( nxf != 0 ) {
    bse		<- sqrt( diag( sig * MMinv ) )

    b		  <- b  [ 2:( nxf + 1 ) ]
    bse		<- bse[ 2:( nxf + 1 ) ]
    bz		<- b / bse
    bp    <- pnorm(bz,lower.tail=FALSE)*2
    b_par		<- data.frame( Estimate = b, SE = bse, z_value = bz, p_value = bp )
    rownames( b_par )<- xfname
  } else {
    b_par		<- NULL
  }

  bz_vc		<- b_vc / bse_vc
  bp_vc   <- pnorm(abs(bz_vc),lower.tail=FALSE)*2
  b_vc		<- data.frame( b_vc )
  bse_vc		<- data.frame( bse_vc )
  bz_vc		<- data.frame( bz_vc )
  bp_vc		<- data.frame( bp_vc )
  names( b_vc )	<- c( "(Intercept)", xname )
  names( bse_vc )	<- c( "(Intercept)", xname )
  names( bz_vc )	<- c( "(Intercept)", xname )
  names( bp_vc )	<- c( "(Intercept)", xname )

  parV		  <- parV * sqrt( sig )
  sf_par		<- data.frame( rbind( parV, moran_vc ) )
  names( sf_par )	<- c( "(Intercept)", xname )
  rownames( sf_par )<- c( "spcomp_SE", "spcomp_Moran.I/max(Moran.I)" )

  re2_vc    <-as.matrix(b_vc)
  e_stat		<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
  rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", lik_nam, "AIC", "BIC")
  vc		<- data.frame(ifelse( sf_par[1,] ==0, 0, 1) )
  names( vc )	<- names( sf_par )
  rownames( vc )	<- "varying coefficients"
  r		<- NULL
  other		<- list( sf_alpha = parR, x_id = x_id, nxf = nxf, xf_id = xf_id, b_s = b_s, b_covs = b_covs,  model = "resf_vc" )
  }
  return( list( b = b_par, s = sf_par, e = e_stat, b_vc = b_vc, bse_vc = bse_vc, z_vc = bz_vc, p_vc = bp_vc,
                pred = pred, resid = resid, vc = vc, r = r, other = other ) )
}






