besf  	<- function( y, x = NULL, coords, method = "reml",
                    covmodel="exp",enum = 200, bsize = 3000, cl=NULL ){

    lik_resf	<- function( par0, ev, M, m, yy, n, nx, ne, emet ){
    	par	<- par0 ^ 2
    	evv	<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    	evSqrt	<- par[ 2 ] * sqrt( evv )
    	Mw	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    	M[  -( 1:nx ), -( 1:nx ) ]	<- Mw + diag( ne )
    	M[  -( 1:nx ),    1:nx   ]	<- M[ -( 1:nx ), 1:nx ] * evSqrt
    	M[     1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ), 1:nx ] )
    	M0				<- M
    	M0[ -( 1:nx ), -( 1:nx ) ]	<- Mw
    	m[  -( 1:nx) ]			<- m[ -( 1:nx ) ] * evSqrt

    	test    <-try( Minv	<- solve( M, tol = 1e-25 ) )
    	if( class( test ) == "try-error" ){
    		loglik  <- Inf
    	} else {
    		b	<- Minv %*% m
    		sse	<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
    		dd	<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
    		if( emet == "reml" ){
    			term1	<- determinant( M )$modulus
    			term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
    		} else if( emet == "ml" ){
    			term1	<- determinant( as.matrix( M[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
    			term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
    		}
    		loglik	<- term1 + term2
    	}
    	return( loglik )
    }

    n		<- length( y )
    if( is.null( x ) ){
    	X	<- as.matrix( rep( 1, n ) )
    	xname	<- "(Intercept)"
    	x_id	<- NULL
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == F ){
    		mode( X00 ) <- "numeric"
    	}
    	x_id	<- apply( X00, 2, sd ) != 0
    	if( sum( x_id ) == 0 ){
    		X	<- as.matrix( rep( 1, n ) )
    		xname	<- "(Intercept)"
    		x_id	<- NULL
    	} else {
    		X0	<- X00[ , x_id ]
    		X	<- as.matrix( cbind( 1, X0 ) )
    		xname	<- c( "(Intercept)", names( as.data.frame( X0 ) ) )
    	}
    }

    nx		<- dim( X )[ 2 ]
    yy    <- sum( y ^ 2 )
    XX		<- crossprod(  X )
    Xy		<- crossprod(  X, y )
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
    ne      <- length(ev)

    nblock  <-max(1,round((n*( nx + ne )/bsize^2)))
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
      X2   <- cbind( X[ids_sub,  ], sf)
      sfsf0<-cbind(crossprod( X2, y[ids_sub] ),crossprod( X2 ))
      sfsf1<-cbind(rep(0,length(sfsf0[,1])),sfsf0)

      clsf<-colSums(sf)
      sfsf1[1:length(clsf),1]<-clsf
      sfsf1
    }

    sfmean <- sfsf[1:sum(ev_full > 1e-08),1]/n
    EX     <- as.matrix( sfsf[, -(1:2)][,1:nx] )
    sf_id  <- c(1:nx,rep(1,sum(ev_full > 1e-08)))
    sfmean2<- c(rep(0,nx), sfmean)
    EX2    <- t(t(EX[,sf_id])*sfmean2)
    XX2    <- XX[sf_id,sf_id]*outer(sfmean2,sfmean2)
    M      <- sfsf[, -(1:2)] - EX2 - t(EX2) + XX2
    m      <- c(sfsf[,2]) - sfmean2*Xy[sf_id]
    yy     <- sum( y ^ 2 )

    res		<- optim( fn = lik_resf, c( 1, 1 ), ev = ev, M = M, m = m, yy = yy,
    			   n = n, nx = nx, ne = ne, emet = method )
    par		<- res$par ^ 2
    loglik<- ( -1 / 2 ) * res$value
    evv		<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    evSqrt<- par[ 2 ] * sqrt( evv )

    M[ -( 1:nx ), -( 1:nx ) ]	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    M[ -( 1:nx ),    1:nx   ]	<-    M[ -( 1:nx ),    1:nx   ] * evSqrt
    M[    1:nx  , -( 1:nx ) ]	<- t( M[ -( 1:nx ),    1:nx   ] )
    diag( M [ -( 1:nx ), -( 1:nx ) ] ) <- diag( M[ -( 1:nx ), -( 1:nx ) ] ) + 1
    MMinv		<- solve( M, tol = 1e-30 )
    m[ -( 1:nx ) ]	<- m[ -( 1:nx ) ] * evSqrt
    b		<- MMinv %*% m
    b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * evSqrt

    nxx   <- nx
    b_s		<- list(NULL)
    b_covs		<- list(NULL)
    #evSqrts	<- list(NULL)
    i     <- NULL
    SF    <- foreach(i = 1:(length(ids)-1), .packages="fields",.combine="+") %dopar%  {

          if(i==1){
            ids_sub<-ids[i]:ids[i+1]
          } else {
            ids_sub<-(ids[i]+1):ids[i+1]
          }

          sf	 <- t( exp( -rdist( coordk, coords[ids_sub,] ) / h ) - Cmean ) %*% sf_ap0[ , ev_full > 1e-08 ]
          sf  <- t( t( sf ) - sfmean )
          sf %*% c( b[-(1:nx)] )
    }
    stopCluster(cl)

    pred	<- X %*% b[1:nx] + SF
    resid	<- y - pred
    SSE		<- sum( resid ^ 2 )
    SSY		<- sum( ( y - mean( y ) ) ^ 2 )
    sig		<- SSE / ( n - nx )
    bse		<- sqrt( sig ) * sqrt( diag( MMinv ) )

    np		<- nx + 1 + 3
    AIC		<- -2 * loglik + np * 2
    BIC		<- -2 * loglik + np * log( n )
    r2_0	<- 1 - SSE / SSY
    r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1)

    bz		<- b[ 1:nx ] / bse[ 1:nx ]
    bp    <- pnorm(abs(bz),lower.tail=FALSE)*2

    b_par	<- data.frame( Estimate = b[ 1:nx ], SE = bse[ 1:nx ], z_value = bz, p_value = bp )
    rownames( b_par ) <- xname

    r_par		<- data.frame( b[ -( 1:nx ) ] )
    names( r_par )	<- "Estimate"
    rownames( r_par )	<- paste( "r", 1:ne, sep = "" )

    sf_moran<-sum(b[ -( 1:nx ) ]^2*ev)/(ev[1]*sum(b[ -( 1:nx ) ]^2))
    par[ 2 ]<- par[ 2 ] * sqrt( sig )
    sf_par	<- data.frame( par = c( par[ 2 ], sf_moran )  )
    names( sf_par )   <- "Estimate"
    rownames( sf_par )<- c( "spcomp_SE", "spcomp_Moran.I/max(Moran.I)" )

    e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    if( method == "reml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "rlogLik", "AIC", "BIC" )
    } else if( method == "ml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
    }
    other	<- list( sf_alpha= par[ 1 ], x_id = x_id, model = "resf", par0 = res$par )

    return( list( b = b_par, s = sf_par, e = e_stat,
    		  r = r_par, sf = SF, pred = pred, resid = resid, other = other ) )
}
