

lsem		<-function( y, x, weig, method = "reml" ){


  lik_llslm	<- function( par, ev, evMax, yy, x, E, Xc, Ec, cy,
                         XX, Xy, EX, Ey,EE, n, nx, ne, emet ){
    ev_inv	<- 1/(evMax - par[1]* ev)

    V_E_diag	<- ev_inv*par[2]

    M00		<- n
    M01		<- Xc
    M02		<- t( t( Ec ) * V_E_diag )
    M11		<- XX
    M12		<- t( EX * V_E_diag )
    M22_0	<- t( t( EE ) * V_E_diag )
    M22		<- t( M22_0 ) * V_E_diag

    M0  <-M <- as.matrix( rbind( cbind( M00, t(M01), t(M02) ),
                                 cbind( M01,  M11 , M12 ),
                                 cbind( M02, t(M12), M22 ) ) )
    diag(M)[-(1:nx)] <- diag(M)[-(1:nx)] + 1
    test    <-try( Minv	<- solve( M, tol = 1e-30 ) )

    m0		<- cy
    m1		<- Xy
    m2		<- t( t( Ey ) * V_E_diag )
    m		<- c( m0, m1, m2 )
    if( class( test )[ 1 ] == "try-error" ){
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

  if( is.null( weig$other$wdum ) ){
    stop( " weig must be defined using the weigen function" )
  }

  n	<- length( y )
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

  ev		<- weig$ev/max(weig$ev)
  evMax	<- 1
  E		<- as.matrix(weig$sf)
  nx		<- dim( X )[ 2 ]
  ne		<- length( ev )
  yy		<- sum( y ^ 2 )
  XX		<- crossprod(  X[,-1] )
  Xy		<- crossprod(  X[,-1], y )
  EX		<- crossprod( weig$sf, X[,-1] )
  Ey		<- crossprod( weig$sf, y )
  EE		<- crossprod( weig$sf )
  Xc		<- colSums(X[,-1])
  Ec		<- colSums(E)
  cy		<- sum(y)
  res	<- optim( fn = lik_llslm, c( 0, 1 ), ev=ev,
                evMax = evMax, yy=yy,XX=XX,Xy=Xy,E=E, x=x,
                Xc = Xc, Ec = Ec, cy = cy,
                EX=EX,Ey=Ey,EE=EE,n=n, nx=nx, ne=ne,
                emet = method, method="L-BFGS-B",
                lower=c(0, 0.001), upper = c(0.995, 5))

  par		<- res$par
  loglik	<- ( -1 / 2 ) * res$value

  ev_inv	<- 1/(evMax - par[1]* ev)

  V_E_diag	<- ev_inv*par[2]

  M00		<- n
  M01		<- Xc
  M02		<- t( t( Ec ) * V_E_diag )
  M11		<- XX
  M12		<- t( EX * V_E_diag )
  M22_0	<- t( t( EE ) * V_E_diag )
  M22		<- t( M22_0 ) * V_E_diag
  M0  <-M <- as.matrix( rbind( cbind( M00, t(M01), t(M02) ),
                               cbind( M01,  M11 , M12 ),
                               cbind( M02, t(M12), M22 ) ) )
  diag(M)[-(1:nx)] <- diag(M)[-(1:nx)] + 1
  Minv	<- solve( M, tol = 1e-30 )

  m0		<- cy
  m1		<- Xy
  m2		<- t( t( Ey ) * V_E_diag )
  m		<- c( m0, m1, m2 )

  b		<- Minv %*% m
  b[ -( 1:nx ) ] <- b[ -( 1:nx ) ] * (ev_inv*par[2])

  XXX		<-as.matrix(cbind(X, weig$sf))
  pred	<-XXX%*%b
  resid	<- y - pred
  SSE		<- sum( resid ^ 2 )
  SSY		<- sum( ( y - mean( y ) ) ^ 2 )
  sig		<- SSE / ( n - nx )
  bse		<- sqrt( sig ) * sqrt( diag( Minv ) )

  np		<- nx + 3
  AIC		<- -2 * loglik + np * 2
  BIC		<- -2 * loglik + np * log( n )
  r2_0	<- 1 - SSE / SSY
  r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1)

  bt		<- b[ 1:nx ] / bse[ 1:nx ]

  Minv_X	<- Minv %*% t( XXX )
  df		<- sum(t(XXX)*Minv_X)
  bp		<- 2 - 2 * pt( abs( bt ), df = n - df )

  b_par	<- data.frame( Estimate = b[ 1:nx ], SE = bse[ 1:nx ], t_value = bt, p_value = bp )
  rownames( b_par ) <- xname

  r_par		<- data.frame( b[ -( 1:nx ) ] )
  names( r_par )	<- "Estimate"
  rownames( r_par )	<- paste( "r", 1:ne, sep = "" )

  par[ 2 ]	<- par[ 2 ] * sqrt( sig )
  sp_par	<- data.frame( par = par )
  rownames( sp_par )  <- c( "sp_lambda", "sp_SE" )
  names( sp_par )     <- "Estimates"
  e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
  if( method == "reml" ){
    rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "rlogLik", "AIC", "BIC" )
  } else if( method == "ml" ){
    rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
  }

  return( list( b = b_par, s = sp_par, e = e_stat, r = r_par, pred = pred, resid = resid ) )
}
