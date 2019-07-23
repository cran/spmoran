resf  	<- function( y, x = NULL, xgroup = NULL, meig, method = "reml" ){

  lik_resf	<- function( par0, ev, M, m, yy, n, nx, nxx, ne, xg_id ){
    par	<- par0 ^ 2
    evv	<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
    evSqrt	<- par[ 2 ] * sqrt( evv )
    if( is.null( xg_id ) == FALSE ){
      for( k in 1:max( xg_id )){
        evSqrt<-c( evSqrt, rep( par[ 2 + k ], sum( xg_id == k ) ) )
      }
      neg  <- ne + length( xg_id )
    } else {
      neg   <-ne
    }
    Mw	<- t( M[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
    M[  -( 1:nx ), -( 1:nx ) ]	<- Mw + diag( neg )
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
      term1	<- determinant( M )$modulus
      term2	<- ( n - nxx ) * ( 1 + log( 2 * pi * dd / ( n - nxx ) ) )
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
    if( is.numeric( X00 ) == FALSE ){
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

  if( is.null( xgroup ) == FALSE ){
    xgroup  <- data.frame(xgroup)
    ng	    <- dim( xgroup )[ 2 ]
    xg_id0	<- 1
    xg_idid2<- NULL
    xg_link_id <- list( NULL )
    for( ff in 1:ng ){
      xg00	  <- factor( xgroup[ , ff ] )
      xg0	    <- sparse.model.matrix( ~ 0 + xg00 )
      xg0s	  <- Matrix::colSums( xg0 )
      xg_idid	<- max( which( xg0s == max( xg0s ) ) )
      Xg0	    <- xg0[ , -xg_idid ]
      #names( Xg0 ) <- paste( names( as.data.frame(xgroup) )[ ff ], "_", levels( xg00 )[ -xg_idid ], sep = "" )
      xg_id1	<- rep( xg_id0,  dim( Xg0 )[2] )
      if( ff == 1 ){
        Xg	  <- Xg0
        xg_id	<- xg_id1
      } else {
        Xg	  <- cbind( Xg, Xg0 )
        xg_id	<- c( xg_id, xg_id1 )
      }
      xg_id0	<- xg_id0 + 1
      xg_idid2<-c(xg_idid2, xg_idid)

      xgroup_datt <- data.frame(id=1:length(xgroup[,ff]),xgroup_sub=xgroup[,ff])
      xg00_datt   <- data.frame(id_b_g = 1:length(levels(xg00)), xgroup_sub=levels(xg00))
      xgroup_datt2<- merge(xgroup_datt,xg00_datt,by="xgroup_sub",all.x=TRUE)
      xg_link_id[[ ff ]]  <- xgroup_datt2[order(xgroup_datt2$id),"id_b_g"]
    }

  } else {
    ng	<- 0
    xg_id<-NULL
  }

  ev		<- meig$ev
  nx		<- dim( X )[ 2 ]
  ne		<- length( ev )
  yy     	<- sum( y ^ 2 )
  XX		<- crossprod(  X )
  Xy		<- crossprod(  X, y )
  nxx <- nx

  if( ng == 0 ){
    EX		<- crossprod( meig$sf, X )
    Ey		<- crossprod( meig$sf, y )
    if( meig$other$fast == 0 ){
      EE	<- diag( ne )
    } else if( meig$other$fast == 1 ){
      EE	<- crossprod( meig$sf )
    }
  } else {
    meig$sf <-as.matrix( cbind(meig$sf, Xg) )
    EX		<- crossprod( meig$sf, X )
    Ey		<- crossprod( meig$sf, y )
    EE	  <- crossprod( meig$sf )
  }

  M		<- as.matrix( rbind( cbind( XX, t( EX ) ), cbind( EX, EE ) ) )
  m		<- c( Xy, Ey )

  par0<- c( 1, 1 )
  if( ng != 0 ) par0  <- c( par0, rep( 1, ng ) )
  res		<- optim( fn = lik_resf, par0, ev = ev, M = M, m = m, yy = yy,
                   n = n, nx = nx, nxx  = nxx, ne = ne, xg_id = xg_id )
  par		<- res$par ^ 2
  loglik	<- ( -1 / 2 ) * res$value
  evv		<- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
  evSqrt	<- par[ 2 ] * sqrt( evv )

  if( is.null( xg_id ) == FALSE){
    for( k in 1:length( par[-c(1:2)])){
      evSqrt<-c( evSqrt, rep( par[ 2 + k ], sum( xg_id == k ) ) )
    }
  }

  sf2		<- t( t( meig$sf ) * evSqrt )
  X2		<- as.matrix( cbind( X, sf2 ) )
  XX2		<- crossprod( X2 )
  diag( XX2 )[ -( 1:nx ) ] <- diag( XX2 )[ -( 1:nx ) ] + 1
  XXinv2	<- solve( XX2 )
  XXinv2_X	<- XXinv2 %*% t( X2 )
  b		<- XXinv2_X %*% y
  b[ -( 1:nx ) ] <- b[ -( 1:nx ) ] * evSqrt
  pred	<- as.matrix( cbind( X, meig$sf ) ) %*% b
  resid	<- y - pred
  SSE		<- sum( resid ^ 2 )
  SSY		<- sum( ( y - mean( y ) ) ^ 2 )
  sig		<- SSE / ( n - nxx )
  bse		<- sqrt( sig ) * sqrt( diag( XXinv2 ) )
  SF		<- meig$sf[, 1:ne ] %*% b[ ( nx + 1 ):( nx + ne ) ]

  np		<- nxx + 1 + 2 + ng
  AIC		<- -2 * loglik + np * 2
  BIC		<- -2 * loglik + np * log( n )
  r2_0	<- 1 - SSE / SSY
  r2		<- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1)

  r		  <- b  [ ( nx + 1 ):( nx + ne ) ]
  rse	  <- bse[ ( nx + 1 ):( nx + ne ) ]
  rt	  <- r / rse
  rpar  <- data.frame( Estimate = r, SE = rse, t_value = rt )
  rownames( rpar )	<- paste( "r", 1:ne, sep = "" )

  Bias  <- 0
  if( ng != 0 ){
    bpar_g <- data.frame( Estimate = b[ -(1:(nx+ne))], SE = bse[ -(1:(nx+ne))], t_value = NA )
    #b_g		 <- b  [ -( 1:( nx + ne ) ) ]

    b_g2   <-list(NULL)
    nulldat<-data.frame(Estimate = 0, SE = NA, t_value = NA)
    for(ggid in 1:ng){
      b_g		   <- bpar_g[xg_id==ggid,]
      xg_idid_0<- xg_idid2[ggid]
      if(xg_idid_0 == 1){
        b_g		 <- rbind(nulldat,b_g)
      } else if( xg_idid_0 == dim(b_g)[1] ){
        b_g		 <- rbind(b_g,nulldat)
      } else {
        b_g		 <- rbind(b_g[1:(xg_idid_0-1),],nulldat,b_g[-(1:(xg_idid_0-1)),])
      }
      bias     <- mean(b_g[xg_link_id[[ggid]],1])
      b_g[,1]  <- b_g[,1] - bias
      Bias     <- Bias + bias

      b_g[is.na(b_g[,2])==FALSE,3]<- b_g[is.na(b_g[,2])==FALSE,1]/b_g[is.na(b_g[,2])==FALSE,2]
      xg00	   <- factor( xgroup[ , ggid ] )
      rownames(b_g) <- paste( names( as.data.frame(xgroup) )[ ggid ], "_", levels( xg00 ), sep = "" )
      b_g2[[ggid]] <- b_g
    }

  } else {
    b_g2  <- NULL
  }

  b[1]  <- b[1] + Bias
  bt		<- b[ 1:nx ] / bse[ 1:nx ]
  df		<- sum(t(X2)*XXinv2_X)
  bp		<- 2 - 2 * pt( abs( bt ), df = n - df )
  b_par	<- data.frame( Estimate = b[ 1:nx ], SE = bse[ 1:nx ], t_value = bt, p_value = bp )
  rownames( b_par ) <- xname

  sf_moran<-sum(r^2*ev)/(ev[1]*sum(r^2))
  par[ 2 ]<- par[ 2 ] * sqrt( sig )
  sf_par	<- data.frame( par = c( par[ 2 ], sf_moran )  )
  names( sf_par )   <- "Estimate"
  rownames( sf_par )<- c( "spcomp_SE", "spcomp_Moran.I/max(Moran.I)" )

  if( ng != 0 ){
    parG		<- t( par[ -( 1:2 )] * sqrt( sig ) )
    bg_par	<- data.frame( parG )
    rownames( bg_par )	<- "group_SE"
    names( bg_par )	<- names( as.data.frame( xgroup ) )
  } else {
    bg_par  <-NULL
  }

  e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
  rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
  other	<- list( sf_alpha= par[ 1 ], x_id = x_id, model = "resf", par0 = res$par, nx = nx, df = df, bias=Bias )

  return( list( b = b_par, b_g = b_g2, s = sf_par, s_g = bg_par, e = e_stat,
                r = rpar, sf = SF, pred = pred, resid = resid, other = other ) )
}
