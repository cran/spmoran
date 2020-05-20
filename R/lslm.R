

lslm	<-function( y, x, weig, method = "reml", boot = FALSE, iter = 200 ){

 lik_llslm	<- function( par, ev, evMax, yy, x, E, Xc, Ec, cy,
    			XX, Xy, EX, Ey,EE, n, nx, ne, emet ){
    ev_inv0	<- 1/(evMax - par[1]* ev)
    ev_inv	<- ev_inv0

    ev_invX0	<- par[1]*ev*ev_inv0
    ev_invX	<- ev_invX0

    V_XE	<- t(ev_invX * EX)
    V_E_diag	<- ev_inv*par[2]

    M11_term2	<- 2 * V_XE %*% EX
    M11_term3_0	<- EE %*% t(V_XE)
    M11_term3	<- V_XE %*% M11_term3_0
    M11		<- XX + M11_term2 + M11_term3

    M12_0	<- t(EX) + t( M11_term3_0 )
    M12		<- t( t( M12_0 ) * V_E_diag )

    M22_0	<- t( t( EE ) * V_E_diag )
    M22		<- t( M22_0 ) * V_E_diag

    M00		<- n
    M01		<- Xc + V_XE %*% Ec
    M02		<-t( t( Ec ) * V_E_diag )


    M0  <-M <- as.matrix( rbind( cbind( M00, t(M01), t(M02) ),
			  cbind( M01,  M11 , M12 ),
    			  cbind( M02, t(M12), M22 ) ) )
    diag(M)[-(1:nx)] <- diag(M)[-(1:nx)] + 1
    test    <-try( Minv	<- solve( M, tol = 1e-30 ) )

    m0		<- cy
    m1		<- Xy + V_XE %*% Ey
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

ev	<- weig$ev/max(weig$ev)
E	<- as.matrix(weig$sf)
nx	<- dim( X )[ 2 ]
ne	<- length( ev )
yy     	<- sum( y ^ 2 )
XX	<- crossprod(  X[,-1] )
Xy	<- crossprod(  X[,-1], y )
EX	<- crossprod( weig$sf, X[,-1] )
Ey	<- crossprod( weig$sf, y )
EE	<- crossprod( weig$sf )
Xc	<- colSums(X[,-1])
Ec	<- colSums(E)
cy	<- sum(y)

evMax	<- 1
res	<- optim( fn = lik_llslm, c( 0, 1 ), ev=ev,
    		evMax = evMax, yy=yy,XX=XX,Xy=Xy,E=E, x=x,
    		Xc = Xc, Ec = Ec, cy = cy,
    		EX=EX,Ey=Ey,EE=EE,n=n, nx=nx, ne=ne,
    		emet = method, method="L-BFGS-B",
    		lower=c(-0.995, 0.001), upper = c(0.995, 10))

    par		<- res$par
    loglik	<- ( -1 / 2 ) * res$value

    ev_inv0	<- 1/(evMax - par[1]* ev)
    ev_inv	<- ev_inv0

    ev_invX0	<- par[1]*ev*ev_inv0
    ev_invX	<- ev_invX0

    V_XE	<- t(ev_invX * EX)
    V_E_diag	<- ev_inv*par[2]

    M11_term2	<- 2 * V_XE %*% EX
    M11_term3_0	<- EE %*% t(V_XE)
    M11_term3	<- V_XE %*% M11_term3_0
    M11		<- XX + M11_term2 + M11_term3

    M12_0	<- t(EX) + t( M11_term3_0 )
    M12		<- t( t( M12_0 ) * V_E_diag )

    M22_0	<- t( t( EE ) * V_E_diag )
    M22		<- t( M22_0 ) * V_E_diag

    M00		<- n
    M01		<- Xc + V_XE %*% Ec
    M02		<-t( t( Ec ) * V_E_diag )


    M0  <-M <- as.matrix( rbind( cbind( M00, t(M01), t(M02) ),
			  cbind( M01,  M11 , M12 ),
    			  cbind( M02, t(M12), M22 ) ) )
    diag(M)[-(1:nx)] <- diag(M)[-(1:nx)] + 1
    Minv	<- solve( M, tol = 1e-30 )

    m0		<- cy
    m1		<- Xy + V_XE %*% Ey
    m2		<- t( t( Ey ) * V_E_diag )
    m		<- c( m0, m1, m2 )

    b		<- Minv %*% m
    b[ -( 1:nx ) ] <- b[ -( 1:nx ) ] * (ev_inv*par[2])

    V_X		<- ev_invX*EX
    x2		<-as.matrix(cbind(1, X[,-1] + E%*%V_X))
    XXX		<-as.matrix(cbind(x2, weig$sf))
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

    par[ 2 ]	<- par[ 2 ]# * sqrt( sig )
    sp_par	<- data.frame( par = par )
    rownames( sp_par ) <- c( "sp_rho", "sp_SE" )
    names( sp_par ) <- c( "Estimates" )
    e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
    if( method == "reml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "rlogLik", "AIC", "BIC" )
    } else if( method == "ml" ){
    	rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", "logLik", "AIC", "BIC" )
    }


    ev_wei	<- par[ 1 ] * ev/( evMax - par[1]* ev )
    VE		<- t( E ) * ev_wei
    DE_vec	<- colSums( t( E ) * VE ) + 1
    IE_vec	<- 1 + E %*% c(ev_wei * Ec) - DE_vec
    DE		<- mean(DE_vec) * b[2:nx]
    IE		<- mean(IE_vec) * b[2:nx]

    if( boot == TRUE ){
    	ev_inv		<- 1/(evMax - par[1]* ev)
    	ev_invE		<- par[2] * ev_inv
    	ev_invX		<- par[1] * ev * ev_inv
    	V_XE		<- t(ev_invX * EX)

    	V_E_diag	<- ev_inv * par[2]
    	xb      <- x2 %*% b[ 1:nx ]

    	sp_par_boot	<- matrix(0, nrow=iter, ncol=2 )
    	b_boot		<- matrix(0, nrow=iter, ncol=nx )
    	DE_boot		<- matrix(0, nrow=iter, ncol=nx - 1 )
    	IE_boot		<- matrix(0, nrow=iter, ncol=nx - 1 )
    	for( it in 1:iter ){
    		u_boot		<- rnorm( n, sd = e_stat[1,1] )
    		r_boot0   <- ev_inv * rnorm( ne )
    		r_boot		<- sd(r_par[,1])/sd(r_boot0) * r_boot0
    		y_boot		<- xb + E %*% r_boot + u_boot

    		yy     		<- sum( y_boot ^ 2 )
    		Xy		<- crossprod(  X[,-1], y_boot )
    		Ey		<- crossprod( weig$sf, y_boot )
    		cy		<- sum(y_boot)

    		res_boot<- optim( fn = lik_llslm, c(0, 1), ev=ev,
    			evMax = evMax, yy=yy,XX=XX,Xy=Xy,E=E, x=x,
    			Xc = Xc, Ec = Ec, cy = cy,
    			EX=EX,Ey=Ey,EE=EE,n=n, nx=nx, ne=ne,
    			emet = method, method="L-BFGS-B",
    			lower=c(-0.995, 0.001), upper = c(0.995, 10))

    		ev_inv	<- 1/(evMax - res_boot$par[1]* ev)
    		ev_invX	<- res_boot$par[1]*ev*ev_inv
    		V_XE	<- t(ev_invX * EX)
    		V_E_diag	<- ev_inv*res_boot$par[2]

    		M11_term2	<- 2 * V_XE %*% EX
    		M11_term3_0	<- EE %*% t(V_XE)
    		M11_term3	<- V_XE %*% M11_term3_0
    		M11		<- XX + M11_term2 + M11_term3

    		M12_0	<- t(EX) + t( M11_term3_0 )
    		M12		<- t( t( M12_0 ) * V_E_diag )

    		M22_0	<- t( t( EE ) * V_E_diag )
    		M22		<- t( M22_0 ) * V_E_diag

    		M00		<- n
    		M01		<- Xc + V_XE %*% Ec
    		M02		<-t( t( Ec ) * V_E_diag )

    		M0  	<-M <- as.matrix( rbind( cbind( M00, t(M01), t(M02) ),
			  cbind( M01,  M11 , M12 ),
    			  cbind( M02, t(M12), M22 ) ) )
    		diag(M)[-(1:nx)] <- diag(M)[-(1:nx)] + 1
    		Minv	<- solve( M, tol = 1e-30 )

    		m0		<- cy
    		m1		<- Xy + V_XE %*% Ey
    		m2		<- t( t( Ey ) * V_E_diag )
    		m		<- c( m0, m1, m2 )
    		b_boot0		<- ( Minv %*% m )[ 1:nx ]
    		b_boot[ it, ]	<- b_boot0

    		sp_par_boot[ it, ]<- res_boot$par
    		ev_wei		<- res_boot$par[ 1 ] * ev/( evMax - res_boot$par[1]* ev )
    		VE		<- t( E ) * ev_wei
    		DE_vec		<- colSums( t( E ) * VE ) + 1
    		IE_vec		<- 1 + E %*% c(ev_wei * Ec) - DE_vec
    		DE_boot[ it, ]<- b_boot0[ -1 ] * mean(DE_vec)
    		IE_boot[ it, ]<- b_boot0[ -1 ] * mean(IE_vec)
    		if( it %% 20 == 0 ){
    			print( paste( "------- Complete:", it, "/", iter, " -------", sep = "" ) )
    			gc();gc()
    		}
    	}

    	probs	    	<- c( 0.025, 0.975 )
    	DE_CI 		<- t( apply( DE_boot, 2, function( x ) quantile( x, probs = probs ) ) )
    	DE_p		<- 1 - abs( apply( DE_boot, 2, function( x ) sum( x > 0 ) ) - iter / 2 ) / ( iter / 2 )
    	DE_res  	<- data.frame( DE , DE_CI, DE_p )
    	names( DE_res ) <- c( "Estimates", "CI_lower", "CI_upper", "p_value" )
    	rownames( DE_res ) <- rownames(b_par)[ -1 ]

    	IE_CI 		<- t( apply( IE_boot, 2, function( x ) quantile( x, probs = probs ) ) )
    	IE_p		<- 1 - abs( apply( IE_boot, 2, function( x ) sum( x > 0 ) ) - iter / 2 ) / ( iter / 2 )
    	IE_res		<- data.frame( IE , IE_CI, IE_p )
    	names( IE_res ) <- c( "Estimates", "CI_lower", "CI_upper", "p_value" )
    	rownames( IE_res ) <- rownames(b_par)[ -1 ]

    	sp_par_boot[,2] <- sp_par_boot[,2] - median(sp_par_boot[,2])+sp_par[2,]
    	sp_par_boot[,2][ sp_par_boot[,2] < 0 ] <- 0
    	s_CI 		<- t( apply( sp_par_boot, 2, function( x ) quantile( x, probs = probs ) ) )
    	sp_par		<- data.frame( sp_par, s_CI )
    	sp_par[ 2, ]    <- sp_par[ 2, ] * sqrt( sig )
    	names( sp_par ) <- c( "Estimates", "CI_lower", "CI_upper" )

    } else {
    	DE_res		<- data.frame( DE )
    	names( DE_res ) <- c( "Estimates" )
    	rownames( DE_res ) 	<- rownames(b_par)[ -1 ]

    	IE_res		<- data.frame( IE )
    	names( IE_res ) <- c( "Estimates" )
    	rownames( IE_res ) 	<- rownames(b_par)[ -1 ]

    	sp_par[ 2, ]    <- sp_par[ 2, ] * sqrt( sig )
    }

    result      <- list( b = b_par, s = sp_par, e = e_stat, de = DE_res, ie = IE_res, r = r_par,
                         pred = pred, resid = resid, call = match.call() )
    class( result ) <- "lslm"
    return( result )
}


print.lslm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Coefficients------------------------------\n")
  print(x$b)
  cat("\n----Spatial effects (lag)---------------------\n")
  print(x$s)
  cat("\n----Effects estimates-------------------------\n")
  if( dim(x$de)[ 2 ] >1 ){
    cat("\nDirect:\n")
    print(x$de)
    cat("\nIndirect:\n")
    print(x$ie)
  } else {
    xx  <- data.frame( x$de, x$ie )
    names(xx)<- c("Direct","Indirect")
    print(xx)
  }
  cat("\n----Error statistics--------------------------\n")
  print(x$e)
}
