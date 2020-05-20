resf_qr <- function( y, x = NULL, meig, tau = NULL, boot = TRUE, iter = 200, cl=NULL ){

  resf_00  	<- function( y, x = NULL, xgroup = NULL, meig, method = "reml" ){

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
      if( class( test )[1] == "try-error" ){
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


    resf_boot  	<-function( tau_sel, q_sel, y, X, M, meig, mod, iter ){

      re_esfb	<-function( y, X, XSF0, M, meig ){

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
          if( class( test )[1] == "try-error" ){
            loglik  <- Inf
          } else {
            b	<- Minv %*% m
            sse	<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
            dd	<- sse + sum( b[ -( 1:nx ) ] ^ 2 )
            if( ( emet == "reml" ) | ( emet == "pls" ) ){
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

        ev	<- meig$ev
        n	<- length( y )
        nx	<- dim( X )[ 2 ]
        ne	<- length( ev )
        yy	<- sum( y^2 )
        m	<- c( crossprod( X, y ), crossprod( meig$sf, y ) )
        res 	<- optim( fn = lik_resf, c( 1, 1 ), ev = ev, M = M, m = m,
                       yy = yy, n = n, nx = nx, ne = ne, emet = "pls" )
        par	<- res$par ^ 2
        loglik 	<- ( -1 / 2 ) * res$value
        evv     <- ev ^ par[ 1 ] * ( sum( ev ) / sum( ev ^ par[ 1 ] ) )
        evSqrt  <- par[ 2 ] * sqrt( evv )
        sf2     <- t( t( meig$sf ) * evSqrt )
        X2     	<- as.matrix( cbind( X, sf2 ) )
        XX2  	<- crossprod( X2 )
        diag( XX2 )[ -( 1:nx ) ] <- diag( XX2 )[ -( 1:nx ) ] + 1
        XXinv2	<- solve( XX2 )
        b	<- XXinv2 %*% t( X2 ) %*% y
        b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * evSqrt
        sig		<- sum( ( y - XSF0 %*% b ) ^ 2 ) / ( n - nx )
        par[ 2 ]	<- par[ 2 ] * sqrt( sig )
        par       <- par[ 2:1 ]
        sf_moran  <- sum( b[ -( 1:nx ) ] ^ 2 * ev )/( ev[ 1 ] * sum( b[ -( 1:nx ) ] ^ 2 ) )
        return(list( b = b[ 1:nx ], par = par, sf_moran = sf_moran ) )
      }

      b	<- mod$b[  , 1 ]
    	s	<- c( mod$s[ 1, ], mod$other$sf_alpha )
    	sd	<- mod$e[ 1, 1 ]
    	ev	<- meig$ev
    	n	<- length( y )
    	ne	<- length( ev )
    	nx	<- dim( X )[ 2 ]

    	evv	<- ev ^ s[ 2 ] * ( sum( ev ) / sum( ev ^ s[ 2 ] ) )
    	evSqrt	<- s[ 1 ] * sqrt( evv )
    	sf2	<- t( t( meig$sf ) * evSqrt )
    	XSF0	<- as.matrix( cbind( X, meig$sf ) )
    	f0	<- density( y )
    	fq0	<- approx( f0$x, f0$y, q_sel )$y

    	res<- foreach( i = 1:iter, .combine = "rbind" ) %dopar% {
    		usim	<- rnorm( n , sd = sd )
    		rsim	<- rnorm( ne, sd = 1  )
    		y_b	<- X %*% b + sf2 %*% rsim + usim
   		f_b	<- density( sample(y, replace=TRUE) )
    		fq_b	<- approx( f_b$x, f_b$y, q_sel )$y
    		RIF_b	<- fq0 / fq_b * ( y_b - q_sel) + q_sel
    		sfres_sim	<- re_esfb( y = RIF_b, X = X, XSF0 = XSF0, M = M, meig = meig )
    		c( sfres_sim$b, c( sfres_sim$par[ 1 ], sfres_sim$sf_moran) )
    	}
    	B	<- as.matrix( res[ , 1:nx ] )
    	S	<- as.matrix( res[ ,( nx + 1): (nx + 2) ] )
    	return( list( B = B, S = S ) )
    }

    indic <- function( y, q_sel, tau_sel ){
    	if( tau_sel == 1 ){
    		res <-ifelse( y < q_sel, 1, 0 )
    	} else {
    		res <-ifelse( y <= q_sel, 1, 0 )
    	}
    	return(res)
    }

    if(is.null( tau ) ) tau <- 1:9 / 10
    q	  <- quantile( y, probs = tau )
    f	  <- density( y )
    fq	<- approx( f$x, f$y, q )$y
    n   <- length( y )
    ne  <- length( meig$ev )

    if( is.null( x ) ){
    	X 	<- as.matrix( rep( 1, n ) )
    	xname 	<- "(Intercept)"
    } else {
    	X00	<- as.matrix( x )
    	if( is.numeric( X00 ) == F ){
    		mode( X00 ) <- "numeric"
    	}
    	xind	<- apply( X00, 2, sd ) != 0
    	X0	<- X00[ , xind ]
    	X	<- as.matrix( cbind( 1, X0 ) )
        if( sum( xind ) == 0 ){
    		xname <- "(Intercept)"
        } else {
    		xname <- c( "(Intercept)", names( as.data.frame( X0 ) ) )
    	}
    }

    if( boot == TRUE ){
    	XX	<- crossprod( X )
    	EX	<- crossprod( meig$sf, X )
    	if( meig$other$fast == 0 ){
    		EE	<- diag( length( meig$ev ) )
    	} else if( meig$other$fast == 1 ){
    		EE	<- crossprod( meig$sf )
    	}
    	M	<- as.matrix( rbind( cbind( XX, t( EX ) ), cbind( EX, EE ) ) )

    	if(is.null(cl)) {
    	  cl <- makeCluster(detectCores(),setup_strategy = "sequential")
    	} else {
    	  cl <- makeCluster(cl,setup_strategy = "sequential")
    	}
    	registerDoParallel(cl)
    }

    SFb  	<- NULL
    SFr  	<- NULL
    SFs		<- NULL
    SFe		<- NULL
    SFb_boot	<-list( NULL )
    SFs_boot	<-list( NULL )
    SFb_boot2 <-list( NULL )
    SFs_boot2 <-list( NULL )
    probs	    <- c( 0.025, 0.975 )
    for( j in 1:length( tau ) ){
      RIF	<- q[ j ] + ( tau[ j ] - indic( y = y, q_sel = q[ j ], tau_sel = tau[ j ] ) ) / fq[ j ]
    	mod	<- resf_00( y = RIF, meig = meig, x = X )
    	SFb	<- cbind( SFb	, mod$b[ , 1 ] )
    	SFr	<- cbind( SFr	, mod$r[ , 1 ] )
    	SFs	<- cbind( SFs	, mod$s[ 1:2, ]	)
    	SFe	<- cbind( SFe	, mod$e[ 1:2, ]	)
    	if( boot == TRUE ){
    		mod_b	<-resf_boot( tau_sel = tau[ j ], q_sel = q[ j ], y = y,
    				     X = X, M = M, meig = meig, mod = mod, iter = iter )
    		SFb_boot_0  <- data.frame( t( mod_b$B ) )
    		SFs_boot_0  <- data.frame( t( mod_b$S ) )
    		SFb_boot2_00 <- t( apply( mod_b$B, 2, function( x ) quantile( x, probs = probs ) ) )
    		SFs_boot2_00 <- t( apply( mod_b$S, 2, function( x ) quantile( x, probs = probs ) ) )
    		SFb_p  <- 1 - abs( apply( mod_b$B, 2, function( x ) sum( x > 0 ) ) - iter / 2 ) / ( iter / 2 )
    		SFb_boot2_0  <- data.frame( mod$b[ , 1 ] , SFb_boot2_00, SFb_p )
    		SFs_boot2_0  <- data.frame( mod$s[ 1:2, ], SFs_boot2_00 )

    		rownames( SFb_boot_0 ) <- rownames( SFb_boot2_0 ) <- xname
    		rownames( SFs_boot_0 ) <- rownames( SFs_boot2_0 ) <- c( "spcomp_SE", "spcomp_Moran.I/max(Moran.I)" )

    		names( SFb_boot_0 ) <- paste( "iter", 1:iter, sep = "" )
    		names( SFs_boot_0 ) <- paste( "iter", 1:iter, sep = "" )
    		SFb_boot[[ j ]]  <- SFb_boot_0
    		SFs_boot[[ j ]]  <- SFs_boot_0

    		names( SFb_boot2_0 ) <- c( "Estimates", "CI_lower", "CI_upper", "p_value" )
    		names( SFs_boot2_0 ) <- c( "Estimates", "CI_lower", "CI_upper" )
    		SFb_boot2[[ j ]]  <- SFb_boot2_0
    		SFs_boot2[[ j ]]  <- SFs_boot2_0

    		gc(); gc()
    	}
    	print( paste( "------- Complete: tau=", tau[ j ], " -------", sep = "" ) )
    }

    SFb		    <- data.frame( SFb )
    SFr		    <- data.frame( SFr )
    SFs		    <- data.frame( SFs )
    SFe		    <- data.frame( SFe )
    rownames( SFb ) <- xname
    rownames( SFr ) <- paste( "r", 1:ne, sep = "" )
    rownames( SFs ) <- c( "spcomp_SE", "spcomp_Moran.I/max(Moran.I)" )
    rownames( SFe ) <- c( "resid_SE", "quasi_adjR2(cond)" )

    tau_name	        <- paste( "tau=", tau, sep= "" )
    names( SFb )      <- names( SFr ) <-names( SFs ) <- names( SFe ) <- tau_name
    if( boot == FALSE ){
      res   <- list( b = SFb, r = SFr, s = SFs, e = SFe, tau = tau, call = match.call() )
    } else {
      tau_name2         <- paste( "--------------- tau=", tau, " ---------------", sep = "")
      names( SFb_boot ) <- names( SFb_boot2 ) <- tau_name2
      names( SFs_boot ) <- names( SFs_boot2 ) <- tau_name2
      res   <- list( b = SFb, r = SFr, s = SFs, e = SFe, B0 = SFb_boot,
                     S0 = SFs_boot, B = SFb_boot2, S = SFs_boot2, tau = tau, call = match.call() )
      stopCluster( cl )
    }

    class(res)<-"resf_qr"
return( res )
}

print.resf_qr <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Coefficients------------------------------\n")
  print(x$b)
  cat("\n----Spatial effects (residuals)---------------\n")
  print(x$s)
  cat("\n----Error statistics--------------------------\n")
  print(x$e)
  invisible(x)
}
