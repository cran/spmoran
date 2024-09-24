esf		<-function( y, x = NULL, vif = NULL, meig, fn = "r2" ){

  if( !is.null(meig$other$coords_z) & meig$other$interact){
    message( paste( "Note: interact=TRUE in meig is ignored because interaction modeling is not implemented" ) )
  }

  ObjEval	<-function( sfmod, fn ){
    if( fn == "r2" ){
      Obj	<- summary( sfmod )$adj.r.squared
    } else if( fn == "aic" ){
      Obj	<- -AIC( sfmod )
    } else if( fn == "bic" ){
      Obj	<- -BIC( sfmod )
    }
    return( Obj )
  }

  n		<- length( y )
  dum		<- 0
  if( is.null( x ) ){
    sfmod	<- lm( y ~ 1 )
    X	<- NULL
    x_id	<- NULL
    nx	<- 1
  } else {
    X00	<- as.matrix( x )
    if( is.numeric( X00 ) == F ){
      mode( X00 ) <- "numeric"
    }
    x_id	<- apply( X00, 2, sd ) != 0
    if( sum( x_id ) == 0 ) {
      sfmod	<- lm( y ~ 1 )
      X	  <- NULL
      x_id<- NULL
      nx	<- 1
    } else {
      X	<- as.matrix( X00[ , x_id ] )
      sfmod	<- lm( y ~ X )
      nx	<- dim( X )[ 2 ] + 1
      vif0	<- max( diag( solve( cor( X ) ) ) )
      if( is.null( vif ) == F ){
        if( vif0 >= vif ){
          message( "Note:" )
          message( paste( "  max(VIF) of the explanatory variables exceeds", vif ) )
          message( "  No eigenvectors are selected" )
          dd	<- data.frame( y, X )
          dum	<- 1
        }
      }
    }
  }

  Xinit	<- X
  sf      	<- meig$sf
  ev      	<- meig$ev
  sf_id     <- rep("s",length(ev))
  if(!is.null(meig$ev_z)){
    sf      <-cbind(sf, meig$sf_z[[1]])
    ev      <-c(ev, meig$ev_z[[1]])
    sf_id   <-c(sf_id, rep("z1", length(meig$ev_z[[1]])))

    if(length(meig$ev_z)==2){
      sf    <-cbind(sf, meig$sf_z[[2]])
      ev    <-c(ev, meig$ev_z[[2]])
      sf_id <-c(sf_id, rep("z2", length(meig$ev_z[[2]])))
    }
  }
  sf        <- as.matrix( sf )

  ne		    <- length( ev )
  sf_list	  <- 1:ne
  Sf_sel_l  <- NULL
  if( dum == 0 ){
    dd	<- data.frame( cbind( y, Xinit ) )
    if( fn != "all" ){
      Obj	<- ObjEval( sfmod = sfmod, fn = fn )
      Obj_old	<- -Inf
      obj	<- Obj
      i	<- 1
      while( ( Obj_old < Obj ) & ( i < ne ) ){
        Obj_old	<- Obj
        Obj_list<- NULL
        for( ii in 1:dim( sf )[ 2 ] ){
          d	<- data.frame( cbind( y, X, sf[ , ii ] ) )
          names(d)[1]<-"y"
          sfmod0	<- lm( y ~ ., d )
          Obj	<- ObjEval( sfmod = sfmod0, fn = fn )
          Obj_list<- c( Obj_list, Obj )
        }
        Obj	<- max( Obj_list )
        if( Obj_old < Obj ){
          obj	<- c( obj, Obj )
          Obj_ind0<- ( 1:length( Obj_list ) ) * ( Obj_list == Obj )
          Obj_ind	<- Obj_ind0[ Obj_ind0 != 0 ][ 1 ]
          sf_sel_l<- sf_list[    Obj_ind ]
          sf_sel	<- sf     [ ,  Obj_ind ]
          sf_list	<- sf_list[   -Obj_ind ]
          sf	<- sf     [ , -Obj_ind ]
          X0	<- cbind( X, sf_sel )
          if( is.null(vif) == F ){
            VIF	<- diag( solve( cor( X0 ) ) )
            if( max( VIF ) < vif ){
              X	<- data.frame( X0 )
              names( X )[ dim( X )[ 2 ] ] <- paste( "sf", sf_sel_l, sep = "" )
              dd	<- data.frame( y, X )
              sfmod	<- lm( y ~ ., dd )
            }
          } else {
            X	<- data.frame( X0 )
            names( X )[ dim( X )[ 2 ] ] <- paste( "sf", sf_sel_l, sep = "" )
            dd	<- data.frame( y, X )
            sfmod	<- lm( y ~ ., dd )
          }
          Sf_sel_l	<- c( Sf_sel_l, sf_sel_l )
          i		<- i + 1
        }# else {
        #	if( i == 1 ){
        #		dd	<- data.frame( y, Xinit )
        #	}
        #}
        if( i %% 50 == 0 ) print( i )
      }
      if( i == ne ){
        Obj_old			<- Obj
        X0			<- cbind( Xinit, sf )#meig$sf
        if( is.null( vif ) ){
          vif_id		<- 0
        } else {
          VIF		<- diag( solve( cor( X0 ) ) )
          vif_id		<- ifelse( max( VIF ) < vif, 0, 1 )
        }

        if( vif_id == 0 ){
          dd0		<- data.frame( y, X0 )
          names( dd0 )[ 1 ]  <- "y"
          names( dd0 )[ ( nx + 1 ):( nx + ne) ]<- paste( "sf", 1:ne, sep = "" )
          sfmod0		<- lm( y ~ ., dd0 )
          Obj		<- ObjEval( sfmod = sfmod0, fn = fn )
          if( Obj_old < Obj ){
            X	<- X0
            dd	<- dd0
            sfmod	<- sfmod0
            Sf_sel_l<- 1:ne
          }
        }
      }
    } else {
      dd	<- data.frame( cbind( y, X ) )
      for( ii in sf_list ){
        X0	<- cbind( X, sf[, ii ] )
        if( is.null( vif ) == FALSE ){
          VIF	<- diag( solve( cor( X0 ) ) )
          if( max( VIF ) < vif ){
            X	<- data.frame( X0 )
            names( X )[ dim( X )[ 2 ] ] <- paste( "sf", ii, sep = "" )
            dd	<- data.frame( y, X )
            Sf_sel_l	<- c( Sf_sel_l, ii )
          }
        } else {
          X	<- data.frame( X0 )
          names( X )[ dim( X )[ 2 ] ] <- paste( "sf", ii, sep = "" )
          dd	<- data.frame( y, X )
          Sf_sel_l	<- c( Sf_sel_l, ii )
        }
      }
      sf    <- X[,-(1:(nx-1))]
      sfmod	<- lm( y ~ ., dd )
    }
  }
  pred	<- predict( sfmod )
  resid	<- residuals ( sfmod )
  df		<- dim( X )[ 2 ] + 1
  sig		<- sum( resid ^ 2 ) / ( n - df )
  if( nx == 1 ){
    b_par	<- summary( sfmod )$coefficients[   1:nx , 1:4 ]
    names( b_par ) <- c( "(Intercept)", "SE","t_value", "p_value" )
  } else {
    b_par	<- data.frame( summary( sfmod )$coefficients[   1:nx , 1:4 ] )
    names( b_par ) <- c( "Estimate", "SE","t_value", "p_value" )
  }

  if( dim( dd )[ 2 ] == nx ){
    r	<- NULL
    nr	<- 0
  } else {
    if( dim( dd )[ 2 ] == ( nx + 1 ) ){
      r	<- data.frame( t( summary( sfmod )$coefficients[ -(1:nx), 1:4 ] ) )
      row.names( r ) <- names( dd )[ dim( dd )[ 2 ] ]
      nr	<- 1
    } else {
      r	<- data.frame( summary( sfmod )$coefficients[ -(1:nx), 1:4 ] )
      nr	<- dim( r )[ 1 ]
    }
    names( r )	<- c( "Estimate", "SE","t_value", "p_value" )
  }

  if( is.null( r ) ){
    SF	<- NULL
  } else {
    if( nx == 1 ){
      SF	<- as.matrix( X ) %*% c(r[,1])
    } else {
      SF	<- as.matrix( X[ , -( 1:( nx - 1 ) ) ] ) %*% c(r[,1])
    }
  }
  vif 	<- data.frame(diag( solve( cor( X ) ) ) )
  names(vif)	<- "VIF"

  len        <- 1:min(length(r[,1]), length(ev[Sf_sel_l]))
  sf_moran   <-sum(r[,1][len]^2*ev[Sf_sel_l][len])/(ev[1]*sum(r[,1][len]^2))
  sf_par	<- data.frame( par = c( sd(SF), sf_moran )  )
  names( sf_par )   <- "Estimate"
  rownames( sf_par )<- c( "SD", "Moran.I/max(Moran.I)" )

  r2		  <- summary( sfmod )$adj.r.squared
  loglik	<- logLik( sfmod )
  AIC		  <- AIC( sfmod )
  BIC		  <- BIC( sfmod )
  e_stat	<- data.frame( stat = c( sqrt( sig ), r2, loglik, AIC, BIC ) )
  rownames( e_stat ) <- c( "resid_SE", "adjR2", "logLik", "AIC", "BIC")
  if( dum ==0) {
    message( paste( "  ", nr, "/", ne, " eigenvectors are selected", sep = "" ) )
  }

  other	<- list( x_id = x_id, sf_id = Sf_sel_l, model = "esf", coords = meig$other$coords )
  result<- list( b = b_par, s = sf_par, r = r, vif = vif, e = e_stat, sf = SF,
                 pred = pred, resid = resid, other = other, call = match.call() )
  class( result ) <- "esf"
  return( result )
}

print.esf <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Coefficients------------------------------\n")
  print(x$b)
  cat("\n----Variance parameter (residuals)---------------\n")
  print(x$s)
  cat("\n----Error statistics--------------------------\n")
  print(x$e)
  invisible(x)
}

