predict0_vc	<- function( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL ){

  {

    af <-function(par,y) par[1]+par[2]*y
    sc <-function(par,y) (y-par[1])/par[2]
    sa <-function(par,y) sinh(par[1]*asinh(y)-par[2])
    bc <-function(par,y, jackup){

      yy     <- y + ifelse(jackup, abs( par[2] ), 0 )
      if(par[1]==0){
        res <- log(yy)
      } else {
        res <- ( yy^par[1] - 1 ) / par[1]
      }
    }
    sal<-function(par,y, noconst=FALSE){
      if(noconst==FALSE){
        par[3]+ par[4]*sinh(par[1]*asinh(y)-par[2])
      } else {
        sinh(par[1]*asinh(y)-par[2])
      }
    }

    d_af<-function(par,y) par[2]
    d_sc <-function(par,y) 1/par[2]
    d_sa<-function(par,y) par[1]*cosh(par[1]*asinh(y)-par[2])/sqrt(1+y^2)
    d_bc<-function(par,y,jackup) (y + ifelse(jackup, abs( par[2] ), 0 ) )^(par[1]-1)#abs(y)^(par-1)
    d_sal<-function(par,y,noconst=FALSE){
      if(noconst==FALSE){
        d_af(par=par[3:4],sa(par=par[1:2],y))*d_sa(par=par[1:2],y)
      } else {
        d_af(par=c(0, 1),sa(par=par[1:2],y))*d_sa(par=par[1:2],y)
      }
    }

    i_af<-function(par,y) (y-par[1])/par[2]
    i_sc <-function(par,y) par[1]+par[2]*y
    i_sa<-function(par,y) sinh(1/par[1]*(asinh(y)+par[2]))
    i_bc<-function(par,y,jackup){
      if(par[1]==0){
        exp(y) - ifelse(jackup, abs( par[2] ), 0 )
      } else {
        (par[1]*y + 1)^(1/par[1]) - ifelse(jackup, abs( par[2] ), 0 )
      }
    }

    i_sal<-function(par,y,noconst=FALSE){
      if(noconst==FALSE){
        i_sa(par=par[1:2],i_af(par=par[3:4],y))
      } else {
        i_sa(par=par[1:2],i_af(par=c(0, 1),y))
      }
    }

    d_sal_k   <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL,y_ms=NULL,z_ms,jackup){
      d       <-1
      if( !is.null( bc_par[1] ) ){
        d     <- d_bc(bc_par,y,jackup=jackup)*d
        y     <-   bc(bc_par,y,jackup=jackup)
      }

      if(k>0){
        d       <- d_sc(par=y_ms,y)*d
        y       <-   sc(par=y_ms,y)

        for(kk in 1:k){
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          d     <-d_sal(par[[kk]],y,noconst=nc_dum)*d
          y     <-  sal(par[[kk]],y,noconst=nc_dum)
        }
      }
      d       <-d_sc(par=z_ms,y)*d
      return(d)
    }
    i_sal_k   <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL,y_ms=NULL,z_ms, jackup){

      y       <-i_sc(par=z_ms,y)
      if(k>0){
        for(kk in k:1) {
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          y <- i_sal(par[[kk]],y,noconst=nc_dum)
        }
      }

      if( !is.null( y_ms ) ){
        y     <-i_sc(par=y_ms,y)
      }

      if( !is.null(bc_par[1]) ){
        y <- i_bc(par=bc_par, y, jackup=jackup)
      }

      return(y)
    }
    sal_k     <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL, jackup){

      if( !is.null(bc_par[1]) ){
        y   <- bc(par=bc_par,y,jackup=jackup)
      }

      if( k > 0){
        y_ms<- c( mean(y), sd(y) )
        y   <- sc(par=y_ms,y)

        for(kk in 1:k) {
          nc_dum<-ifelse((noconst_last==TRUE)&(kk==k),TRUE,FALSE)
          y <- sal(par[[kk]],y,noconst=nc_dum)
        }
      } else {
        y_ms<- NULL
      }

      z_ms  <- c( mean(y), sd(y) )
      y     <- sc(par=z_ms,y)

      return(list(y=y,y_ms=y_ms, z_ms=z_ms))
    }

    ######## Negative log-likelihood (SAL distribution)
    NLL_sal<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){
      n   <-length(y)
      par2<-list(NULL)
      for(kk in 1:k){
        if(noconst_last & (kk==k)){
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+2)]

        } else {
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par[(4*(kk-1)+4)]<- abs( par[(4*(kk-1)+4)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+4)]
        }
      }

      if(tr_nonneg){
        np_b  <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp<- sum( log( d_abs ) )
      nll <- ee/2 - comp
    }

    NLL_sal2<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){
      n   <-length(y)
      par2<-list(NULL)
      for(kk in 1:k){
        if(noconst_last & (kk==k)){
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+2)]
        } else {
          par[(4*(kk-1)+1)]<- abs( par[(4*(kk-1)+1)] )
          par[(4*(kk-1)+4)]<- abs( par[(4*(kk-1)+4)] )
          par2[[kk]]        <- par[(4*(kk-1)+1):(4*(kk-1)+4)]
        }
      }

      if(tr_nonneg){
        np_b    <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp<- sum( log( d_abs ) )
      loglik<- ee/2 - comp
      return(list(z=z, b=b, loglik=loglik,comp=comp,y_ms=y_ms,z_ms=z_ms))
    }

    ######## Negative log-likelihood (BC distribution)
    NLL_bc <-function(par,y,M,Minv,m0,jackup){
      par[2]<-abs(par[2])
      z0  <- bc(par=par,y=y,jackup=jackup)
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]#z0#
      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,y ) )#abs( d_bc( par=par,y ) )#
      comp<- sum( log( d_abs ) )
      nll <- ee/2 - comp
    }
    NLL_bc2<-function(par,y, M, Minv,m0,jackup){
      par[2]<-abs(par[2])
      z0  <- bc(par=par,y=y,jackup=jackup)
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]#z0#
      m   <- t(m0)%*%z
      b   <- Minv%*%m
      ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      d_abs<-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,y ) )#abs( d_bc( par=par,y ) )#
      comp<- sum( log( d_abs ) )
      loglik<- ee/2 - comp
      return(list(z=z, b=b, loglik=loglik,comp=comp,z_ms=z_ms))#
    }

    lik_resf_vc_tr <- function( par, k, tr_nonneg=FALSE,noconst_last=TRUE, X,
                                evSqrt, y0, n, nx, emet,M0,Minv,term1,null_dum3=NULL,jackup ){
      tr_par       <- par
      if(k > 0){
        tr_par2<-list(NULL)
        for(kk in 1:k){
          if(noconst_last & (kk==k)){
            tr_par2[[kk]]        <- tr_par[(4*(kk-1)+1):(4*(kk-1)+2)]
            tr_par2[[kk]][1]     <- abs(tr_par2[[kk]][1])
          } else {
            tr_par2[[kk]]        <- tr_par[(4*(kk-1)+1):(4*(kk-1)+4)]
            tr_par2[[kk]][c(1,4)]<- abs(tr_par2[[kk]][c(1,4)])
          }
        }

        if(tr_nonneg){
          np_b  <-length(tr_par)
          bc_par<-tr_par[(np_b-1):np_b]
          #if(bc_par[1] < -0.5) bc_par[1] <- -0.5###################### added 2020/12/14
          bc_par[2]<-abs(bc_par[2])
        } else {
          bc_par<-NULL
        }
        z0      <- sal_k(par=tr_par2,y=y0,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup)
        z       <- z0$y
        z_ms    <- z0$z_ms
        y_ms    <- z0$y_ms
        d_abs   <- abs( d_sal_k( par=tr_par2 ,y=y0, k=k,noconst_last=noconst_last,
                                 bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup ) )
        comp    <- (-2)*sum( log( d_abs ))

        if(max(abs(z_ms)) > 10^10){######### deleted
          comp  <- 10^100
        }

      } else if(tr_nonneg){
        z0      <- bc(par=tr_par,y=y0,jackup=jackup)
        z_ms    <- c(mean(z0),sd(z0))
        z       <- (z0-z_ms[1])/z_ms[2]
        d_abs   <- abs( d_bc(par=tr_par,y0,jackup=jackup)*d_sc(par=z_ms,y0) )
        comp    <- (-2)*sum( log( d_abs ) )
      } else {
        z       <- y0
        comp    <- 0
      }

      m	        <- crossprod( X[,null_dum3], z )
      zz        <- sum( z^2 )

      ########################## reml
      m[-(1:nx)]		<- m[ -( 1:nx ) ] * evSqrt[null_dum3[-( 1:nx )]]
      b		<- Minv %*% m
      sse		<- zz - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
      dd		<- abs(sse) + sum( b[ -( 1:nx ) ] ^ 2 )
      if( emet == "reml" ){
        term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
      } else if( emet == "ml" ){
        term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
      }
      loglik		<- term1 + term2 + comp

      if(is.nan(loglik)) loglik <-10^100
      if(is.na(loglik))  loglik <-10^100

      return( loglik[ 1 ] )
    }

    lm_cw<-function(y, M, Minv, m0, k=2,noconst_last=TRUE,tr_nonneg=FALSE,jackup){

      if(k != 0){
        par00<-rep(c(1,0,0,1),k)
        #lower<-rep(c(1e-10,-3,-Inf,1e-10),k)###### check!!!
        #upper<-rep(c(Inf, 3, Inf, Inf),k)

        if(noconst_last){
          par00<-par00[1:(4*k-2)]
          #lower<-lower[1:(4*k-2)]
          #upper<-upper[1:(4*k-2)]
        }

        if(tr_nonneg){
          par00<-c(par00,1, 0.001)
          #lower<-c(lower, -1, 0)
          #upper<-c(upper,  5, 10)
        }
        res    <-optim(par=par00,fn=NLL_sal, y=y,M=M,Minv=Minv,m0=m0,k=k,
                       noconst_last=noconst_last,tr_nonneg=tr_nonneg,jackup=jackup)
        #method="L-BFGS-B",lower=lower,upper=upper)
        est0   <-res$par
        est    <-list(NULL)
        for(kk in 1:k){
          if(noconst_last & (kk==k)){
            est0[(4*(kk-1)+1)]<- abs( est0[(4*(kk-1)+1)] )
            est[[kk]]    <- est0[(4*(kk-1)+1):(4*(kk-1)+2)]
          } else {
            est0[(4*(kk-1)+1)]<- abs( est0[(4*(kk-1)+1)] )
            est0[(4*(kk-1)+4)]<- abs( est0[(4*(kk-1)+4)] )
            est[[kk]]    <- est0[(4*(kk-1)+1):(4*(kk-1)+4)]
          }
        }

        if(tr_nonneg){
          np_be      <- length(est0)
          bc_par     <- est0[(np_be-1):np_be]
          if(jackup) bc_par[2]  <- abs(bc_par[2])
          #bc_par     <- est0[length(est0)]
          est[[kk+1]]<- bc_par
        } else {
          bc_par<-NULL
        }
        res2  <-NLL_sal2(par=unlist(est), y=y, M=M, Minv=Minv, m0=m0,k=k,
                         noconst_last=noconst_last, tr_nonneg=tr_nonneg,jackup=jackup)
        b     <-res2$b
        z     <-res2$z
        loglik<-res2$loglik
        comp  <-res2$comp
        y_ms  <-res2$y_ms
        z_ms  <-res2$z_ms
        #if(tr_nonneg){
        #  est <- est[1:(length(est)-2)]
        #}

      } else if(tr_nonneg){
        #res   <-optimize(f=NLL_bc, interval=c(-5,5),y=y,M=M,Minv=Minv,m0=m0 )#c(-5,5)
        #bc_par<-res$minimum
        #est   <-NULL

        res   <-optim(par=c(1,0.001),fn=NLL_bc,y=y,M=M,Minv=Minv,m0=m0,jackup=jackup )# interval=c(-5,5)
        bc_par<-res$par
        bc_par[2]<-abs(bc_par[2])
        est   <-NULL

        res2  <-NLL_bc2(par=bc_par, y=y, M=M, Minv=Minv, m0=m0,jackup=jackup )
        b     <-res2$b
        z     <-res2$z
        loglik<-res2$loglik
        comp  <-res2$comp
        y_ms  <-NULL
        z_ms  <-res2$z_ms# this part is different from Kaihatsuchushi

      } else {
        est    <- NULL
        res2   <- NULL
        b      <- NULL
        z      <- y
        loglik <- NULL
        comp   <- 0
        bc_par <- NULL
        y_ms  <-NULL
        z_ms  <-NULL
      }

      return(list(vpar=est,bc_par=bc_par,b=b,z=z,loglik=loglik,comp=comp,
                  y_ms=y_ms,z_ms=z_ms))
    }
  }

  n<-length(mod$other$coords[,1])
  if( ( length(mod$other$evSqrts_n)==1 )&( is.null(mod$other$evSqrt_n[[1]]))&(length(mod$other$evSqrts) > 1 )){
    mod$other$evSqrts_n <-list(NULL)
    for(pp in 2:length(mod$other$evSqrts)){
      mod$other$evSqrts_n<-c(mod$other$evSqrts_n, list(NULL))
    }
  }


  if( !is.null( mod$other$x ) ){
    X1      <- as.matrix( as.matrix( mod$other$x ) [,mod$other$x_id] )
  }

   n0        <- dim(meig0$sf)[1]
   if( !is.null( x0 ) ){
     if(dim( as.matrix( x0 ) )[2] != length( mod$other$x_id )){
       stop("x and x0 must have the same number of columns")
     }
     if( n0 == 1){
       X0      <- as.matrix( x0 )[ mod$other$x_id]
     } else {
       X0      <- as.matrix( as.matrix( x0 )[,mod$other$x_id] )
     }
   }

   if( !is.null( mod$other$xconst ) ){
     X1const <- as.matrix( as.matrix( mod$other$xconst ) [,mod$other$xf_id] )
   }


   if( !is.null( xconst0 ) ){
     if( dim( as.matrix( xconst0 ) )[2] != length( mod$other$xf_id ) ){
       stop("xconst and xconst0 must have the same number of columns")
     }

     if( n0 == 1 ){
       X0const <- as.matrix( xconst0 )[ mod$other$xf_id]
     } else {
       X0const <- as.matrix( as.matrix( xconst0 )[,mod$other$xf_id] )
     }
   }

  XX1_0	      <- list(NULL)
  XX1	        <- NULL
  nvc_x     <- mod$other$nvc_x
  if( is.logical( nvc_x[1] ) == FALSE ){
    X1_nvc  <- as.matrix( as.matrix(X1)[ , nvc_x ] )

    if( n0 == 1 ){
      X0_nvc  <- as.matrix(X0)[ nvc_x ]
    } else {
      X0_nvc  <- as.matrix( as.matrix(X0)[ , nvc_x ] )
    }
    xxname	<- names( as.data.frame( X1 ) )[ nvc_x ]
    nnsv    <- length( xxname )
    np_xx   <-apply( X1_nvc,2,function( x ) length( unique( x )))
    np_xx   <-ifelse( np_xx < mod$other$nvc_num/0.7, round( np_xx * 0.7 ) ,mod$other$nvc_num )

    np_xx_max <-round( n/nnsv ) - 2
    np_xx[ np_xx > np_xx_max ] <-np_xx_max
    np_xx[ np_xx < 2 ] <- 2

    #XB_n    <-NULL
    B_n     <-list(NULL)
    for( ii in 1:dim( X1_nvc )[ 2 ] ){
      if( np_xx[ ii ] <= 2 ){
        B_n[[ii+1]]<- 0
        np_xx[ ii ]<- 0
      } else {

        test<-TRUE
        iiii<-0
        while( test ){
          kkk       <- np_xx[ ii ]-iiii
          knots     <- seq(min( X1_nvc[ , ii ] ),max( X1_nvc[ , ii ] ),len=kkk+2)[2:(kkk+1)]
          testt     <- try(XX1_00    <- ns( X1_nvc[ ,ii ], knots = knots ), silent=TRUE)#replaced
          test      <- class(testt)[1] == "try-error"
          iiii      <- iiii+1
        }

        XX1_0       <- cbind( X1_nvc[,ii], XX1_00)
        if( n0 == 1 ){
          XX0_0     <- cbind( X0_nvc[ ii], predict( XX1_00, newx = X0_nvc[ ii ] ) )
        } else {
          XX0_0     <- cbind( X0_nvc[,ii], predict( XX1_00, newx = X0_nvc[,ii ] ) )
        }

        #if( !is.na( mod$other$sel_basis_n[[ ii ]][ 1 ] ) ){
        #}

        if( n0 == 1 ){
          n_ids <-(1:length(c(XX0_0)))[(1:length(c(XX0_0))) %in% mod$other$sel_basis_n[[ii]]]
          for( j in length(c(XX0_0)) ){
            XX0_0[ j]<- ( XX0_0[ j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])
          }
        } else {
          n_ids <-(1:dim(XX0_0)[2])[(1:dim(XX0_0)[2]) %in% mod$other$sel_basis_n[[ii]]]
          for( j in  n_ids){
            XX0_0[,j]<- ( XX0_0[,j] - mean(XX1_0[,j]) )/sd(XX1_0[,j])
          }
        }

        if( n0 == 1 ){
          B_n[[ii+1]]<- XX0_0[ mod$other$sel_basis_n[[ii]] ]
        } else {
          B_n[[ii+1]]<- XX0_0[ ,mod$other$sel_basis_n[[ii]] ]
        }

      }
    }
  }

  XXconst_0	  <- list( NULL )
  XXconst	    <- NULL
  nvc_xconst <- mod$other$nvc_xconst
  if( is.logical( nvc_xconst[ 1 ] ) ==FALSE ){
    X1const_nvc<- as.matrix( X1const[ , nvc_xconst ] )
    if( n0 == 1 ){
      X0const_nvc<- c( X0const[ nvc_xconst ] )
    } else {
      X0const_nvc<- as.matrix( X0const[ , nvc_xconst ] )
    }

    xxfname	  <- names( as.data.frame( X1const ) )[ nvc_xconst ]
    nnxf      <- length( xxfname )
    np_xxconst<-apply( X1const_nvc,2,function( x ) length( unique( x )))
    np_xxconst<-ifelse( np_xxconst < mod$other$nvc_num/0.7, round( np_xxconst * 0.7 ) ,mod$other$nvc_num )
    np_xxconst_max <-round( n/nnxf ) - 2
    np_xxconst[ np_xxconst > np_xxconst_max ] <-np_xxconst_max
    np_xxconst[ np_xxconst < 2 ] <- 2

    #XB_c    <-NULL
    B_c     <-list(NULL)
    for( ii in 1:dim( X1const_nvc )[ 2 ] ){
     if( np_xxconst[ ii ] <= 2 ){
        B_c[[ii+1]]     <- 0
        np_xxconst[[ii]]<- 0
     } else {

       test<-TRUE
       iiii<-0
       while(test){
         kkk      <- np_xxconst[ ii ]-iiii
         knots    <-seq(min( X1const_nvc[ ,ii ] ),max( X1const_nvc[ ,ii ] ),len=kkk+2)[2:(kkk+1)]
         testt<- try(XX1const_00<- ns( X1const_nvc[ , ii], knots = knots ), silent=TRUE)
         test <- class(testt)[1] == "try-error"
         iiii <- iiii+1
       }

       XX1const_0   <- cbind( X1const[,ii] , XX1const_00)
       if( n0 == 1 ){
         XX0const_0 <- cbind( X0const_nvc[ ii], predict( XX1const_00, newx= X0const_nvc[ ii]) )
       } else {
         XX0const_0 <- cbind( X0const_nvc[,ii], predict( XX1const_00, newx= X0const_nvc[,ii]) )
       }

       if( !is.na( mod$other$sel_basis_c[[ ii ]][ 1 ] ) ){
         XX1const_0<- XX1const_0[,mod$other$sel_basis_c[[ii]]]
         if( n0 == 1){
           XX0const_0<- XX0const_0[ mod$other$sel_basis_c[[ii]]]
         } else {
           XX0const_0<- XX0const_0[,mod$other$sel_basis_c[[ii]]]
         }
       }

       if( n0 == 1 ){
         for( j in 1:length(XX0const_0)){
           XX0const_0[ j]<- ( XX0const_0[ j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])
         }
       } else {
         for( j in 1:dim(XX0const_0)[2]){
           XX0const_0[,j]<- ( XX0const_0[,j] - mean(XX1const_0[,j]))/sd(XX1const_0[,j])
         }
       }

       B_c[[ii]]    <- XX0const_0
       #XB_c         <- cbind( XB_c , X0const_nvc[, ii ] * B_c[[ii]] )
     }
    }
  }

  n	  <- length( mod$other$y )
  n0	<- length( meig0$sf[ , 1 ] )
  ne	<- length( mod$other$b_s[[ 1 ]][ -1 ] )
  nsv	<- sum( mod$other$x_id ) + 1
  nxf	<- sum( mod$other$xf_id )
  meig0$sf<- as.matrix( meig0$sf[ , 1:ne ])

  #for( i in 1:length(mod$other$evSqrts_n)){
  #  if(length( mod$other$evSqrts_n[[i]] ) == 1) mod$other$evSqrts_n[[i]] <- NULL
  #}
  #for( i in 1:length(mod$other$evSqrts)){
  #  if(length( mod$other$evSqrts[[i]] ) == 1) mod$other$evSqrts[[i]] <- NULL
  #}

  b_vc	<- matrix(0, nrow = n0, ncol = nsv )
  bse_vc<- matrix(0, nrow = n0, ncol = nsv )
  bt_vc	<- matrix(0, nrow = n0, ncol = nsv )
  bp_vc	<- matrix(0, nrow = n0, ncol = nsv )
  for( i in 1:nsv ){
    evSqrts    <-mod$other$evSqrts[[ i ]]
    if(length( evSqrts ) == 1) evSqrts <- NULL

    if(i ==1){
      evSqrts_n<- NULL
    } else {
      evSqrts_n<- mod$other$evSqrts_n[[ i ]]
      if(length( evSqrts_n ) == 1) evSqrts_n <- NULL
    }

    evSqrt_i        <-c(evSqrts, evSqrts_n)
    if( length( evSqrt_i ) <= 1 ){
      b_vc[ , i ]	  <- mod$other$b_s[[ i ]][ 1 ]
      if( i ==1 ) b_vc[ , i ]<-b_vc[ , i ] + mod$other$Bias      #### added

      bse_vc[ , i ]	<- sqrt( mod$other$b_covs[[ i ]] )
      bt_vc[ , i ]	<- b_vc[ , i ] / bse_vc[ , i ]
      bp_vc[ , i ]	<- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    } else {
      if( n0 == 1 ){
        if( is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + c(meig0$sf) %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( meig0$sf * evSqrts )
        } else if( !is.null( evSqrts_n ) & is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + c(B_n[[i]]) %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( B_n[[i]] * evSqrts_n )
        } else if( !is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          basis2     <- c( meig0$sf, B_n[[i]] )
          evSvqr2    <- c( evSqrts, evSqrts_n )
          sf2		     <- t( basis2 * evSvqr2 )
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + basis2 %*% c(  mod$other$b_s[[ i ]][ -1 ] )
        }

      } else {
        if( is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + meig0$sf %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( meig0$sf ) * evSqrts )
        } else if( !is.null( evSqrts_n ) & is.null( evSqrts ) ){
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + B_n[[i]] %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( B_n[[i]] ) * evSqrts_n )
        } else if( !is.null( evSqrts_n ) & !is.null( evSqrts ) ){
          basis2     <- cbind( meig0$sf, B_n[[i]] )
          evSqrts2    <- c( evSqrts, evSqrts_n )
          b_vc[ , i ]<- mod$other$b_s[[ i ]][ 1 ] + basis2 %*% c( mod$other$b_s[[ i ]][ -1 ] )
          sf2		     <- t( t( basis2 ) * evSqrts2 )
        }
      }

      if( i ==1 ) b_vc[ , i ] <- b_vc[ , i ] + mod$other$Bias      #### added

      x_sf		        <- as.matrix( cbind( 1, sf2 ) )
      bse_vc[ , i ]	  <- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs[[ i ]] %*% t( x_sf ) ) ) )
      bt_vc[ , i ]	  <- b_vc[ , i ] / bse_vc[ , i ]
      bp_vc[ , i ]	  <- 2 - 2 * pt( abs( bt_vc[ , i ] ), df = n - mod$other$df )
    }
  }

  if( is.null( mod$c_vc ) ){
    c_vc  <- cse_vc <- ct_vc <- cp_vc <- NULL
  } else {
    xc_vc	<- 0
    c_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    cse_vc<- matrix(0, nrow = n0, ncol = nnxf )
    ct_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    cp_vc	<- matrix(0, nrow = n0, ncol = nnxf )
    for( i in 1:nnxf ){
      evSqrts_c<- mod$other$evSqrts_c[[ i ]]
      if(length( evSqrts_c ) == 1) evSqrts_c <- NULL

      if( length( mod$other$evSqrts_c[[ i ]] ) <= 1 ){#!= ne
        c_vc[ , i ]	  <- mod$other$b_c[[ i ]][ 1 ]
        cse_vc[ , i ]	<- sqrt( mod$other$b_covs_c[[ i ]] )
        ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
        cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$df )
      } else {
        if( n0 == 1 ){
          c_vc[ , i ]	<- mod$other$b_c[[ i ]][1] + c(B_c[[ i ]]) %*% c( mod$other$b_c[[ i ]][ -1 ])
          sf2		<- t( B_c[[ i ]] * mod$other$evSqrts_c[[ i ]] )
        } else {
          c_vc[ , i ]	<- mod$other$b_c[[ i ]][1] + B_c[[ i ]] %*% c( mod$other$b_c[[ i ]][ -1 ] )
          sf2		<- t( t( B_c[[ i ]] ) * mod$other$evSqrts_c[[ i ]] )
        }
        x_sf		      <- as.matrix( cbind( 1, sf2 ) )
        cse_vc[ , i ]	<- sqrt( colSums( t( x_sf ) * ( mod$other$b_covs_c[[ i ]] %*% t( x_sf ) ) ) )
        ct_vc[ , i ]	<- c_vc[ , i ] / cse_vc[ , i ]
        cp_vc[ , i ]	<- 2 - 2 * pt( abs( ct_vc[ , i ] ), df = n - mod$other$df )
      }
    }
  }


  if( !is.null( mod$b_g ) ){
    if( is.null(xgroup0) ){
      message( "Note: Group effects are ignored because xgroup0 is missing")
      b_g0     <- NULL
    } else {
      ng       <- length(mod$b_g)
      xgroup0  <- data.frame( xgroup0 )
      b_g0     <- xgroup0
      for(ggid in 1:ng){
        b_g_sub        <- as.data.frame( mod$b_g[[ggid]] )
        xg00_0	       <- factor( xgroup0[ , ggid ] )
        xg00_nam0      <- names( xgroup0 )[ ggid ]
        xg00_nam       <-ifelse( xg00_nam0 == "xgroup0", "xgroup", xg00_nam0 )
        xgroup0_nam    <- paste( xg00_nam, "_", xg00_0, sep = "" )
        xgroup0_dat    <- data.frame(id=1:length(xgroup0_nam),xg_nam=xgroup0_nam)
        xgroup0_dat2   <- merge(xgroup0_dat, b_g_sub, by.x = "xg_nam", by.y = "row.names", all.x = TRUE )
        b_g0[,ggid]    <- xgroup0_dat2[order(xgroup0_dat2$id), "Estimate" ]
      }
      xg_names         <- names( xgroup0 )
      if( length( xg_names ) == 1 ){
        xg_names       <- ifelse( xg_names == "xgroup0", "xgroup", paste("xgroup_",xg_names,sep="" ) )
      } else {
        xg_names       <- paste( "xgroup_", xg_names, sep = "" )
      }
      names(b_g0)      <- xg_names

    }
  } else {
    b_g0 <- NULL
  }
  b_g    <- b_g0


  xb_vc	<- 0
  if( is.null( x0 ) == FALSE ){
    if( length( mod$other$x_id ) == 1 ){
      x0	<- as.matrix( cbind( 1, x0[  mod$other$x_id ] ))
    } else {
      x0	<- as.matrix( cbind( 1, x0[, mod$other$x_id ] ))
    }
    if( is.numeric( x0 ) == FALSE ){
      mode( x0 )	<- "numeric"
    }
  }

  if( is.null( mod$other$xf_id ) ){
    xb_const	<- 0
    if( !is.null( xconst0 ) ) message( "Note: xconst0 is ignored" )
  } else {
    if( is.null( xconst0 ) ){
      xb_const	<- 0
    } else {
      if( length( mod$other$xf_id ) == 1 ){
        xconst0	<- as.matrix( xconst0[  mod$other$xf_id ] )
      } else {
        xconst0	<- as.matrix( xconst0[, mod$other$xf_id ] )
      }
      if( is.numeric( xconst0 ) == FALSE ){
        mode( xconst0 )	<- "numeric"
      }
    }
  }

  if( n0 == 1 ){
    xb_vc   <- sum(x0[-1]*b_vc[-1])
    sf_pred <- b_vc[ 1]
  } else {
    xb_vc   <- rowSums(as.matrix(x0[,-1]*b_vc[,-1]))
    sf_pred <- b_vc[,1]
  }

  if( is.null(xconst0)){
    xb_const  <- 0
  } else {
    if( !is.null( mod$c_vc ) ){
      if( n0 == 1 ){
        xb_const  <- sum( xconst0 * c_vc )
      } else {
        xb_const  <- rowSums( xconst0 * c_vc )
      }
    } else {
      if( !is.null( mod$c ) ){
        if( is.null( dim( mod$c ) ) ){
          xb_const<- xconst0 * mod$c[  1 ]
        } else {
          xb_const<- xconst0 %*% mod$c[, 1 ]
        }
      } else {
        xb_const  <- 0
      }
    }
  }

  xconst_problem <- !is.null( mod$other$xf_id ) & is.null( xconst0 )
  if( ( !is.null( x0 ) ) & ( xconst_problem == FALSE ) ){
    	sf_pred	<- sf_pred - mean( mod$b_vc[,1] )
    	xb	    <- xb_vc + xb_const + mean( mod$b_vc[,1] )
    	pred0	<- xb + sf_pred
    	if( is.null( b_g0 ) == FALSE ){
    		na_b_g0   <- is.na( b_g0 )
    		if( sum( na_b_g0 ) > 0 ){
    		  b_g0[ na_b_g0 ] <- 0
    		  message( "Note: NAs are given to the groups that are not in xgroup")
    		}
    		pred0   <- pred0 + rowSums( b_g0 )
    	}

    	y0         <- mod$other$y
    	tr_num     <- mod$other$tr_num
    	tr_nonneg  <- mod$other$tr_nonneg
    	tr_par     <- mod$tr_par
    	tr_bpar    <- mod$tr_bpar$Estimate
    	y_added    <- mod$other$y_added
    	jackup     <- mod$other$jackup
    	noconst_last<-TRUE
    	if( tr_num > 0 ){######## transfer this part to prediction functions
    	  z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,jackup=jackup)
    	  z_ms     <- z0$z_ms
    	  y_ms     <- z0$y_ms

    	  dif      <- 1/d_sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)
    	  pred     <- i_sal_k(par=tr_par,y=pred0,k=tr_num,noconst_last=noconst_last,
    	                      bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup) - y_added
    	  if( tr_nonneg ) pred[ pred < 0 ] <- 0

    	} else if( tr_nonneg ){
    	  z0      <- bc(par=tr_bpar,y=y0,jackup=jackup)
    	  z_ms    <- c(mean(z0),sd(z0))
    	  y_ms     <- NULL

    	  dif      <- 1 / ( d_bc(par=tr_bpar,y=y0,jackup=jackup)/z_ms[2] )
    	  pred0b   <- pred0*z_ms[2] + z_ms[1]
    	  pred     <- i_bc(par=tr_bpar,y=pred0b,jackup=jackup) - y_added
    	  pred[ pred < 0 ] <- 0
    	} else {
    	  dif      <- 1
    	  z_ms     <- NULL
    	  y_ms     <- NULL
    	  pred     <- pred0
    	}

    	if( is.null(z_ms) ){
    	  if( !is.null(xgroup0) ){
    	    res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred, "group" = b_g )
    	  } else {
    	    res	<- data.frame( "pred" = pred, "xb" = xb, "sf_residual" = sf_pred )
    	  }

    	} else {
    	  if( !is.null(xgroup0) ){
    	    res	<- data.frame( "pred" = pred, "pred_trans" = pred0, "xb" = xb, "sf_residual" = sf_pred, "group" = b_g )
    	  } else {
    	    res	<- data.frame( "pred" = pred, "pred_trans" = pred0, "xb" = xb, "sf_residual" = sf_pred )
    	  }
    	}

  } else if( is.null( x0 ) & ( xconst_problem == FALSE )){
    #sf_pred	<- sf_pred - mean( mod$b_vc[,1] )
    #xb	    <- xb_vc + xb_const + mean( mod$b_vc[,1] )
    #pred	<- xb + sf_pred
    #if( is.null( b_g0 ) == FALSE ){
    #  na_b_g0   <- is.na( b_g0 )
    #  if( sum( na_b_g0 ) > 0 ){
    #    b_g0[ na_b_g0 ] <- 0
    #    message( "Note: NAs are given to the groups that are not in xgroup")
    #  }
    #  pred   <- pred + rowSums( b_g0 )
    #}
    res	<- NULL
    message( "Note: y is not predicted because x0 and/or xconst0 is missing" )
  } else {
    	res	<- NULL
    	message( "Note: y is not predicted because x0 and/or xconst0 is missing" )
  }

  b_vc	  <- as.data.frame( b_vc )
  bse_vc	<- as.data.frame( bse_vc )
  bt_vc	  <- as.data.frame( bt_vc )
  bp_vc	  <- as.data.frame( bp_vc )
  names( b_vc )	 <- names( mod$b_vc )
  names( bse_vc )<- names( mod$b_vc )
  names( bt_vc ) <- names( mod$b_vc )
  names( bp_vc ) <- names( mod$b_vc )

  c_vc	  <- as.data.frame( c_vc )
  cse_vc	<- as.data.frame( cse_vc )
  ct_vc	  <- as.data.frame( ct_vc )
  cp_vc	  <- as.data.frame( cp_vc )

  if( dim(c_vc)[1] > 0){
    names( c_vc )	 <- names( mod$c_vc )
    names( cse_vc )<- names( mod$c_vc )
    names( ct_vc ) <- names( mod$c_vc )
    names( cp_vc ) <- names( mod$c_vc )
  } else {
    c_vc <- cse_vc <- ct_vc <- cp_vc <- NULL
  }

  result	<- list( pred = res, b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc,
                  c_vc = c_vc, cse_vc = cse_vc, ct_vc = ct_vc, cp_vc = cp_vc )

 return( result )
}

