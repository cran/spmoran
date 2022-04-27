
resf_vc	  <- function( y, x, xconst = NULL, xgroup = NULL, weight = NULL, offset = NULL,
                       x_nvc = FALSE, xconst_nvc = FALSE,
                       x_sel = TRUE, x_nvc_sel = TRUE, xconst_nvc_sel = TRUE, nvc_num = 5,
                       meig, method = "reml", penalty = "bic", maxiter = 30, nongauss = NULL ){


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

    d_af <-function(par,y) par[2]
    d_sc <-function(par,y) rep(1/par[2], length(y))###
    d_sa <-function(par,y) par[1]*cosh(par[1]*asinh(y)-par[2])/sqrt(1+y^2)
    d_bc <-function(par,y,jackup) (y + ifelse(jackup, abs( par[2] ), 0 ) )^(par[1]-1)#abs(y)^(par-1)
    d_sal<-function(par,y,noconst=FALSE){
      if(noconst==FALSE){
        par[4]*d_sa(par=par[1:2],y)#d_af(par=par[3:4],sa(par=par[1:2],y))
      } else {
        d_sa(par=par[1:2],y)# d_af(par=c(0, 1),sa(par=par[1:2],y))*
      }
    }

    i_af<-function(par,y) (y-par[1])/par[2]
    i_sc <-function(par,y) par[1]+par[2]*y
    i_sa<-function(par,y) sinh(1/par[1]*(asinh(y)+par[2]))
    i_bc<-function(par,y,jackup){#, ulim=NULL
      y_pre <- y
      if(par[1]==0){
        y   <-exp(y) - ifelse(jackup, abs( par[2] ), 0 )
      } else {
        y   <- (par[1]*y + 1)^(1/par[1]) - ifelse(jackup, abs( par[2] ), 0 )
      }

      #if( !is.null(ulim) ){
      #  ulim_id <- which(y_pre > ulim)
      #  if( length( ulim_id > 1) ){
      #    y_lower    <-sort(unique(y[-ulim_id]),decreasing=TRUE)[2:1]
      #    y_pre_lower<-sort(unique(y_pre[-ulim_id]),decreasing=TRUE)[2:1]
      #    grad_lower <-diff(y_lower)/diff(y_pre_lower)
      #    y[ulim_id] <-y_lower[2]+(y_pre[ulim_id] - y_pre_lower[2])*grad_lower
      #  }
      #}

      if( sum( is.na(y) )>0 ) y[is.na(y)]<-min(y[!is.na(y)]/2)
      return(y)
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
    i_sal_k   <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL,y_ms=NULL,z_ms, jackup,
                         y_added2=0){#,ulim=NULL
      y_pre   <-y
      y       <-i_sc(par=z_ms,y)
      y       <-y - y_added2
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

      #if( !is.null(ulim) ){
      #  ulim_id <- which(y_pre > ulim)
      #  if( length( ulim_id > 1) ){
      #    y_ulim     <- y[ulim_id]
      #    y_pre_ulim <-y_pre[ulim_id]
      #    y_lower    <-sort(y[-ulim_id],decreasing=TRUE)[2:1]
      #    y_pre_lower<-sort(y_pre[-ulim_id],decreasing=TRUE)[2:1]
      #    grad_lower <-diff(y_lower)/diff(y_pre_lower)
      #    y[ulim_id] <-y_lower[2]+(y_pre_ulim - y_pre_lower[2])*grad_lower
      #  }
      #}

      if( sum( is.na(y) )>0 ) y[is.na(y)]<-min(y[!is.na(y)]/2)
      return(y)
    }
    sal_k     <-function(par,y,k=2,noconst_last=TRUE,bc_par=NULL, jackup,y_added2=0){

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
      y     <- y + y_added2

      z_ms  <- c( mean(y), sd(y) )
      y     <- sc(par=z_ms,y)

      return(list(y=y,y_ms=y_ms, z_ms=z_ms))
    }

    ######## Negative log-likelihood (SAL distribution)
    NLL_sal<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,y_nonneg=FALSE,jackup,weight=NULL,plim,
                      bc_value=NULL,y_added2=0, y_type){
      n   <-length(y)
      par[par > plim]  <-  plim
      par[par < -plim] <- -plim

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

      if(y_nonneg){
        np_b  <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
        if( !is.null( bc_value ) ){
          bc_par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
        } else {
          bc_par[1]<- ifelse( bc_par[1] < -5, -5, ifelse( bc_par[1] > 5, 5, bc_par[1]) )
        }

      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup,
                   y_added2=y_added2)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      #if(1<0){
      steepen <-0
      if( y_type=="count" ){
        yz0    <-cbind(y,z)

        yz_ord <-yz0[order(yz0[,1]),]
        yz     <-yz_ord[!duplicated(yz_ord[,1]),]
        uni_y  <-yz[,1];duni_y<-diff(uni_y)
        uni_z  <-yz[,2];duni_z<-diff(uni_z)
        duni_n <-length(duni_y)
        steep1a<-max( (duni_z/duni_y)[-duni_n]/(duni_z/duni_y)[-1] )
        steep1b<-max( (duni_z/duni_y)[-1]/(duni_z/duni_y)[-duni_n] )

        yz_ord <-yz0[order(yz0[,2]),]
        yz     <-yz_ord[!duplicated(yz_ord[,2]),]
        uni_y  <-yz[,1];duni_y<-diff(uni_y)
        uni_z  <-yz[,2];duni_z<-diff(uni_z)
        duni_n <-length(duni_y)
        steep1c<-max( (duni_y/duni_z)[-duni_n]/(duni_y/duni_z)[-1] )
        steep1d<-max( (duni_y/duni_z)[-1]/(duni_y/duni_z)[-duni_n] )

        steep1 <-max( c(steep1a, steep1b) )#steep1c, steep1d

        if( length(duni_y) >=4 ){
          steep2   <- (duni_z[2]/duni_y[2]) / (duni_z[3]/duni_y[3])
          steep3   <- (duni_z[3]/duni_y[3]) / (duni_z[4]/duni_y[4])
          dsteep123<- c(steep2/steep1, steep3/steep2)
          dsteepen <- max(dsteep123)/min(dsteep123)
          dsteep   <- (steep2 > steep1 & steep3 < steep2)|(steep2 < steep1 & steep3 > steep2)
        } else {
          dsteep   <- FALSE
        }

        duni_y_ex_last<-diff( uni_z[(duni_n):(duni_n+1)] )
        duni_y_last   <-duni_y[duni_n]
        dsteep_last   <-duni_y_last/duni_y_ex_last

        steepen<-ifelse( steep1 < 10,0,10^10 + steep1^3) +
          ifelse( dsteep==FALSE,0,10^10 + dsteepen^3) +
          ifelse( dsteep_last < 1.5, 0, 10^10 + dsteep_last^2)
      }
      #}

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      if( !is.null(weight)){
        ee	<- sum( weight*z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      } else {
        ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      }
      d_abs<-abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp<- sum( log( d_abs ) )
      nll <- n * ( 1 + log( 2 * pi * ee / n ) ) -2*comp + steepen
    }

    NLL_sal2<-function(par,y,M,Minv,m0,k=2,noconst_last=TRUE,y_nonneg=FALSE,jackup,weight=NULL,
                       plim, bc_value=NULL,y_added2=0){
      par[par > plim]  <-  plim
      par[par < -plim] <- -plim

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

      if(y_nonneg){
        np_b    <-length(par)
        bc_par<-par[(np_b-1):np_b]
        bc_par[2]<-abs(bc_par[2])
        if( !is.null( bc_value ) ){
          bc_par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
        } else {
          bc_par[1]<- ifelse( bc_par[1] < -5, -5, ifelse( bc_par[1] > 5, 5, bc_par[1]) )
        }

      } else {
        bc_par<-NULL
      }

      z0  <- sal_k(par=par2,y=y,k=k,noconst_last=noconst_last,bc_par=bc_par,jackup=jackup,
                   y_added2=y_added2)
      z   <- z0$y
      z_ms<- z0$z_ms
      y_ms<- z0$y_ms

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      if( !is.null(weight)){
        ee	<- sum( weight*z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      } else {
        ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      }
      d_abs <- abs( d_sal_k(par=par2,y,k=k,noconst_last=noconst_last,bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
      comp  <- sum( log( d_abs ) )
      #loglik<- ee/2 - comp
      nll   <- n * ( 1 + log( 2 * pi * ee / n ) ) + (-2)*comp
      loglik<- (-1/2)*nll
      return(list(z=z, b=b, loglik=loglik,comp=comp,y_ms=y_ms,z_ms=z_ms))
    }

    ######## Negative log-likelihood (BC distribution)
    NLL_bc <-function(par,y,M,Minv,m0,jackup,weight=NULL,bc_value=NULL,y_added2=0){
      n     <-length(y)
      par[2]<-abs(par[2])
      if( !is.null( bc_value ) ){
        par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
      } else {
        par[1]<- ifelse( par[1] < -5, -5, ifelse( par[1] > 5, 5, par[1]) )
      }
      z0  <- bc(par=par,y=y,jackup=jackup) + y_added2
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]

      #if(1<0){
      steepen <-0
      if( y_type=="count" ){
        yz0   <-cbind(y,z)

        yz_ord<-yz0[order(yz0[,1]),]
        yz    <-yz_ord[!duplicated(yz_ord[,1]),]
        uni_y <-yz[,1];duni_y<-diff(uni_y)
        uni_z <-yz[,2];duni_z<-diff(uni_z)
        duni_n<-length(duni_y)
        steep1a<-max( (duni_z/duni_y)[-duni_n]/(duni_z/duni_y)[-1] )
        steep1b<-max( (duni_z/duni_y)[-1]/(duni_z/duni_y)[-duni_n] )

        yz_ord<-yz0[order(yz0[,2]),]
        yz    <-yz_ord[!duplicated(yz_ord[,2]),]
        uni_y <-yz[,1];duni_y<-diff(uni_y)
        uni_z <-yz[,2];duni_z<-diff(uni_z)
        duni_n<-length(duni_y)
        steep1c<-max( (duni_y/duni_z)[-duni_n]/(duni_y/duni_z)[-1] )
        steep1d<-max( (duni_y/duni_z)[-1]/(duni_y/duni_z)[-duni_n] )

        steep1 <-max( c(steep1a, steep1b) )#,steep1c, steep1d

        if( length(duni_y) >=4 ){
          steep2   <- (duni_z[2]/duni_y[2]) / (duni_z[3]/duni_y[3])
          steep3   <- (duni_z[3]/duni_y[3]) / (duni_z[4]/duni_y[4])
          dsteep123<- c(steep2/steep1, steep3/steep2)
          dsteepen <- max(dsteep123)/min(dsteep123)
          dsteep   <- (steep2 > steep1 & steep3 < steep2)|(steep2 < steep1 & steep3 > steep2)
        } else {
          dsteep   <- FALSE
        }

        duni_y_ex_last<-diff( uni_z[(duni_n):(duni_n+1)] )
        duni_y_last   <-duni_y[duni_n]
        dsteep_last   <-duni_y_last/duni_y_ex_last

        steepen<-ifelse( steep1 < 10,0,10^10 + steep1^3) +
          ifelse( dsteep==FALSE,0,10^10 + dsteepen^3) +
          ifelse( dsteep_last < 1.5, 0, 10^10 + dsteep_last^2)
      }
      #}

      m   <- t(m0)%*%z
      b   <- Minv%*%m
      if( !is.null(weight)){
        ee	<- sum( weight*z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      } else {
        ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      }
      d_abs<-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,z0 ) )# y -> z0 ####### Fix this ############################
      comp<- sum( log( d_abs ) )
      nll   <- n * ( 1 + log( 2 * pi * ee / n ) )  -2 * comp + steepen
    }

    NLL_bc2<-function(par,y, M, Minv,m0,jackup,weight=NULL,bc_value=NULL,y_added2=0){
      n     <-length(y)
      par[2]<-abs(par[2])
      if( !is.null( bc_value ) ){
        par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
      } else {
        par[1]<- ifelse( par[1] < -5, -5, ifelse( par[1] > 5, 5, par[1]) )
      }

      z0  <- bc(par=par,y=y,jackup=jackup) + y_added2
      z_ms<- c(mean(z0),sd(z0))
      z   <- (z0-z_ms[1])/z_ms[2]#z0#
      m   <- t(m0)%*%z
      b   <- Minv%*%m
      if( !is.null(weight)){
        ee	<- sum( weight*z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      } else {
        ee	<- sum( z ^ 2 ) - 2 * t( b ) %*% m + t( b ) %*% M %*% b
      }
      d_abs <-abs( d_bc( par=par,y,jackup=jackup )*d_sc( par=z_ms,z0 ) )# y->z0 ###### Fix this ############################
      comp  <- sum( log( d_abs ) )
      nll   <- n * ( 1 + log( 2 * pi * ee / n ) ) -2*comp
      loglik<- (-1/2)*nll
      return(list(z=z, b=b, loglik=loglik,comp=comp,z_ms=z_ms))#
    }

    lik_resf_vc_tr <- function( par, k, y_nonneg=FALSE,noconst_last=TRUE, X, weight=NULL, bc_value=NULL,
                                evSqrt, y0, n, nx, emet,M0,Minv,term1,null_dum3=NULL,jackup, plim,
                                y_added2,y_type){
      tr_par                 <-  par
      tr_par[tr_par > plim]  <-  plim
      tr_par[tr_par < -plim] <- -plim

      steepen <-0
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

        if(y_nonneg){
          np_b  <-length(tr_par)
          bc_par<-tr_par[(np_b-1):np_b]
          #if(bc_par[1] < -0.5) bc_par[1] <- -0.5###################### added 2020/12/14
          bc_par[2]<-abs(bc_par[2])
          if( !is.null( bc_value ) ){
            bc_par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
          } else {
            bc_par[1]<- ifelse( bc_par[1] < -5, -5, ifelse( bc_par[1] > 5, 5, bc_par[1]) )
          }

        } else {
          bc_par<-NULL
        }
        z0      <- sal_k(par=tr_par2,y=y0,k=k,noconst_last=noconst_last,bc_par=bc_par,
                         jackup=jackup,y_added2 = y_added2)
        z       <- z0$y
        z_ms    <- z0$z_ms
        y_ms    <- z0$y_ms

        #if(1<0){
        steepen <-0
        if( y_type=="count" ){
          yz0   <-cbind(y0,z)

          yz_ord<-yz0[order(yz0[,1]),]
          yz    <-yz_ord[!duplicated(yz_ord[,1]),]
          uni_y <-yz[,1];duni_y<-diff(uni_y)
          uni_z <-yz[,2];duni_z<-diff(uni_z)
          duni_n<-length(duni_y)
          steep1a<-max( (duni_z/duni_y)[-duni_n]/(duni_z/duni_y)[-1] )
          steep1b<-max( (duni_z/duni_y)[-1]/(duni_z/duni_y)[-duni_n] )

          yz_ord<-yz0[order(yz0[,2]),]
          yz    <-yz_ord[!duplicated(yz_ord[,2]),]
          uni_y <-yz[,1];duni_y<-diff(uni_y)
          uni_z <-yz[,2];duni_z<-diff(uni_z)
          duni_n<-length(duni_y)
          steep1c<-max( (duni_y/duni_z)[-duni_n]/(duni_y/duni_z)[-1] )
          steep1d<-max( (duni_y/duni_z)[-1]/(duni_y/duni_z)[-duni_n] )
          steep1 <-max( c(steep1a, steep1b) )#,steep1c, steep1d

          if( length(duni_y) >=4 ){
            steep2   <- (duni_z[2]/duni_y[2]) / (duni_z[3]/duni_y[3])
            steep3   <- (duni_z[3]/duni_y[3]) / (duni_z[4]/duni_y[4])
            dsteep123<- c(steep2/steep1, steep3/steep2)
            dsteepen <- max(dsteep123)/min(dsteep123)
            dsteep   <- (steep2 > steep1 & steep3 < steep2)|(steep2 < steep1 & steep3 > steep2)
          } else {
            dsteep   <- FALSE
          }

          duni_y_ex_last<-diff( uni_z[(duni_n):(duni_n+1)] )
          duni_y_last   <-duni_y[duni_n]
          dsteep_last   <-duni_y_last/duni_y_ex_last

          steepen<-ifelse( steep1 < 10,0,10^10 + steep1^3) +
            ifelse( dsteep==FALSE,0,10^10 + dsteepen^3) +
            ifelse( dsteep_last < 1.5, 0, 10^10 + dsteep_last^2)
        }
        #}

        d_abs   <- abs( d_sal_k( par=tr_par2 ,y=y0, k=k,noconst_last=noconst_last,
                                 bc_par=bc_par,y_ms=y_ms,z_ms=z_ms,jackup=jackup ) )
        comp    <- (-2)*sum( log( d_abs ))

        if( is.na(max(abs(z_ms))) ){
          comp  <- 10^100
        } else if(max(abs(z_ms)) > 10^10){######### deleted
          comp  <- 10^100
        }

      } else if(y_nonneg){
        if( !is.null( bc_value ) ){
          tr_par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
        } else {
          tr_par[1]<- ifelse( tr_par[1] < -5, -5, ifelse( tr_par[1] > 5, 5, tr_par[1]) )
        }

        z0      <- bc(par=tr_par,y=y0,jackup=jackup) + y_added2
        z_ms    <- c(mean(z0),sd(z0))
        z       <- (z0-z_ms[1])/z_ms[2]
        d_abs   <- abs( d_bc(par=tr_par,y0,jackup=jackup)*d_sc(par=z_ms,z0) )#y0 ->z0###### Fix this ############################
        comp    <- (-2)*sum( log( d_abs ) )
      } else {
        z       <- y0
        comp    <- 0
      }

      if( !is.null(weight)){
        m	        <- crossprod( X[,null_dum3], sqrt( weight )*z )#X is sqrt(weight)*X
        zz        <- sum( weight*z^2 )
      } else {
        m	        <- crossprod( X[,null_dum3], z )
        zz        <- sum( z^2 )
      }

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
      loglik		<- term1 + term2 + comp + steepen

      if(is.nan(loglik)) loglik <-10^100
      if(is.na(loglik))  loglik <-10^100

      return( loglik[ 1 ] )
    }

    lm_cw<-function(y, M, Minv, m0, k=2,noconst_last=TRUE,y_nonneg=FALSE,jackup,weight=NULL,
                    plim, bc_value=NULL,y_added2=0){

      if(k != 0){
        par00<-rep(c(1,0,0,1),k)
        #lower<-rep(c(1e-10,-3,-Inf,1e-10),k)###### check!!!
        #upper<-rep(c(Inf, 3, Inf, Inf),k)

        if(noconst_last){
          par00<-par00[1:(4*k-2)]
          #lower<-lower[1:(4*k-2)]
          #upper<-upper[1:(4*k-2)]
        }

        if(y_nonneg){
          if(!is.null(bc_value)) {
            par00<-c(par00,bc_value, 0.001)
          } else {
            par00<-c(par00,1, 0.001)
          }
          #lower<-c(lower, -1, 0)
          #upper<-c(upper,  5, 10)
        }
        res    <-optim(par=par00,fn=NLL_sal, y=y,M=M,Minv=Minv,m0=m0,k=k,weight=weight,
                       noconst_last=noconst_last,y_nonneg=y_nonneg,jackup=jackup,
                       plim=plim,bc_value=bc_value,y_added2=y_added2,y_type=y_type)
        #method="L-BFGS-B",lower=lower,upper=upper)
        est0   <-res$par
        est0[est0 >  plim] <-  plim
        est0[est0 < -plim] <- -plim

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

        if(y_nonneg){
          np_be      <- length(est0)
          bc_par     <- est0[(np_be-1):np_be]
          if(jackup) bc_par[2]  <- abs(bc_par[2])
          if( !is.null( bc_value ) ){
            bc_par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
          } else {
            bc_par[1]<- ifelse( bc_par[1] < -5, -5, ifelse( bc_par[1] > 5, 5, bc_par[1]) )
          }
          est[[kk+1]]<- bc_par
        } else {
          bc_par<-NULL
        }
        res2  <-NLL_sal2(par=unlist(est), y=y, M=M, Minv=Minv, m0=m0,k=k,weight=weight,
                         noconst_last=noconst_last, y_nonneg=y_nonneg,jackup=jackup,
                         plim=plim, bc_value = bc_value,y_added2=y_added2)
        b     <-res2$b
        z     <-res2$z
        loglik<-res2$loglik
        comp  <-res2$comp
        y_ms  <-res2$y_ms
        z_ms  <-res2$z_ms
        #if(y_nonneg){
        #  est <- est[1:(length(est)-2)]
        #}

      } else if(y_nonneg){
        #res   <-optimize(f=NLL_bc, interval=c(-5,5),y=y,M=M,Minv=Minv,m0=m0 )#c(-5,5)
        #bc_par<-res$minimum
        #est   <-NULL

        if(!is.null(bc_value)){
          bc_par0<-c(bc_value, 0.001)
        } else {
          bc_par0<-c(1, 0.001)
        }
        res   <-optim(par=bc_par0,fn=NLL_bc,y=y,M=M,Minv=Minv,m0=m0,jackup=jackup,
                      weight=weight, bc_value=bc_value,y_added2=y_added2)# interval=c(-5,5)
        bc_par<-res$par
        bc_par[2]<-abs(bc_par[2])
        if( !is.null( bc_value ) ){
          bc_par[1]<- ifelse( bc_value < -1, -1, ifelse( bc_value > 5, 5, bc_value) )
        } else {
          bc_par[1]<- ifelse( bc_par[1] < -5, -5, ifelse( bc_par[1] > 5, 5, bc_par[1]) )
        }

        est   <-NULL

        res2  <-NLL_bc2(par=bc_par, y=y, M=M, Minv=Minv, m0=m0,jackup=jackup,weight=weight,
                        bc_value=bc_value,y_added2=y_added2 )
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
        M0_sub	<- as.matrix(M0[ id != par0_sel, id != par0_sel ])
        term2		<- determinant( M0_sub )$modulus
        Msub_00	<- M[ id == par0_sel, id == par0_sel ]
        Msub_01	<- M[ id == par0_sel, id != par0_sel ]
        if(sum(id == par0_sel)==1){
          term3_0	<- Msub_00 - c( t( Msub_01 ) %*% solve( M0_sub, tol = 1e-30 ) %*% Msub_01 )
        } else {
          term3_0	<- Msub_00 - Msub_01 %*% solve( M0_sub, tol = 1e-30 ) %*% t( Msub_01 )
        }
      }
      return(list(term2 = term2, term3_0 = term3_0))
    }

    Mdet_f	  	<- function( evSqrt, id, term2, term3_0, par0_sel ){
      term1		  <- sum( log( evSqrt ) ) * 2
      if( dim( as.matrix( term3_0 ) )[1] == 1 ){
        term3_0 <- term3_0 + 1/evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
        term3_0 <- as.matrix(term3_0)
      } else {
        diag( term3_0 ) <- diag( term3_0 ) + 1/evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
      }

      term3	  <- determinant( term3_0 )$modulus
      if( is.null( term2 )  ){
        Mdet    <- term1 + term3
      } else {
        Mdet	  <- term1 + term2 + term3
      }
      return(Mdet)
    }

  }

  lik_resf_vc		<- function( par0, par0_est, par0_id, par0_sel, ev, M, M0inv, M0inv_01, M0inv_00,
                            m, yy, b_01, b_02, n, nx, nsv, nnxf, nnsv, ng, emet, term2, term3_0, null_dum2,
                            id, tr_comp=0 ){
    par		    <- par0 ^ 2
    par_est		<- par0_est ^ 2
    par[ par0_id == par0_sel ]  <- par_est
    evSqrt	  <- NULL
    for( i in ( 1:nsv )[ null_dum2[ 1:nsv ] == 0 ] ){
      evv  	  <- ev ^ par[ nsv + i ] * sum( ev ) / sum( ev ^ par[ nsv + i ] )
      evSqrt	<- c( evSqrt, par[ i ] * sqrt( evv ) )
    }

    if( ng != 0 ){
      for( j in (1:ng)[ null_dum2[ (( nsv + 1 ):( nsv + ng )) ] == 0 ] ){
        xgg	  <- rep( 1, sum( id == nsv + j ) )
        evSqrt<- c( evSqrt, par[ 2 * nsv + j ] * xgg )
      }
    }

    if( nnxf != 0 ){
      for( i2 in (1:nnxf)[ null_dum2[ (nsv+ng+1):(nsv+ng+nnxf) ] == 0 ] ){
        xgg	  <- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
        evSqrt<- c( evSqrt, par[ 2 * nsv + ng + i2 ] * xgg )
      }
    }

    if( nnsv != 0 ){
      for( i2 in (1:nnsv)[ null_dum2[ (nsv+ng+nnxf+1):(nsv+ng+nnxf+nnsv) ] == 0 ] ){
        xgg	  <- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
        evSqrt<- c( evSqrt, par[ 2 * nsv + ng + nnxf + i2 ] * xgg )
      }
    }

    Mdet		<- Mdet_f( id = id, par0_sel=par0_sel, term2 = term2, term3_0 = term3_0, evSqrt = evSqrt )
    M2		  <- M
    for( j in 1:max( id ) ){
      if( sum( id == j )==1){
        M[ id == j, id == j ]       <- M[ id == j, id == j ] + 1/evSqrt[ id[ -c( 1:nx ) ] == j ] ^ 2
      } else {
        diag( M[ id == j, id == j ])<- diag( M[ id == j, id == j ] ) + 1/evSqrt[ id[ -c( 1:nx ) ] == j ] ^ 2
      }
    }

    if( dim(as.matrix(M0inv_00))[1]==1){
      M0inv_00	      <- M0inv_00 + evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
    } else {
      diag(M0inv_00)	<- diag( M0inv_00 ) + evSqrt[ id[ id != 0 ] == par0_sel ] ^ 2
    }

    b_02_b		<- solve( M0inv_00, tol = 1e-30 ) %*% b_02
    b_02		<- M0inv_01 %*% b_02_b
    b		<- b_01 - b_02
    sse		<- yy - 2 * t( b ) %*% m + t( b ) %*% M2 %*% b
    dd		<- abs(sse) + sum( ( b[ -( 1:nx ) ] / evSqrt ) ^ 2 )
    if( emet == "reml" ){
      loglik	<- Mdet + ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) ) - 2*tr_comp
    } else if( emet == "ml" ){
      loglik	<- Mdet + n * ( 1 + log( 2 * pi * dd / n ) ) - 2*tr_comp
    }
    return( loglik )
  }

  lik_resf_vc0	<- function( par0, ev, M, m, yy, n, nx, nsv, ng, nnxf, nnsv, emet, null_dum4, tr_comp=0 ){
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
    if( inherits(test, "try-error") ){#class(test)[1]=="try-error"
      loglik  	<- Inf
    } else {
      b		<- Minv %*% m
      sse		<- yy - 2 * t( b ) %*% m + t( b ) %*% M0 %*% b
      dd		<- abs(sse) + sum( b[ -( 1:nx ) ] ^ 2 )
      if( emet == "reml" ){
        term1	<- determinant( M )$modulus
        term2	<- ( n - nx ) * ( 1 + log( 2 * pi * dd / ( n - nx ) ) )
      } else if( emet == "ml" ){
        term1	<- determinant( as.matrix( M[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
        term2	<- n * ( 1 + log( 2 * pi * dd / n ) )
      }
      loglik		<- term1 + term2 -2*tr_comp
    }
    return( loglik[ 1 ] )
  }

  sel_basis_vif <-function(test0, vif_max=15){
    tres3<-TRUE
    sel  <-1:dim(test0)[2]
    while(tres3){
      suppressWarnings(cor_test0<-cor(test0))
      testt2<-try(vif <- diag(solve(cor_test0)), silent=TRUE)
      if( inherits(testt2, "try-error" ) ){#class(testt2)=="try-error"
        suppressWarnings(corr_0<-cor(test0))
        diag(corr_0)<-0
        cor_max <-apply(corr_0,2,max)
        test0   <-test0[, -(which(cor_max == max(cor_max))[1])]
        sel     <-sel[-(which(cor_max == max(cor_max))[1])]

      } else {
        if( max(vif)>vif_max ){
          test0<-test0[,-(which(vif==max(vif))[1])]
          sel     <-sel[-(which(vif==max(vif))[1])]
        } else {
          tres3<-FALSE
        }
      }
    }
    return(list(basis=test0,sel=sel))
  }


  n     	   <- length( y )
  resf_flag  <- ifelse( is.null(x), 1, 0)
  if( !is.null(nongauss)){
    y_type   <- nongauss$y_type
    tr_num   <- nongauss$tr_num
    y_nonneg <- nongauss$y_nonneg
    if( y_nonneg ){
      if( min( y ) < 0 ) stop(" y must be non-negative when y_nonneg = TRUE")
    }

  } else {
    y_type   <- "continuous"
    tr_num   <- 0
    y_nonneg <- FALSE
  }

  y_added <- 0
  y_added2<- 0
  if(y_type  == "count"){
    bc_value <- NULL
    y_nonneg <- FALSE
    y_org    <- y
    if( is.null(offset) ){
      y      <- log( y_org + 0.5 ) + y_added2
    } else {
      if( min( offset ) <= 0 ) stop( "offset must be positive" )
      y      <- log( ( y_org+0.5 )/offset ) + y_added2
    }

    zrat     <- sum(y_org==0)/length(y_org)
    y_added  <- 0#0.5
    y_added2 <- 0#-(1+0.5*zrat)/(y_org+0.5)
    #if(tr_num==0) y<- y + y_added2
    #y<- y + y_added2
    #y_added2<-0

    if( is.null(weight) ){
      weight0<- 1
      weight <- y_org + 0.5
    } else {
      weight0<- weight
      weight <- weight0*(y_org+0.5)
    }

  } else if( y_type =="continuous"){
    if( !is.null(offset) ) {
      message("Note: offset is available only for count data (y_type = count). It is ignored")
      offset <- NULL
    }

    if(min(y) <= 0.00000001){
      if(y_nonneg){
        y_added <- 0.00000001
        y       <- y + y_added
        #jackup  <- TRUE
      }
    }

  } else {
    stop( "y_type must be continuous or count" )
  }

  jackup  <- FALSE
  noconst_last <- TRUE

  if( !is.null(weight) ){
    if( min(weight) <= 0) stop( "weight must be positive" )

    if(length(weight)==1){
      w_scale <- 1/weight
    } else {
      w_scale <- n/sum(weight)
    }

    weight  <-c( weight*w_scale )
  } else {
    w_scale <- 1
  }

  nvc_x       <- x_nvc
  nvc_xconst  <- xconst_nvc

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
    allnvc_xconst     <- !xconst_nvc_sel
    allnvc_xconst_flag<-0
  } else {
    allnvc_xconst<- xconst_nvc_sel
    allnvc_xconst_flag<-1
  }

  if( x_nvc[1]      == FALSE ) x_nvc_sel = FALSE
  if( xconst_nvc[1] == FALSE ) xconst_nvc_sel = FALSE

  if( method == "reml" ){
    lik_nam	<- "rlogLik"
  } else if( method == "ml" ){
    lik_nam	<- "logLik"
  }

  if( n > 150000 ){
    message( paste( "Note: besf_vc function is available for large samples. see help(besf_vc)" ) )
  }

  if(is.null(xconst)) nvc_xconst <- FALSE

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
      X1	  <- NULL
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

  XX1_0	      <- list(NULL)
  XX1	        <- NULL
  if( is.logical( nvc_x[ 1 ] )&( nvc_x[ 1 ] == TRUE) ) nvc_x     <- 1:nx0

  sel_basis_n <-list( NULL )
  if( is.logical( nvc_x[ 1 ] ) == FALSE ){
    X1_nvc  <- as.matrix(X1[ , nvc_x ])
    xxname	<- names( as.data.frame( X1 ) )[ nvc_x ]
    nnsv    <- length( xxname )
    np_xx     <-apply( X1_nvc,2,function( x ) length( unique( x )))
    np_xx     <-ifelse( np_xx < nvc_num/0.7, round( np_xx * 0.7 ) ,nvc_num )
    np_xx_max <-round( n/nnsv ) - 2
    np_xx[ np_xx > np_xx_max ] <-np_xx_max
    np_xx[ np_xx < 2 ] <- 2

    for( ii in 1:dim( X1_nvc )[ 2 ] ){
      test<-TRUE
      iiii<-0
      while(test){
        kkk  <- np_xx[ ii ]-iiii
        knots<-seq(min( X1_nvc[ , ii ] ),max( X1_nvc[ , ii ] ),len=kkk+2)[2:(kkk+1)]
        testt<- try(XX1_00<- ns( X1_nvc[ ,ii ], knots = knots ), silent=TRUE)#replaced
        test <- inherits(testt, "try-error")#class(testt)[1] == "try-error"
        iiii <- iiii+1
      }

      XX1_00      <- cbind( X1_nvc[,ii], XX1_00)
      sel_basis_n[[ ii ]]    <- 1:dim(XX1_00)[2]
      testt2<-try(XX1_00_mod <- sel_basis_vif( XX1_00, vif_max=15),silent=TRUE)
      if( !inherits(testt2, "try-error") ){#class(testt2)[1] != "try-error"
        XX1_00               <- XX1_00_mod$basis
        sel_basis_n[[ ii ]]  <- XX1_00_mod$sel
        np_xx[ ii ]          <- dim(XX1_00)[2]
      } else {
        np_xx[ ii ]          <- 0
      }
      if( np_xx[ ii ] > 2 ){
        XX1_0[[ii]]<- scale( XX1_00 )
        XX1	      <- cbind( XX1, X1_nvc[, ii ] * XX1_0[[ii]] )
      } else {
        XX1_0[[ii]]<- 0
        np_xx[ ii ]<- 0
      }
    }
  } else {
    nnsv      <- 0
    np_xx     <- NULL
  }

  XXconst_0	  <- list( NULL )
  XXconst	    <- NULL
  if( is.logical( nvc_xconst[ 1 ] )&( nvc_xconst[ 1 ] == TRUE) ) nvc_xconst<- 1:nxf

  sel_basis_c <-list( NULL )
  if( is.logical( nvc_xconst[ 1 ] ) ==FALSE ){
    Xconst_nvc<- as.matrix( Xconst[ , nvc_xconst ] )
    xxfname	  <- names( as.data.frame( Xconst ) )[ nvc_xconst ]
    nnxf      <- length( xxfname )
    np_xxconst       <-apply( Xconst_nvc,2,function( x ) length( unique( x )))
    np_xxconst       <-ifelse( np_xxconst < nvc_num/0.7, round( np_xxconst * 0.7 ) ,nvc_num )
    np_xxconst_max <-round( n/nnxf ) - 2
    np_xxconst[ np_xxconst > np_xxconst_max ] <-np_xxconst_max
    np_xxconst[ np_xxconst < 2 ] <- 2

    for( ii in 1:dim( Xconst_nvc )[ 2 ] ){
      test<-TRUE
      iiii<-0
      while(test){
        kkk      <- np_xxconst[ ii ]-iiii
        knots    <-seq(min( Xconst_nvc[ ,ii ] ),max( Xconst_nvc[ ,ii ] ),len=kkk+2)[2:(kkk+1)]
        #testt<-try(XXconst_00<- scale( poly( Xconst_nvc[ , ii ], np_xxconst[ ii ]-iiii )[ , -1 ] ), silent=TRUE)
        testt<- try(XXconst_00<- ns( Xconst_nvc[ , ii], knots = knots ), silent=TRUE)
        test <- inherits(testt, "try-error")#class(testt)[1] == "try-error"
        iiii <- iiii+1
      }

      XXconst_00<- cbind( Xconst_nvc[,ii], XXconst_00)
      sel_basis_c[[ ii ]]    <- 1:dim(XXconst_00)[2]
      testt2<-try(XXconst_00_mod     <- sel_basis_vif( XXconst_00, vif_max=15),silent=TRUE)
      if( !inherits(testt2, "try-error") ){#class(testt2)[1] != "try-error"
        XXconst_00       <- XXconst_00_mod$basis
        sel_basis_c[[ii]]<- XXconst_00_mod$sel
        np_xxconst[ ii ] <- dim(XXconst_00)[2]
      } else {
        np_xxconst[ ii ]<- 0
      }
      if( np_xxconst[ ii ] > 2 ){
        XXconst_0[[ii]]<- scale( XXconst_00 )
        XXconst	       <- cbind( XXconst, Xconst_nvc[, ii ] * XXconst_0[[ii]] )
        #np_xx[ ii ]<- np_xx[ ii ]  - iiii
      } else {
        XXconst_0[[ii]]<- 0
        np_xxconst[ ii ]<- 0
      }
    }
  } else {
    nnxf         <- 0
    np_xxconst   <- NULL # added
  }

  if( ( nx0 > 0 ) & ( nxf > 0 ) ){
    dup      <-duplicated(cbind(X1,Xconst),MARGIN=2)
    if(sum(dup)>0){
      stop( "x and xconst cannot have the same column" )
    }
  }

  nsv	  <- ifelse( is.null( X1 ),1, dim( X1 )[ 2 ] + 1 )
  nev0	<- min( round( n /nsv ) - 2, length( meig$ev ))
  nev0  <- max(nev0, 3)
  meig$sf	<- meig$sf[ , 1 : nev0 ]
  meig$ev	<- meig$ev[   1 : nev0 ]

  if( !is.null( xgroup ) ){
    xgroup  <- data.frame(xgroup)
    ng	    <- dim( xgroup )[ 2 ]

    if( sum(is.na(xgroup)) > 0) stop("NA is not allowed in xgroup")

    if( ng == 1 ){
      xg_num<- length( unique( xgroup[,1] ) ) - ng
    } else {
      xg_num<- sum( apply(xgroup,2,function(x) length(unique(x)))) - ng
    }

    nev0g   <- round( ( n - xg_num - nxf )/nsv ) - 2
    if( nev0 > nev0g ){
      if( nev0g < 2 ){
        stop("Error: Too many groups. Reduce variables in xgroup")
      } else {
        message( paste( "Note: ", nev0 - nev0g, " eigenvectors are omitted to stablize the estimates",sep=""))
        nev0    <- nev0g
        meig$sf	<- meig$sf[ , 1 : nev0 ]
        meig$ev	<- meig$ev[   1 : nev0 ]
      }
    }
  }

  X2	    <- meig$sf
  ev	    <- meig$ev
  meig0   <- X2
  id	    <- c( rep( 0, nsv ), rep( 0, nxf ), rep( 1, length( ev ) ) )
  if( nsv >= 2 ){
    for( i in 1:( nsv - 1 ) ){
      X2  <- cbind( X2, X1[, i ] * meig0 )
      id	<- c( id, rep( i+1, length( ev )))
    }
  }

  if( is.null( xgroup ) == FALSE ){
    xg_id0	<- nsv + 1
    xg_idid2<- NULL
    xg_link_id <- list( NULL )
    xg_levels  <- list(NULL)
    for( ff in 1:ng){
      xg00	<- factor( xgroup[ , ff ] )
      xg_levels[[ff]] <- levels( xg00 )

      xg0	  <- sparse.model.matrix( ~ 0 + xg00 )
      xg0s	<- Matrix::colSums( xg0 )
      xg_idid	<- max( which( xg0s == max( xg0s ) ) )
      Xg0	  <- as.matrix( xg0[ , -xg_idid ] )
      xg_id1	<- rep( xg_id0, dim( Xg0 )[2] )
      if( ff == 1 ){
        Xg	<- Xg0
        xg_id	<- xg_id1
      } else {
        Xg	<- cbind( Xg, Xg0 )
        xg_id	<- c( xg_id, xg_id1 )
      }
      xg_id0	<- xg_id0 + 1
      xg_idid2<-c(xg_idid2, xg_idid)

      xgroup_datt <- data.frame(id=1:length(xgroup[,ff]),xgroup_sub=xgroup[,ff])
      xg00_datt   <- data.frame(id_b_g = 1:length(levels(xg00)), xgroup_sub=levels(xg00))
      xgroup_datt2<- merge(xgroup_datt,xg00_datt,by="xgroup_sub",all.x=TRUE)
      xg_link_id[[ ff ]]  <- xgroup_datt2[order(xgroup_datt2$id),"id_b_g"]
    }
    X2	         <- cbind( X2, Xg )
    id	         <- c(id, xg_id )

  } else {
    ng	<- 0
    xg_levels <- NULL
  }

  pp_id    <-max( id )
  if( is.logical( nvc_xconst[ 1 ] ) == FALSE ){
    xxc_id   <-NULL
    for(pp in 1:length(np_xxconst)){
      xxc_id <- c(xxc_id, rep(pp_id + pp, np_xxconst[pp]))
    }
    id	   <- c( id, xxc_id )
    X2     <- cbind( X2, XXconst )
    pp_id  <- pp_id + length(np_xxconst)
  }

  if( is.logical( nvc_x[ 1 ] ) == FALSE ){
    xx_id   <-NULL
    for(pp in 1:length(np_xx)){
      xx_id <- c(xx_id, rep(pp_id + pp, np_xx[pp]))
    }
    id	   <- c( id, xx_id )
    X2     <- cbind( X2, XX1 )
  }

  if( is.null( Xconst )&is.null( X1 )){
    X0	<- as.matrix( rep(1,n) )
  } else {
    X0	<- as.matrix( cbind( 1, Xconst, X1 ) )
  }
  X	  <- as.matrix( cbind( X0, X2 ) )

  if( !is.null(weight) ) {
    X0_org<- X0
    X2_org<- X2
    X_org <- X
    X0    <- sqrt(weight) * X0_org
    X2    <- sqrt(weight) * X2_org
    X     <- sqrt(weight) * X_org
  } else {
    X0_org<- NULL
    X2_org<- NULL
    X_org <- NULL
  }

  Mo 	<- crossprod( X0 )
  nx	<- dim( X )[ 2 ] - dim( X2 )[ 2 ]
  M   <- crossprod( X )

  y0  <- y
  plim<- 10
  if( ( tr_num > 0 )|( y_nonneg == TRUE ) ){
    Moinv  <- solve(Mo)
    if( !is.null(weight) ){
      m0   <- weight * X0_org
    } else {
      m0   <- X0
    }
    cv_init<-lm_cw(y=y0, M=Mo, Minv=Moinv,m0=m0,k=tr_num,noconst_last=noconst_last,
                   y_nonneg=y_nonneg,jackup=jackup, weight = weight, plim=plim, bc_value=bc_value,
                   y_added2=y_added2)
    y      <-cv_init$z
    tr_par <-cv_init$vpar[1:tr_num]
    tr_bpar<-cv_init$bc_par
    if( y_nonneg == TRUE ){
      if( jackup ){
        tr_bpar[2]<-abs(tr_bpar[2])
      } else {
        tr_bpar[2]<-0
      }
    }

    tr_npar<-length(unlist(c(tr_par,tr_bpar)))
    #if( jackup==FALSE ) tr_npar<- tr_npar - 1
    tr_comp<-cv_init$comp
    tr_par0<-NULL

  } else {
    tr_par <- NULL
    tr_bpar<- NULL
    tr_comp<- 0
    tr_npar<- 0
  }

  if( !is.null(weight) ) {
    w_y <- weight * y
    mo	<- crossprod( X0_org, w_y )
    m	  <- crossprod( X_org , w_y )
    yy  <- sum( y * w_y )
    parVmax_sq <- sqrt( sd( y ) / sd( y - X0_org %*% ( solve( Mo ) %*% mo ) ) )

  } else {
    mo	<- crossprod( X0, y )
    m	  <- crossprod( X, y )
    yy  <- sum( y ^ 2 )
    parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )
  }


  if( penalty == "aic" ){
    pen	<- 2
  } else if( penalty == "bic" ){
    pen	<- log( n )
  }

  par0	  <- rep( 1, 2 * nsv )
  par0[1] <- parVmax_sq / 3
  par0_est<- c( 1, 1 )
  par0_id	<- rep( 1:nsv, 2 )
  null_dum<- rep( 0, nsv )

  if( is.null( xgroup ) == FALSE ){
    par0	<- c( par0, rep( parVmax_sq/3, ng ) )
    par0_id	<- c( par0_id , max(par0_id) + ( 1 : ng ) )
    null_dum<- c( null_dum, rep( 0, ng ) )
  }
  if( is.logical( nvc_xconst[ 1 ] ) == FALSE ){
    par0	  <- c( par0, rep( 1, nnxf ) )
    par0_id	<- c( par0_id , max(par0_id) + ( 1:nnxf ) )
    null_dum0<- rep( 0, nnxf )
    null_dum0[np_xxconst==0]<-1
    null_dum<- c( null_dum, null_dum0 )
  }
  if( is.logical( nvc_x[ 1 ] ) == FALSE ){
    par0	<- c( par0, rep( 1, nnsv ) )
    par0_id	<- c( par0_id , max(par0_id) + ( 1:nnsv ) )
    null_dum0<- rep( 0, nnsv )
    null_dum0[np_xx==0]<-1
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
    for( par0_sel in (1:( nsv+ ng + nnxf + nnsv ))[id_exist]){
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
            if( sum(id_sub)==1 ){
              M0[ id_sub, id_sub ] <- M0[ id_sub, id_sub ] + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
            } else {
              diag( M0[ id_sub, id_sub ] )<-
                diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
            }
          } else {
            id_sub	<- ( id == ( j2 + nsv ) )
            g_add   <- (1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2)/100

            if( sum(id_sub)==1 ){
              M0[ id_sub, id_sub ]        <- M0[ id_sub, id_sub ] + g_add
            } else {
              diag( M0[ id_sub, id_sub ] )<- diag( M0[ id_sub, id_sub ] ) + g_add
            }
          }
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
      id_omit2_0<- which( id_omit1 != 0 )
      id_omit2<-id_omit2_0[2:(nsv + 1)]
      id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0

      sstop <-FALSE
      if( ( ng != 0 )&( nsv < par0_sel )&( par0_sel <= nsv + ng ) ){
        if( rcond( MM0 ) < 1e-30 ){
          par0[ par0_id == par0_sel ]<- 0
          null_dum[ par0_sel ]		   <- 1
          loglik 	<- ( -1 / 2 ) * c( res_ref )
          res_old	<- res_ref
          score   <- res_ref+ pen * (nx+sum(par0!=0)+1 + tr_npar)
          if( !is.null(weight)) score <- score - sum( log( weight ))

          message( paste( "Note:", "group effect", par0_sel - nsv, "is omitted to stablize the estimates", sep = " " ) )
          gmess   <- gmess + 1
          warn    <- TRUE
          sstop   <- TRUE
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
            M	<- M [ id_omit1 == 0, id_omit1 == 0 ]
            M0	<- M0[ id_omit1 == 0, id_omit1 == 0 ]
            m	<- m [ id_omit1 == 0 ]
            id	<- id[ id_omit1 == 0 ]
            X	<- X [, id_omit1 == 0 ]
            if( !is.null( weight ) ){
              X_org	<- X_org[, id_omit1 == 0 ]
            }

            null_dum3<- null_dum3[ id_omit1 == 0 ]

            MM	<- M [ null_dum3, null_dum3 ]
            MM0	<- M0[ null_dum3, null_dum3 ]
            mm	<- m [ null_dum3 ]
            idd	<- id[ null_dum3 ]

            id_omit1<- c( diff( id ), 1)
            id_omit2_0<- which( id_omit1 != 0 )
            id_omit2<-id_omit2_0[2:(nsv + 1)]
            id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0
            n_omit0	<- n_omit0 + 1
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

        ev	    <- ev     [   1:sum( id == 1 ) ]
        meig$sf	<- meig$sf[ , 1:sum( id == 1 ) ]
        term2	  <- Mdet0$term2
        term3_0	<- Mdet0$term3_0

        if( min( par0[ par0_id == par0_sel ] ) >= 1e-5 ){
          par00	<- par0[ par0_id == par0_sel ]
        } else {
          if(par0_sel<=nsv){
            par00_id<- max( which( c(apply(Par[,par0_id == par0_sel], 1, min)) != 0 ) )
          } else {
            par00_id<- max( which( Par[,par0_id == par0_sel] != 0 ) )
          }

          par00	<- Par[ par00_id, par0_id == par0_sel ]
        }

        if( n_omit0 > 0 ){
          res_old <- lik_resf_vc0( par0, ev = ev, M = MM, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                                   n = n, nx = nx, nsv = nsv, ng = ng, emet = method, null_dum4 = null_dum2, tr_comp=tr_comp )
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
                            term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                            n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                            par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2,
                            lower = llim, upper = ulim, method = omethod, tr_comp=tr_comp )
          res_int <- res

        } else if( par0_sel <= nsv ){
          res    	<- optim( fn = lik_resf_vc, par00, par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                            M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                            term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                            n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                            par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2, tr_comp=tr_comp )

        } else if( par0_sel <= nsv+ng ){
          res    	<- optimize( f = lik_resf_vc, lower = parVmax_sq/1000, upper = parVmax_sq, par00,
                               par0 = par0, ev = ev, M = MM, M0inv = M0inv,
                               M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                               term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                               n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                               par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2, tr_comp=tr_comp )
          res$value	<- c( res$objective )
          res$par	<- c( res$minimum )

        } else {
          res    	<- optimize( f = lik_resf_vc, lower = parVmax_sq/10^5, upper = 10^10, par00, par0 = par0, ev = ev,
                               M = MM, M0inv = M0inv,
                               M0inv_01 = M0inv_01, M0inv_00 = M0inv_00, b_01 = b_01, b_02 = b_02,
                               term2 = term2, term3_0 = term3_0, m = mm, yy = yy,nnxf=nnxf,nnsv=nnsv,
                               n = n, nx = nx, nsv = nsv, ng = ng, emet = method, id = idd,
                               par0_sel = par0_sel, par0_id = par0_id, null_dum2 = null_dum2, tr_comp=tr_comp )
          res$value	<- c( res$objective )
          res$par	<- c( res$minimum )
        }

        if( ( iter > 3 ) & ( res$value > res_old ) ){
          loglik	<- ( - 1 / 2 ) * res_old

        } else {
          if( par0_sel != 1 ){
            MM_ref 	<- MM [ idd != par0_sel, idd != par0_sel ]
            mm_ref	<- mm [ idd != par0_sel ]
            null_dum4<- null_dum2
            null_dum4[ par0_sel ] <- 1
            res_ref <- lik_resf_vc0( par0, ev = ev, M = MM_ref, m = mm_ref, yy = yy,nnxf=nnxf,nnsv=nnsv,
                                     n = n, nx = nx, nsv = nsv, ng = ng, emet = method, null_dum4 = null_dum4, tr_comp=tr_comp )
          } else {
            res_ref	<- Inf
          }

          np_add       <- length(par00)
          flag         <- 0
          if( par0_sel <= nsv ){
            if( allsvc_flag == 0 ){
              allvc  <- allsvc
            } else {
              flag   <-1
              allvc  <- c(TRUE, allsvc_id)[ par0_sel ]
            }
          } else if( ( nsv < par0_sel )&( par0_sel <= nsv + ng )){
            allvc  <- TRUE
          } else if( ( nsv +ng < par0_sel )&( par0_sel <= nsv + ng + nnxf )){
            if( allnvc_xconst_flag == 0 ){
              allvc  <- allnvc_xconst
            } else {
              flag   <-1
              allvc  <- allnvc_xconst_id[ par0_sel - nsv - ng ]
            }
          } else {
            if( allnvc_x_flag == 0 ){
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
            score   <- res_ref+ pen * (nx+sum(par0!=0)+1+tr_npar)
            if( !is.null( weight ) ) score <- score - sum( log( weight ))

          } else {      ####### accept
            par0[ par0_id == par0_sel ]	<- res$par
            null_dum[ par0_sel ]		<- 0
            res_old	<- res$value
            loglik 	<- ( -1 / 2 ) * res$value
            score   <- res_old + pen * (nx+sum(par0!=0)+1+tr_npar)
            if( !is.null( weight ) ) score <- score - sum( log( weight ))
          }
        }
      }

      LL0	<- c( LL0, loglik )
      if( !is.null(x) ) print( paste( par0_sel, "/", nsv + ng + nnxf + nnsv, sep = "" ) )
      #print(score)
    }

    ########################################### TWGP: main
    if((tr_num != 0)|( y_nonneg )){ ####### with transformation
      evSqrt	<- NULL
      par0_sq	<- par0 ^ 2
      for( i in ( 1:nsv )[ null_dum[ 1:nsv ] == 0 ] ){
        evv  	  <- ev ^ par0_sq[ nsv + i ] * sum( ev ) / sum( ev ^ par0_sq[ nsv + i ] )
        evSqrt	<- c( evSqrt, par0_sq[ i ] * sqrt( evv ) )
      }

      M0	<- M
      for( j in ( 1:nsv )[ null_dum[ 1:nsv ] == 0 ] ){
        id_sub	<- ( id == j )
        diag( M0[ id_sub, id_sub ] )<-
          diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == j ] ^ 2
      }

      if( ng != 0 ){
        for( i2 in (1:ng)[ null_dum[ (( nsv + 1 ):( nsv + ng )) ] == 0 ] ){
          xgg	<- rep( 1, sum( id == ( nsv + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ 2 * nsv + i2 ] * xgg )
        }
        for( j2 in (1:ng)[ null_dum[ (( nsv + 1 ):( nsv + ng )) ] == 0 ] ){
          id_sub	<- ( id == ( j2 + nsv ) )
          if( sum(id_sub)==1 ){
            M0[ id_sub, id_sub ] <- M0[ id_sub, id_sub ] + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
          } else {
            diag( M0[ id_sub, id_sub ] )<-
              diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
          }
        }
      }

      if( nnxf != 0 ){
        for( i2 in (1:nnxf)[ null_dum[ (nsv+ng+1):(nsv+ng+nnxf) ] == 0 ] ){
          xgg	<- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng ) + i2 ] * xgg )
        }

        for( j2 in (1:nnxf)[ null_dum[ (nsv+ng+1):(nsv+ng+nnxf) ] == 0 ] ){
          id_sub	<- ( id == ( j2 + nsv + ng ) )
          diag( M0[ id_sub, id_sub ] )<-
            diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng) ] ^ 2
        }
      }

      if( nnsv != 0 ){
        for( i2 in (1:nnsv)[ null_dum[ (nsv+ng+nnxf+1):(nsv+ng+nnxf+nnsv) ] == 0 ] ){
          xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
          evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng + nnxf ) + i2 ] * xgg )
        }

        for( j2 in (1:nnsv)[ null_dum[ (nsv+ng+nnxf+1):(nsv+ng+nnxf+nnsv) ] == 0 ] ){
          id_sub	<- ( id == ( j2 + nsv + ng + nnxf) )
          diag( M0[ id_sub, id_sub ] )<-
            diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng + nnxf) ] ^ 2
        }
      }

      null_dum3<- c( rep( 0, nx ), null_dum[ id ]) == 0
      MM	<- M [ null_dum3, null_dum3 ]
      MM0	<- M0[ null_dum3, null_dum3 ]
      mm	<- m [ null_dum3 ]
      idd	<- id[ null_dum3 ]
      id_omit1<- c( diff( id ), 1)
      id_omit2_0<- which( id_omit1 != 0 )                           ### changed
      id_omit2<-id_omit2_0[2:(nsv + 1)]                             ### changed
      id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0### changed

      er_dum	<- TRUE
      n_omit0	<- 0
      while( er_dum == TRUE ){
        try1	<- try( M0inv	<- solve( MM0,  tol = 1e-30 ), silent = TRUE )
        er_dum<- inherits( try1, "try-error")#( class(try1)[1] =="try-error" )

        if( er_dum == TRUE ){
          M	<- M [ id_omit1 == 0, id_omit1 == 0 ]
          M0	<- M0[ id_omit1 == 0, id_omit1 == 0 ]
          m	  <- m [ id_omit1 == 0 ]
          id	<- id[ id_omit1 == 0 ]
          X	  <- X [, id_omit1 == 0 ]
          if( !is.null( weight ) ) {
            X_org	 <- X_org[, id_omit1 == 0 ]
          }

          null_dum3<- null_dum3[ id_omit1 == 0 ]
          evSqrt   <- evSqrt[ id_omit1[ -( 1:nx ) ] == 0 ]

          MM	<- M [ null_dum3, null_dum3 ]
          MM0	<- M0[ null_dum3, null_dum3 ]
          mm	<- m [ null_dum3 ]
          idd	<- id[ null_dum3 ]

          id_omit1<- c( diff( id ), 1)
          id_omit2_0<- which( id_omit1 != 0 )                           ### changed
          id_omit2<-id_omit2_0[2:(nsv + 1)]                             ### changed
          id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0### changed
          n_omit0	<- n_omit0 + 1
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

      ev	    <- ev[   1:sum( id == 1 ) ]
      meig$sf	<- meig$sf[ , 1:sum( id == 1 ) ]

      MM[ -( 1:nx ), -( 1:nx ) ]	<- t(MM[ -( 1:nx ), -( 1:nx ) ] * evSqrt ) * evSqrt
      MM[ -( 1:nx ),    1:nx   ]	<-   MM[ -( 1:nx ),    1:nx   ] * evSqrt
      MM[    1:nx  , -( 1:nx ) ]	<- t(MM[ -( 1:nx ),    1:nx   ] )
      MM0		<- MM
      diag( MM[ -( 1:nx ), -( 1:nx ) ] ) <- diag( MM[ -( 1:nx ), -( 1:nx ) ] ) + 1
      MMinv	<- solve( MM, tol = 1e-30 )
      if( method == "reml" ){
        term1	<- determinant( MM )$modulus
      } else if( method == "ml" ){
        term1	<- determinant( as.matrix( MM[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
      }
    }

    if(tr_num != 0){

      if(is.null(tr_par0)){
        tr_par0  <-rep(c(1,0,0,1),tr_num)
        if(noconst_last) tr_par0<-tr_par0[1:(4*tr_num-2)]
        if(y_nonneg)    tr_par0<-c(tr_par0,1,0.01)
      }

      res      <-optim(par=tr_par0,fn=lik_resf_vc_tr, k=tr_num, y_nonneg=y_nonneg,noconst_last=noconst_last,
                       evSqrt=evSqrt, y0=y0, n=n, nx=nx, emet=method, X=X,M0=MM0,Minv=MMinv,term1=term1,
                       null_dum3=null_dum3,jackup=jackup, plim=plim, y_type=y_type,
                       y_added2=y_added2,weight=weight, control=list(maxit=100))#,
      #method="L-BFGS-B",lower=lower,upper=upper)

      score_tr   <-res$value + pen * (nx+sum(par0!=0)+1+tr_npar)
      if( !is.null( weight )) score_tr <-score_tr - sum( log( weight ))

      if(score_tr< score){
        res_old	 <- res$value
        loglik	 <- (-1/2)*c( res$value )
        score    <-score_tr
        tr_est0  <-res$par#c(1,0,0,1)
        tr_est0[tr_est0 >  plim] <-  plim
        tr_est0[tr_est0 < -plim] <- -plim

        tr_par   <-list(NULL)
        for(kk in 1:tr_num){
          if(noconst_last & (kk==tr_num)){
            tr_est0[(4*(kk-1)+1)]<- abs( tr_est0[(4*(kk-1)+1)] )
            tr_par[[kk]]         <- tr_est0[(4*(kk-1)+1):(4*(kk-1)+2)]
          } else {
            tr_est0[(4*(kk-1)+1)]<- abs( tr_est0[(4*(kk-1)+1)] )
            tr_est0[(4*(kk-1)+4)]<- abs( tr_est0[(4*(kk-1)+4)] )
            tr_par[[kk]]         <- tr_est0[(4*(kk-1)+1):(4*(kk-1)+4)]
          }
        }

        if(y_nonneg){
          np_b      <-length(tr_est0)
          tr_bpar   <- tr_est0[(np_b - 1):np_b]
          if( jackup  ){
            tr_bpar[2]<- abs(tr_bpar[2])
          } else {
            tr_bpar[2]<- 0
          }
          if( !is.null( bc_value ) ){
            tr_bpar[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
          } else {
            tr_bpar[1]<- ifelse( tr_bpar[1] < -5, -5, ifelse( tr_bpar[1] > 5, 5, tr_bpar[1]) )
          }

        } else {
          tr_bpar   <- NULL
        }

        z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,bc_par=tr_bpar,jackup=jackup,y_added2=y_added2)
        y        <- z0$y
        z_ms     <- z0$z_ms
        y_ms     <- z0$y_ms

        d_abs   <- abs( d_sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
                                bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup) )
        tr_comp <- sum( log( d_abs ) )
        if( !is.null(weight)){
          w_y <- weight * y
          mo	<- crossprod( X0_org, w_y )
          m	  <- crossprod( X_org, w_y )
          yy  <- sum( y * w_y )
          parVmax_sq <- sqrt( sd( y ) / sd( y - X0_org %*% ( solve( Mo ) %*% mo ) ) )
        } else {
          mo	<- crossprod( X0, y )
          m	  <- crossprod( X, y )
          yy  <- sum( y ^ 2 )
          parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )
        }

        #print(tr_comp)
      }

    } else if(y_nonneg){ ####### without sal layer (only boxcox)

      if(!is.null(bc_value)){
        bc0_par<-c(bc_value,0.01)
      } else {
        bc0_par<-c(1,0.01)
      }
      res     <-optim(par=bc0_par,fn=lik_resf_vc_tr, k=tr_num, y_nonneg=y_nonneg,noconst_last=noconst_last,
                      evSqrt=evSqrt, y0=y0, n=n, nx=nx, emet=method, X=X,M0=MM0,Minv=MMinv,term1=term1,
                      null_dum3=null_dum3,jackup=jackup, weight=weight, plim=plim,y_added2=y_added2,
                      y_type=y_type)#interval=c(-5, 5),
      if( !is.null( bc_value ) ){
        res$par[1]<- ifelse( bc_value < -5, -5, ifelse( bc_value > 5, 5, bc_value) )
      } else {
        res$par[1]<- ifelse( res$par[1] < -5, -5, ifelse( res$par[1] > 5, 5, res$par[1]) )
      }

      score_tr  <-res$value + pen * (nx+sum(par0!=0)+1+tr_npar)
      if( !is.null(weight )) score_tr <- score_tr - sum( log( weight ))

      if(score_tr < score){
        res_old	<- res$value
        loglik	<- (-1/2)*c( res$value )
        score   <- score_tr
        tr_bpar <- res$par#minimum
        if(jackup){
          tr_bpar[2]<-abs(tr_bpar[2])
        } else {
          tr_bpar[2]<-0
        }
        tr_comp <- sum( log(d_bc(par=tr_bpar,y=y0,jackup=jackup)) )
        y       <- bc(par=tr_bpar, y=y0,jackup=jackup ) + y_added2

        if( !is.null(weight) ) {
          w_y <- weight * y
          mo	<- crossprod( X0_org, w_y )
          m	  <- crossprod( X_org , w_y )
          yy  <- sum( y * w_y )
          parVmax_sq <- sqrt( sd( y ) / sd( y - X0_org %*% ( solve( Mo ) %*% mo ) ) )

        } else {
          mo	<- crossprod( X0, y )
          m	  <- crossprod( X, y )
          yy  <- sum( y ^ 2 )
          parVmax_sq <- sqrt( sd( y ) / sd( y - X0 %*% ( solve( Mo ) %*% mo ) ) )
        }
      }
    } else {

      res_old	  <- res$value
      loglik 	  <- ( -1 / 2 ) * res_old
      score     <- res_old + pen * (nx+sum(par0!=0)+1+tr_npar)
      if( !is.null(weight)) score     <- score - sum( log( weight ))

    }

    ###############################################

    if( iter > 1 ){
      if( sum( n_omit ) == 0 ){
        obj 	<- abs( loglik - loglik0 )
      } else {
        obj	<- Inf
      }
    } else {
      obj   <- Inf
    }

    Par	    <- rbind( Par, par0 )
    loglik0	<- loglik
    LL	    <- c( LL, loglik )
    #if( !is.null(x) ) print( round( loglik, 3 ) )
    if( !is.null(x) ) print( paste( toupper( penalty ), ": ", round( score, 3 ), sep = "" ) )
    iter	  <- iter + 1
  }


  weight_lik<-weight

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

  eevSqrt		<- evSqrt[ null_dum3[ -( 1:nx )] ]
  MM		    <- M [ null_dum3, null_dum3 ]
  MM[ -( 1:nx ), -( 1:nx ) ]	<- t( MM[ -( 1:nx ), -( 1:nx ) ] * eevSqrt ) * eevSqrt
  MM[ -( 1:nx ),    1:nx   ]	<-    MM[ -( 1:nx ),    1:nx   ] * eevSqrt
  MM[    1:nx  , -( 1:nx ) ]	<- t( MM[ -( 1:nx ),    1:nx   ] )
  diag( MM [ -( 1:nx ), -( 1:nx ) ] ) <- diag( MM[ -( 1:nx ), -( 1:nx ) ] ) + 1

  MMinv	  <- solve( MM, tol = 1e-30 )
  MM0	    <- M0[ null_dum3, null_dum3 ]
  idd	    <- id[ null_dum3 ]

  if( method == "reml" ){
    term1	<- determinant( MM )$modulus
  } else if( method == "ml" ){
    term1	<- determinant( as.matrix( MM[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
  }

  mm		          <- m [ null_dum3 ]
  mm[ -( 1:nx ) ]	<- mm[ -( 1:nx ) ] * eevSqrt
  b               <- MMinv %*% mm
  b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * eevSqrt

  X3        <- X
  X3[, null_dum3 ][,-( 1:nx )]<- t(t(X[, null_dum3 ][,-( 1:nx )])* eevSqrt)
  nxx       <- nx

  ############################################# Poisson estimator
  #ulim      <- NULL
  if( y_type == "count" ){
    pred0_lm   <- X_org[ , null_dum3 ] %*% b

    if( tr_num > 0 ){
      z0_p     <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
                        bc_par=tr_bpar,jackup=jackup,y_added2=y_added2)#
      z_ms     <- z0_p$z_ms
      y_ms     <- z0_p$y_ms
      #ulim     <- max(z0_p$y)
      pred0_ex <- i_sal_k(par=tr_par,y=pred0_lm,k=tr_num,noconst_last=noconst_last,
                          bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)#y_added2=y_added2,ulim=ulim
    } else {
      pred0_ex  <-pred0_lm
      #ulim      <- max(bc(par=tr_bpar,y=y0,jackup=jackup))
      #pred0_ex  <- i_bc(par=tr_bpar,y=pred0_lm,jackup=jackup)#,ulim=ulim
    }

    pred0_ex   <- exp(pred0_ex)
    if( !is.null( offset ) ) pred0_ex <- offset*pred0_ex

    #########
    weight00    <- c(  weight0 * pred0_ex )
    w_scale_test<- n/sum(weight00)
    weight_test <- c( weight00*w_scale_test )

    y_p        <- sqrt( weight_test ) * (pred0_lm + (y_org - pred0_ex)/pred0_ex)
    X_test     <- sqrt( weight_test ) * X_org
    M_test     <- crossprod( X_test )
    m_test	   <- crossprod( X_test , y_p )

    MM_test		  <- M_test[ null_dum3, null_dum3 ]
    MM_test[ -( 1:nx ), -( 1:nx ) ]	<- t( MM_test[ -( 1:nx ), -( 1:nx ) ] * eevSqrt ) * eevSqrt
    MM_test[ -( 1:nx ),    1:nx   ]	<-    MM_test[ -( 1:nx ),    1:nx   ] * eevSqrt
    MM_test[    1:nx  , -( 1:nx ) ]	<- t( MM_test[ -( 1:nx ),    1:nx   ] )
    diag( MM_test [ -( 1:nx ), -( 1:nx ) ] ) <- diag( MM_test[ -( 1:nx ), -( 1:nx ) ] ) + 1
    MMinv_test	  <- solve( MM_test, tol = 1e-30 )
    mm_test		            <- m_test [ null_dum3 ]
    mm_test[ -( 1:nx ) ]	<- mm_test[ -( 1:nx ) ] * eevSqrt
    b_test                <- MMinv_test %*% mm_test            #u
    b_test[ -( 1:nx ) ]	  <- b_test[ -( 1:nx ) ] * eevSqrt#Vu

    X3_test        <- X_test
    X3_test[, null_dum3 ][,-( 1:nx )]<- t(t(X_test[, null_dum3 ][,-( 1:nx )])* eevSqrt)

    pred0_lm_test  <- X_org[ , null_dum3 ] %*% b_test

    if( tr_num > 0 ){
      pred0_ex_test <- i_sal_k(par=tr_par,y=pred0_lm_test,k=tr_num,noconst_last=noconst_last,
                               bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)#y_added2=y_added2,ulim=ulim
    } else {
      pred0_ex_test <- pred0_lm_test
      #pred0_ex_test <- i_bc(par=tr_bpar,y=pred0_lm_test,jackup=jackup)# - y_added,ulim=ulim
    }

    pred0_ex_test   <- exp(pred0_ex_test)
    if( !is.null( offset ) ) pred0_ex_test <- offset*pred0_ex_test

    ######## devance
    nz_dum     <- (y_org > 0) & ( pred0_ex_test[,1] > 0) & ( pred0_ex[,1] > 0)
    dev1        <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/pred0_ex[nz_dum])) - 2*sum(y_org - pred0_ex)
    dev1_test   <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/pred0_ex_test[nz_dum])) - 2*sum(y_org - pred0_ex_test)
    if(dev1 >= dev1_test){
      X    <- X_test
      X3   <- X3_test
      M    <- m_test
      m    <- m_test
      MM   <- MM_test
      mm   <- mm_test
      MMinv<- MMinv_test
      weight<-weight_test
      w_scale<-w_scale_test
      b    <- b_test
      pred0_ex<-pred0_ex_test
      dev1 <- dev1_test
    }

    ######## null deviance
    y_mean    <- mean( y_org[nz_dum] )
    dev0      <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/y_mean )) - 2*sum(y_org - y_mean )
    dev_rat0  <- (sum( dev0 ) - sum( dev1 ))/sum( dev0 )
    dev_rat   <- ifelse( dev_rat0 < 0, 0, dev_rat0)

    ######## mean squared error
    rmse      <- sqrt( mean((y_org - pred0_ex)^2) )
  }


  if( y_type == "continuous"){
    if( !is.null( weight ) ){
      pred0	  <- X_org[ , null_dum3 ] %*% b
      resid	  <- y - pred0
      SSE     <- sum( resid ^ 2 )
      SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
      sig_org <- SSE /( n - nxx )

      wsq_y   <- sqrt(weight)*y
      resid_w <- wsq_y - X[ , null_dum3 ] %*% b
      SSE     <- sum( resid_w ^ 2 )
      SSY		  <- sum( ( wsq_y - mean( wsq_y ) ) ^ 2 )
      sig		  <- SSE /( n - nxx )
      b_cov		<- sig * MMinv
      b_cov2  <- b_cov[ b!=0, b!=0 ]
      pred0_se<- sqrt( colSums( t( X3[ , null_dum3 ][ ,b!=0 ] ) *
                                  ( b_cov2 %*% t( X3[ , null_dum3 ][ ,b!=0 ] ) ) )  + sig )
      pred0_se<- pred0_se / sqrt( weight )

    } else {
      pred0	  <- X[ , null_dum3 ] %*% b
      resid	  <- y - pred0
      SSE		  <- sum( resid ^ 2 )
      SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
      sig		  <- sig_org <- SSE /( n - nxx )
      b_cov		<- sig * MMinv
      b_cov2  <- b_cov[ b!=0, b!=0 ]
      pred0_se<- sqrt( colSums( t( X3[ , null_dum3 ][, b!=0 ] ) *
                                  ( b_cov2 %*% t( X3[ , null_dum3 ][, b!=0 ] ) ) ) + sig )
    }
  } else if( y_type == "count" ){
    pred0	  <- X_org[ , null_dum3 ] %*% b
    resid	  <- y - pred0
    sig_org <- sum( ( y_org - pred0_ex )^2 / pred0_ex )/( n - nxx )

    wsq_y   <- sqrt(weight)*y
    resid_w <- wsq_y - X[ , null_dum3 ] %*% b
    SSE     <- sum( resid_w ^ 2 )
    SSY		  <- sum( ( wsq_y - mean( wsq_y ) ) ^ 2 )
    sig     <- SSE /( n - nxx )

    b_cov		<- sig * MMinv
    b_cov2  <- b_cov[ b!=0, b!=0 ]

    ######## (Gaussian approx)
    pred0_se<- sqrt( colSums( t( X3[ , null_dum3 ][ ,b!=0 ] ) *
                                ( b_cov2 %*% t( X3[ , null_dum3 ][ ,b!=0 ] ) ) ) + sig )
    pred0_se<- pred0_se
  }

  bse		  <- sqrt( diag( b_cov ) )
  nev		  <- length( ev )

  b_vc_s0	  <- NULL
  b_vc_n0	  <- NULL
  B_vc_s	  <- NULL
  B_vc_n  	<- NULL
  if( nnsv == 0 ){
    b_vc		  <- matrix( 0, nrow = n, ncol = nsv )
    bse_vc		<- matrix( 0, nrow = n, ncol = nsv )
    moran_vc  <- rep( 0, nsv )
    bb		    <- list(NULL)
    bb_cov		<- list(NULL)
    evSqrts		<- list(NULL)
    evSqrts_n <- list(NULL)
    j0		<- 1
    for( j in 1:nsv ){
      bid0		<- ifelse( j == 1, 1, nxf + j )
      if( null_dum[ j ] == 0 ){
        bid0_vc		    <- which( idd == j )
        bid		        <- c( bid0, bid0_vc )
        moran_vc[ j ] <- sum(b[ bid0_vc ]^2*ev)/(ev[1]*sum(b[ bid0_vc ]^2))
        b_vc[ , j ]	  <- b[ bid0 ] + meig$sf %*% b[ bid0_vc ]
        b_cov_sub	    <- b_cov[ bid, bid ]
        sf2		        <- t( t( meig$sf ) * evSqrt2[[ j ]] )
        x_sf		      <- as.matrix( cbind( 1, sf2 ) )
        bse_vc[ , j ]	<- sqrt( colSums( t( x_sf ) * ( b_cov_sub %*% t( x_sf ) ) ) )
        bb[[ j ]]	    <- c(b[ bid0 ], b[ bid0_vc ] )
        bb_cov[[ j ]]	<- b_cov_sub
        evSqrts[[ j ]]<- evSqrt2[[ j ]]
        j0		        <- j0 + 1
      } else {
        b_vc[ , j ]	  <- b[ bid0 ]
        moran_vc[ j ] <- NA
        bse_vc[ , j ]	<- sqrt( b_cov[ bid0, bid0 ] )
        bb[[ j ]]	    <- b[ bid0 ]
        bb_cov[[ j ]]	<- b_cov[ bid0, bid0 ]
        evSqrts[[ j ]]<- 0
      }
    }
  } else {
    b_vc_s0		<- matrix( 0, nrow = n, ncol = nsv )
    b_vc_n0	  <- matrix( 0, nrow = n, ncol = nsv )
    bse_vc_s0	<- matrix( 0, nrow = n, ncol = nsv )
    bse_vc_n0	<- matrix( 0, nrow = n, ncol = nsv )
    b_vc		  <- matrix( 0, nrow = n, ncol = nsv )
    bse_vc		<- matrix( 0, nrow = n, ncol = nsv )
    moran_vc  <- rep( 0, nsv )
    bb		    <- list(NULL)
    bb_cov		<- list(NULL)
    evSqrts	  <- list(NULL)
    evSqrts_n <- list(NULL)
    j_nvc     <- 1
    for( j in 1:nsv ){
      bid0		<- ifelse( j == 1, 1, nxf + j )
      bse_vc_n0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] )## added
      if( null_dum[ j ] == 0 ){
        bid0_vc		    <- which( idd == j )
        b_vc_s0[ , j ]<- meig$sf %*% b[ bid0_vc ]
        moran_vc[ j ] <- sum(b[ bid0_vc ]^2*ev)/(ev[1]*sum(b[ bid0_vc ]^2))

        x_sf_s        <- t( t( meig$sf ) * evSqrt2[[ j ]] )
        x_sf_ss       <- as.matrix( cbind( 1, x_sf_s ))              #### added
        b_cov_sub_s0  <- b_cov[ c(bid0, bid0_vc), c(bid0, bid0_vc) ] #### changed
        bse_vc_s0[ , j ]<- sqrt( colSums( t( x_sf_ss ) * ( b_cov_sub_s0 %*% t( x_sf_ss ) ) ) )

      } else {
        bid0_vc       <- NULL
        bse_vc_s0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] )
        x_sf_s        <- NULL
        moran_vc[ j ] <- NA
      }

      if((( 1:nsv ) %in% ( nvc_x + 1 ))[ j ] ){
        if( ( null_dum[ nsv + ng + nnxf + j - 1 ] == 0 )&(sum(idd == nsv + ng + nnxf + j - 1) > 0) ){ # changed
          bid0_nvc      <- which(idd == nsv + ng + nnxf + j - 1 ) # changed
          b_vc_n0[ , j ]<- XX1_0[[ j - 1 ]] %*% b[ bid0_nvc ]

          x_sf_n        <- t( t( XX1_0[[ j - 1 ]] ) * evSqrt2[[ nsv + ng + nnxf + j - 1 ]] )
          x_sf_nn       <- as.matrix( cbind( 1, x_sf_n ))                 #### added
          b_cov_sub_n0  <- b_cov[ c(bid0, bid0_nvc), c(bid0, bid0_nvc) ] #### changed
          bse_vc_n0[ , j ]<- sqrt( colSums( t( x_sf_nn ) * ( b_cov_sub_n0 %*% t( x_sf_nn ) ) ) )
          j_nvc_id      <-1

        } else {
          bid0_nvc <- NULL
          bse_vc_n0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] )
          x_sf_n   <- NULL
          j_nvc_id <-0
        }
        #j_nvc_id   <-0
      } else {
        bid0_nvc   <- NULL
        x_sf_n     <- NULL
        j_nvc_id <-0
      }

      if( (null_dum[ j ] == 0)|(null_dum[ nsv + ng +nnxf + j - 1 ] == 0 )){#minus 1 because line 804
        b_vc[ , j ]   <- b[ bid0 ] + b_vc_s0[ , j ] + b_vc_n0[ , j ]
        bid		        <- c( bid0, bid0_vc, bid0_nvc )
        b_cov_sub	    <- b_cov[ bid, bid ]
        x_sf		      <- as.matrix( cbind( 1,x_sf_s, x_sf_n))
        bse_vc[ , j ]	<- sqrt( colSums( t( x_sf ) * ( b_cov_sub %*% t( x_sf ) ) ) )
        bb[[ j ]]	    <- b[ bid ]
        bb_cov[[ j ]]	<- b_cov_sub
        evSqrts[[ j ]]<- evSqrt2[[ j ]]
        if( j_nvc_id == 1 ){
          evSqrts_n[[ j ]]<- evSqrt2[[ nsv + ng + nnxf + j - 1 ]]
        } else {
          evSqrts_n[[ j ]]<- 0
        }
        if( null_dum[ j ] == 1 ) evSqrts[[ j ]]<- 0

        b_vc_s0[ , j ]<-b[ bid0 ] +b_vc_s0[ , j ]   #### added
        b_vc_n0[ , j ]<-b[ bid0 ] +b_vc_n0[ , j ]   #### added
      } else {
        b_vc[ , j ]	  <- b_vc_s0[ , j ]  <- b_vc_n0[ , j ]  <- b[ bid0 ] ##### added
        bse_vc[ , j ]	<- bse_vc_s0[ , j ]<- bse_vc_n0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] ) #### added
        bb[[ j ]]	  <- b[ bid0 ]
        bb_cov[[ j ]]	<- b_cov[ bid0, bid0 ]
        evSqrts  [[ j ]]	<- 0
        evSqrts_n[[ j ]]	<- 0
      }
    }
  }

  if( nnxf != 0 ){
    bf_vc		  <- matrix( 0, nrow = n, ncol = nxf )
    bfse_vc		<- matrix( 0, nrow = n, ncol = nxf )
    bf_s		  <- list(NULL)
    bf_covs		<- list(NULL)
    evSqrts_c	<- list(NULL)
    for( j in which( (1:nxf) %in% nvc_xconst ) ){ # check!
      #bid0		        <- 1 + j
      if( (null_dum[ nsv + ng + j ] == 0)&( sum(idd == nsv + ng + j)>0 ) ){
        bid0_nvc      <- which( idd == nsv + ng + j )
        bf_vc[ , j ]  <- b[ 1 + j ] + XXconst_0[[ j ]] %*% b[ bid0_nvc ]#bid0
        x_sf_n        <- t( t( XXconst_0[[ j ]] ) * evSqrt2[[ nsv + ng + j ]] )
        bid           <-c(1 + j, bid0_nvc)
        b_cov_sub	    <- b_cov[ bid, bid ]
        x_sf		      <- as.matrix( cbind( 1,x_sf_n ))
        bfse_vc[ , j ]<- sqrt( colSums( t( x_sf ) * ( b_cov_sub %*% t( x_sf ) ) ) )

      } else {
        bid            <- 1 + j#bid0
        bf_vc  [ , j ] <- b[ bid ]
        b_cov_sub      <- sqrt( b_cov[ bid, bid ])
        bfse_vc[ , j ] <- b_cov_sub
        bid0_nvc       <- NULL
        x_sf_n         <- NULL
      }

      bf_s[[ j ]]	    <- b[ bid ]
      bf_covs[[ j ]]	<- b_cov_sub
      evSqrts_c[[ j ]]<- evSqrt2[[ nsv + ng  + j ]]
      if( null_dum[ nsv + ng + j ] == 1 ) evSqrts_c[[ j ]]<- 0  #added
    }
  } else {
    bf_vc		  <- NULL
    bfse_vc		<- NULL
    bf_s		  <- NULL
    bf_covs		<- NULL
    evSqrts_c	<- NULL
  }

  parV		<- par2[   1:nsv  ]
  parR		<- par2[ (nsv + 1):( 2 * nsv ) ]

  Bias  <-0
  if( ng != 0 ){
    id2   <- id[( id %in% which(parV!=0) )|( id > nsv )|( id == 0 ) ]
    b_g		<- b  [ id2 > nsv ]
    bse_g <- bse[ id2 > nsv ]
    bt_g	<- b_g / bse_g
    bpar_g		<- data.frame( Estimate = b_g, SE = bse_g, t_value = bt_g )

    glist <-unique(id2[id2 > nsv])
    bpar_g2<-list(NULL)
    nulldat<-data.frame(Estimate = 0, SE = NA, t_value = NA)
    for(ggid in 1:ng){
      bpar_g2_0		<- bpar_g[id2[id2 > nsv] == glist[ggid],]
      xg_idid_0   <- xg_idid2[ggid]
      if(xg_idid_0 == 1){
        bpar_g2_0		<- rbind(nulldat,bpar_g2_0)
      } else if( xg_idid_0 == dim(bpar_g2_0)[1] ){
        bpar_g2_0		<- rbind(bpar_g2_0,nulldat)
      } else {
        bpar_g2_0		<- rbind(bpar_g2_0[1:(xg_idid_0-1),],nulldat,bpar_g2_0[-(1:(xg_idid_0-1)),])
      }
      bias  <- mean(bpar_g2_0[xg_link_id[[ggid]],1])
      bpar_g2_0[,1]<- bpar_g2_0[,1] - bias
      Bias  <- Bias + bias

      bpar_g2_0[is.na(bpar_g2_0[,2])==FALSE,3]<- bpar_g2_0[is.na(bpar_g2_0[,2])==FALSE,1]/bpar_g2_0[is.na(bpar_g2_0[,2])==FALSE,2]
      xg00	<- factor( xgroup[ , ggid ] )
      rownames(bpar_g2_0) <- paste( names( as.data.frame(xgroup) )[ ggid ], "_", levels( xg00 ), sep = "" )
      bpar_g2[[ggid]] <- bpar_g2_0
    }

  } else {
    bpar_g2  <- NULL
  }

  if( ng != 0 ){
    parG	<- par2[ (2 * nsv + 1 ):( 2 * nsv + ng ) ]
  } else {
    parG  <- NULL
  }
  if(nnxf != 0 ){
    parNxf<- par2[ (2 * nsv +ng + 1 ):( 2 * nsv + ng + nnxf) ]
  } else {
    parNxf<- NULL
  }
  if(nnsv != 0 ){
    parNsv<- par2[ (2 * nsv + ng + nnxf + 1 ):( 2 * nsv + ng + nnxf + nnsv ) ]
  } else {
    parNsv<- NULL
  }

  nsv2		<- sum( parV != 0 )
  nnxf2   <- sum( parNxf != 0 )
  nnsv2   <- sum( parNsv != 0 )
  Xm		  <- X[ , null_dum3 ]
  Xm[ , -( 1:nx ) ] <- t( t( Xm[ , -( 1:nx ) ] ) * eevSqrt )
  np		  <- nxx + nsv2 * 2 + ng + nnxf2 + nnsv2 + 1 + tr_npar

  if( !is.null(weight) ){
    loglik<-loglik + 0.5*sum( log( weight_lik ))########## check!!
    Xm_org    <- X_org[ , null_dum3 ]
    Xm_org[ , -( 1:nx ) ] <- t( t( Xm_org[ , -( 1:nx ) ] ) * eevSqrt )
  }
  AIC		  <- -2 * loglik + np * 2
  BIC		  <- -2 * loglik + np * log( n )
  r2_0		<- 1 - SSE / SSY
  r2		  <- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1 )

  if( is.null( x ) ){
    b_vc[,1]<- b_vc[,1] + Bias
    b_c		  <- b[ 1:( nxf + 1 ) ]
    b_c[1]  <- b_c[1] + Bias
    r       <- b[( nx + 1 ):( nx + nev) ]
    bse_c	  <- bse[ 1:( nxf + 1 ) ]
    bt_c	  <- b_c / bse_c
    if( !is.null(weight)){
      df	    <- sum( t( Xm_org ) * ( MMinv %*% t( weight * Xm_org  ) ) )
    } else {
      df	    <- sum( t( Xm ) * ( MMinv %*% t( Xm ) ) )
    }
    bp_c	  <- 2 - 2 * pt( abs( bt_c ), df = n - df )
    b_par		<- data.frame( Estimate = b_c, SE = bse_c, t_value = bt_c, p_value = bp_c )
    rownames( b_par )<- c("(Intercept)", xfname)

  } else {
    if( nxf != 0 ) {
      b_c		<- b  [ 2:( nxf + 1 ) ]
      bse_c	<- bse[ 2:( nxf + 1 ) ]
      bt_c	<- b_c / bse_c
      if( !is.null( weight ) ){
        df	    <- sum( t( Xm_org ) * ( MMinv %*% t( weight * Xm_org  ) ) )
      } else {
        df	    <- sum( t( Xm ) * ( MMinv %*% t( Xm ) ) )
      }
      bp_c	<- 2 - 2 * pt( abs( bt_c ), df = n - df )
      b_par	<- data.frame( Estimate = b_c, SE = bse_c, t_value = bt_c, p_value = bp_c )
      rownames( b_par )<- xfname
    } else {
      if( !is.null( weight ) ){
        df	    <- sum( t( Xm_org ) * ( MMinv %*% t( weight * Xm_org  ) ) )
      } else {
        df	    <- sum( t( Xm ) * ( MMinv %*% t( Xm ) ) )
      }
      b_par	<- NULL
    }

    b_vc[,1]<- b_vc[,1] + Bias
    r       <- NULL
  }

  bt_vc		<- as.matrix( b_vc / bse_vc )
  bp_vc		<- 2 - 2 * pt( abs( bt_vc ), df = n - df )
  b_vc		<- data.frame( b_vc )
  bse_vc		<- data.frame( bse_vc )
  bt_vc		<- data.frame( bt_vc )
  bp_vc		<- data.frame( bp_vc )
  names( b_vc )	  <- c( "(Intercept)", xname )
  names( bse_vc )	<- c( "(Intercept)", xname )
  names( bt_vc )	<- c( "(Intercept)", xname )
  names( bp_vc )	<- c( "(Intercept)", xname )

  if(is.null(b_vc_s0) ==FALSE ){
    bt_vc_s0		    <- as.matrix( b_vc_s0 / bse_vc_s0 )
    bp_vc_s0        <- 2 - 2 * pt( abs( bt_vc_s0 ), df = n - df )

    b_vc_s0		      <- data.frame( b_vc_s0 )
    bse_vc_s0       <- data.frame( bse_vc_s0 )
    bt_vc_s0		    <- data.frame( bt_vc_s0 )
    bp_vc_s0		    <- data.frame( bp_vc_s0 )
    names( b_vc_s0 )	<- c( "(Intercept)", xname )
    names( bse_vc_s0 )<- c( "(Intercept)", xname )
    names( bt_vc_s0 )	<- c( "(Intercept)", xname )
    names( bp_vc_s0 )	<- c( "(Intercept)", xname )

    B_vc_s <-list(b_vc_s0, bse_vc_s0, bt_vc_s0, bp_vc_s0)
  }
  if(is.null(b_vc_n0) ==FALSE ){
    bt_vc_n0		    <- as.matrix( b_vc_n0 / bse_vc_n0 )
    bp_vc_n0        <- 2 - 2 * pt( abs( bt_vc_n0 ), df = n - df )

    b_vc_n0		      <- data.frame( b_vc_n0 )
    bse_vc_n0       <- data.frame( bse_vc_n0 )
    bt_vc_n0		    <- data.frame( bt_vc_n0 )
    bp_vc_n0		    <- data.frame( bp_vc_n0 )
    names( b_vc_n0 )	<- c( "(Intercept)", xname )
    names( bse_vc_n0 )<- c( "(Intercept)", xname )
    names( bt_vc_n0 )	<- c( "(Intercept)", xname )
    names( bp_vc_n0 )	<- c( "(Intercept)", xname )

    B_vc_n <-list(b_vc_n0, bse_vc_n0, bt_vc_n0, bp_vc_n0)
  }

  if( nnxf != 0 ) {
    bft_vc		<- as.matrix( bf_vc / bfse_vc )
    bfp_vc		<- 2 - 2 * pt( abs( bft_vc ), df = n - df )
    bf_vc		<- data.frame( bf_vc )
    bfse_vc		<- data.frame( bfse_vc )
    bft_vc		<- data.frame( bft_vc )
    bfp_vc		<- data.frame( bfp_vc )
    names( bf_vc )	  <- xxfname
    names( bfse_vc )	<- xxfname
    names( bft_vc )	<- xxfname
    names( bfp_vc )	<- xxfname
  } else {
    bf_vc    <- NULL
    bfse_vc  <- NULL
    bft_vc   <- NULL
    bfp_vc   <- NULL
  }

  parV		        <- parV * sqrt( sig )
  sf_par		      <- data.frame( rbind( parV, moran_vc ) )
  names( sf_par )	<- c( "(Intercept)", xname )
  rownames( sf_par )<- c( "random_SE", "Moran.I/max(Moran.I)" )

  s_g             <- NULL
  if( ng != 0 ){
    parG		<- t( parG * sqrt( sig ) )
    bg_par	<- data.frame( parG )
    rownames( bg_par )	<- "ramdom_SE"
    names( bg_par )	<- names( as.data.frame(xgroup) )
    s_g     <- bg_par
  }

  vc		          <- data.frame(ifelse( sf_par[1,] ==0, 0, 1) )
  names( vc )	    <- names( sf_par )
  rownames( vc )	<- "Spatial"

  if( is.null( parNsv ) ==FALSE ){
    vc  <- rbind( vc, c(0, ifelse( parNsv ==0, 0, 1 )))
    rownames( vc )[2]	<- "Non-spatial"

    parN		<- t( parNsv * sqrt( sig ) )
    bn_par	<- data.frame( parN )
    rownames( bn_par )	<- "random_SE"
    names( bn_par )	<- xxname
    s_n    <- bn_par

    totv   <-as.vector(sf_par[1,]+c(0,s_n))
    s_rat  <-as.vector(sf_par[1,])/totv
    n_rat  <-as.vector(c(0,s_n))/totv
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



  ################ pred_quantile
  pquant        <- c(0.01, 0.025, 0.05, seq(0.1,0.9,0.1), 0.95, 0.975, 0.99)
  pquant2       <- seq(0.001,0.999,0.001)

  pq_dat0       <- NULL
  for(pq in pquant){
    pq_dat0     <- cbind(pq_dat0,qnorm(pq,pred0,pred0_se))
  }
  pq_dat0       <- as.data.frame(pq_dat0)
  names(pq_dat0)<- paste("q",pquant,sep="")

  if( y_nonneg==TRUE ){
    tr_bpar          <- data.frame(tr_bpar[1])
    names(tr_bpar)   <- "Estimates"
    rownames(tr_bpar)<- "Lambda (Box-Cox)"
  } else {
    tr_bpar          <- NULL
  }

  if( tr_num > 0 ){
    z0       <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
                      bc_par=tr_bpar$Estimates,jackup=jackup,y_added2=y_added2)#
    z_ms     <- z0$z_ms
    y_ms     <- z0$y_ms
    dif      <- 1/d_sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
                          bc_par=tr_bpar$Estimates,y_ms=y_ms,z_ms=z_ms,jackup=jackup)
    pred     <- i_sal_k(par=tr_par,y=pred0,k=tr_num,noconst_last=noconst_last,
                        bc_par=tr_bpar$Estimates,y_ms=y_ms,z_ms=z_ms,jackup=jackup,y_added2=y_added2)
    # - y_added,ulim=ulim

    pq_dat   <- pq_dat0
    for(pq in 1:ncol(pq_dat0)){
      ptest<-try(pq_pred<- i_sal_k(par=tr_par,y=pq_dat0[,pq],k=tr_num,noconst_last=noconst_last,
                                   bc_par=tr_bpar$Estimates,y_ms=y_ms,z_ms=z_ms,jackup=jackup,
                                   y_added2=y_added2))#,  - y_added, ulim=ulim
      if( !inherits(ptest, "try-error") ){#class(ptest)!="try-error"
        pq_dat[,pq]     <- pq_pred
      } else {
        pq_dat[,pq]     <- NA
      }
    }

    min_s        <- min(y) - (max(y)-min(y))/10
    max_s        <- max(y) + (max(y)-min(y))/10
    if( y_type == "continuous"){
      snorm       <- seq(min_s, max_s, len=100000)
      snorm_d     <- dnorm( snorm )
      snorm_tr    <- i_sal_k(par=tr_par,y=snorm,k=tr_num,noconst_last=noconst_last,
                             bc_par=tr_bpar$Estimates,y_ms=y_ms,z_ms=z_ms,jackup=jackup) - y_added#, ulim=ulim,y_added2=y_added2
      y_probfun00 <- cbind(snorm_tr,snorm_d)
      y_probfun00 <- y_probfun00[is.finite(y_probfun00[,1]),]
      y_prob_sk      <- y_probfun00
      y_bin          <- seq(0.8*min(y0), 1.2*max(y0), len=1000)
      suppressWarnings(approx_fun     <- approxfun(x=y_probfun00[,1],y=y_probfun00[,2], rule=2))
      y_probfun      <- cbind( y_bin, approx_fun(y_bin))

      pred        <- data.frame(pred = pred, pred_transG = pred0, pred_transG_se = pred0_se)

    } else if( y_type =="count" ){
      snorm       <- seq(min_s, max_s, len=1000000)
      snorm_d     <- dnorm( snorm )
      snorm_d     <- c(pnorm( snorm )[1],diff(pnorm( snorm )))

      snorm_tr    <- i_sal_k(par=tr_par,y=snorm,k=tr_num,noconst_last=noconst_last,
                             bc_par=tr_bpar$Estimates,y_ms=y_ms,z_ms=z_ms,jackup=jackup)# - y_added#,y_added2=y_added2
      if( !is.null(offset)){
        snorm_tr_pred <- median( offset )*exp( snorm_tr )
        snorm_tr      <- snorm_tr_pred * exp(- (1+0.5*zrat)/( snorm_tr_pred + 0.5 ) )
      } else {
        snorm_tr_pred <- exp( snorm_tr )
        snorm_tr      <- snorm_tr_pred * exp(- (1+0.5*zrat)/( snorm_tr_pred + 0.5) )
      }
      y_probfun00 <- cbind( round( snorm_tr ), snorm_d )
      y_probfun00 <- y_probfun00[is.finite(y_probfun00[,1]),]

      y_bin          <- min( y_probfun00[,1] ):min(max(y_org*1.5),max( y_probfun00[,1] ))
      y_probfun00_id <- get.knnx(y_bin,y_probfun00[,1],k=1)$nn.index
      y_bin2         <- y_bin[y_probfun00_id]
      y_probfun      <- aggregate(y_probfun00[,2],list(y_bin2),sum)
      y_probfun[,2]  <- y_probfun[,2]/sum(y_probfun[,2])
      y_prob_sk      <- y_probfun
      y_probfun      <- y_probfun[( y_probfun[,1] >= 0.8*min(y_org) )&( y_probfun[,1] <= 1.2*max(y_org) ),]

      dif            <- dif*y_org
      pq_dat         <- exp(pq_dat)
      pred           <- data.frame(pred = exp(pred), pred_transG = pred0)
    }

  } else if( y_nonneg ==TRUE ){
    y        <- bc(par=tr_bpar$Estimates,y=y0,jackup=jackup)# +
    z_ms     <- NULL
    y_ms     <- NULL

    dif      <- 1 / ( d_bc(par=tr_bpar$Estimates,y=y0,jackup=jackup) )
    pred     <- i_bc(par=tr_bpar$Estimates,y=pred0,jackup=jackup) - y_added# - y_added2,ulim=ulim

    pq_dat   <- pq_dat0
    for(pq in 1:ncol(pq_dat0)){
      ptest  <- try(suppressWarnings(pq_pred<- i_bc(par=tr_bpar$Estimates,y=pq_dat0[,pq],jackup=jackup) - y_added))#,ulim=ulim
      if( !inherits(ptest, "try-error")){#class(ptest)!="try-error"
        pq_dat[,pq]       <-pq_pred
      } else {
        pq_dat[,pq]       <-NA
      }
    }

    min_s     <- min(y) - (max(y)-min(y))/10
    max_s     <- max(y) + (max(y)-min(y))/10
    snorm     <- seq(min_s, max_s, len=10000)
    snorm_d   <- dnorm(snorm,mean=mean(y),sd=sd(y))
    tr_bpar00 <- tr_bpar
    snorm_tr     <- i_bc(par=tr_bpar00$Estimates,y=snorm,jackup=jackup)#exp(  - y_added
    y_probfun00  <- cbind(snorm_tr,snorm_d)
    y_probfun00  <- y_probfun00[is.finite(y_probfun00[,1]),]
    y_prob_sk    <- y_probfun00

    #suppressWarnings(approx_fun<- approxfun(x=y_probfun00[,1],y=y_probfun00[,2], rule=2))
    #y_bin     <- seq(max(0, 0.8*min(y0)), 1.2*max(y0), len=1000)
    #y_probfun <- cbind( y_bin, approx_fun(y_bin))

    y_bin          <- seq(0.8*min(y0), 1.2*max(y0), len=1000)
    suppressWarnings(approx_fun     <- approxfun(x=y_probfun00[,1],y=y_probfun00[,2], rule=2))
    y_probfun      <- cbind( y_bin, approx_fun(y_bin))

    pred        <- data.frame(pred = pred, pred_transG = pred0, pred_transG_se = pred0_se)

  } else {
    dif          <- 1
    z_ms         <- NULL
    y_ms         <- NULL
    tr_par       <- NULL
    tr_bpar      <- NULL
    min_s        <- min(y) - (max(y)-min(y))/10
    max_s        <- max(y) + (max(y)-min(y))/10
    if( y_type == "continuous" ){
      pq_dat   <- pq_dat0
      pred     <- data.frame(pred = pred0, pred_se=pred0_se)
      snorm    <- seq(min_s, max_s, len=10000)
    } else if( y_type == "count" ){
      dif      <- dif*y_org
      pq_dat   <- exp( pq_dat0 )
      pred     <- data.frame(pred =exp( pred0 ), pred_transG = pred0 )
      snorm    <- seq(min_s, max_s, len=1000000)
    }

    snorm_d  <- dnorm(snorm,mean=mean(pred0),sd=sd(pred0))
    if( y_type=="count" ){

      if( !is.null(offset)){
        snorm_tr_pred <- median( offset )*exp( snorm )
        snorm_tr      <- snorm_tr_pred * exp(- (1+0.5*zrat)/( snorm_tr_pred + 0.5 ) )
      } else {
        snorm_tr_pred <- exp( snorm )
        snorm_tr      <- snorm_tr_pred * exp(- (1+0.5*zrat)/( snorm_tr_pred + 0.5) )
      }

      y_probfun   <- cbind( round( snorm_tr ), snorm_d )
      y_probfun    <- cbind(snorm_tr,snorm_d)
      y_probfun00  <- y_probfun[is.finite(y_probfun[,1]),]
      y_bin        <- min( y_probfun00[,1] ):max( y_probfun00[,1] )
      y_probfun00_id<-get.knnx(y_bin,y_probfun00[,1],k=1)$nn.index
      y_bin2       <- y_bin[y_probfun00_id]
      y_probfun    <- aggregate(y_probfun00[,2],list(y_bin2),sum)
      y_prob_sk    <- y_probfun

      y_probfun[,2]<- y_probfun[,2]/sum(y_probfun[,2])
      y_probfun    <- y_probfun[( y_probfun[,1] >= 0.8*min(y_org) )&( y_probfun[,1] <= 1.2*max(y_org) ),]

    } else {
      y_probfun00  <- cbind(snorm,snorm_d)
      y_probfun00  <- y_probfun00[is.finite(y_probfun00[,1]),]
      y_prob_sk    <- y_probfun00

      suppressWarnings(approx_fun   <- approxfun(x=y_probfun00[,1],y=y_probfun00[,2], rule=2))
      y_bin        <- seq(0.8*min(y0), 1.2*max(y0), len=1000)
      y_probfun    <- cbind( y_bin, approx_fun(y_bin))
    }
  }

  y_probfun        <- data.frame(y_probfun)
  names(y_probfun) <- c("Value","Intensity")

  cumrat<- cumsum(y_prob_sk[,2])/sum(y_prob_sk[,2])
  m1_est<- weighted.mean(x=y_prob_sk[,1],w=y_prob_sk[,2])
  m2_est<- weighted.mean(x=(y_prob_sk[,1] - m1_est)^2,w=y_prob_sk[,2])
  m3_est<- weighted.mean(x=(y_prob_sk[,1] - m1_est)^3,w=y_prob_sk[,2])
  m4_est<- weighted.mean(x=(y_prob_sk[,1] - m1_est)^4,w=y_prob_sk[,2])
  if( (( tr_num > 0 )|( y_nonneg ==TRUE ))&( y_type == "continuous" ) ){
    Skew  <- m3_est/m2_est^1.5
    Kurt  <- m4_est/m2_est^2 - 3
    skew_kurt<- data.frame(Estimates = c(Skew, Kurt))
    row.names( skew_kurt )<-c("skewness","excess kurtosis")
  } else if(y_type == "count"){
    Skew      <- m3_est/m2_est^1.5
    Kurt      <- m4_est/m2_est^2 - 3
    skew_kurt<- data.frame(Estimates= c(Skew, Kurt))
    row.names( skew_kurt )<-c("skewness","excess kurtosis")
  } else {
    Skew  <- 0
    Kurt  <- 0
    skew_kurt<- data.frame(Estimates= c(Skew, Kurt))
    row.names( skew_kurt )<-c("skewness","excess kurtosis")
  }

  if( y_type=="count" ){
    if( !is.null( offset )) {
      pred[,1] <-pred[,1]*offset
      pq_dat   <-round(pq_dat*offset)
    } else {
      pq_dat   <-round(pq_dat)
    }
  }

  if( y_type=="continuous" ){
    e_stat		 <- data.frame( stat = c( sqrt( sig_org ), r2, loglik, AIC, BIC ) )
    rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", lik_nam, "AIC", "BIC")

    ######################## Approximate inference for the NULL model
    if( is.null( weight ) ){
      mod_NULL_id0<-" )"
    } else {
      if( length(weight)==1 ) weight <-rep( weight, n )
      mod_NULL_id0<-", weights = weight )"
    }

    if( is.null( x ) & is.null( xconst ) ){
      mod_NULL   <- lm(y0 ~ 1, weights = weight )
      mod_NULL_id<- paste("lm( y ~ 1", mod_NULL_id0, sep="")
    } else if( is.null( x ) ){
      lm_dat     <- data.frame(y0, xconst)
      mod_NULL   <- lm(y0 ~ ., data=lm_dat, weights = weight )
      if( resf_flag ){
        mod_NULL_id<- paste("lm( y ~ x", mod_NULL_id0, sep="")
      } else {
        mod_NULL_id<- paste("lm( y ~ xconst", mod_NULL_id0, sep="")
      }
    } else if( is.null( xconst ) ){
      lm_dat     <- data.frame(y0, x)
      mod_NULL   <- lm(y0 ~ ., data=lm_dat, weights = weight )
      mod_NULL_id<- paste("lm( y ~ x", mod_NULL_id0, sep="")
    } else {
      lm_dat     <- data.frame(y0, x, xconst)
      mod_NULL   <- lm(y0 ~ ., data = lm_dat, weights = weight )
      mod_NULL_id<- paste("lm( y ~ x + xconst", mod_NULL_id0, sep="")
    }

    e_stat_N     <- data.frame( stat = c( logLik(mod_NULL), AIC(mod_NULL), BIC(mod_NULL) ) )
    rownames( e_stat_N ) <- c( lik_nam, "AIC", "BIC")
    e_stat_NULL  <- list(NULL)
    e_stat_NULL[[1]]<- e_stat_N
    e_stat_NULL[[2]]<- mod_NULL_id
    r2_devrat<-r2

  } else {
    lik_nam2  <- paste( "Gaussian ",lik_nam, " approximating the model",sep="")
    #if( tr_num==0 & y_nonneg==FALSE ){
    #  e_stat		<- data.frame( stat = c( sig_org, rmse, 100*dev_rat, loglik, AIC, BIC ) )
    #  rownames( e_stat ) <- c( "dispersion parameter","RMSE","deviance explained (%)", lik_nam2, "AIC", "BIC")
    #} else {
    e_stat		<- data.frame( stat = c( rmse, 100*dev_rat, loglik, AIC, BIC ) )
    rownames( e_stat ) <- c("RMSE","deviance explained (%)", lik_nam2, "AIC", "BIC")
    #}

    ######################## Approximate Gaussian likelihood for the NULL model
    if( is.null(offset) ){
      y_org2   <- log( y_org + 0.5) + y_added2
      mod_NULL_id0<-", family = poisson )"
    } else {
      y_org2<- log( ( y_org+0.5 )/offset ) +y_added2
      mod_NULL_id0<-", offset = log( offset ), family = poisson )"
    }

    if( length( weight0 ) > 1 ){
      mod_NULL_id0<- paste( ", weights = weight", mod_NULL_id0, sep="" )
    }

    if( is.null( x ) & is.null( xconst ) ){
      mod_NULL   <- lm(y_org2 ~ 1, weights = weight_lik )
      mod_NULL_id<- paste("glm( y ~ 1", mod_NULL_id0, sep="")
    } else if( is.null( x ) ){
      mod_NULL   <- lm(y_org2 ~ as.matrix( xconst ), weights = weight_lik )
      if( resf_flag ){
        mod_NULL_id<- paste("glm( y ~ x", mod_NULL_id0, sep="")
      } else {
        mod_NULL_id<- paste("glm( y ~ xconst", mod_NULL_id0, sep="")
      }
    } else if( is.null( xconst ) ){
      mod_NULL   <- lm(y_org2 ~ as.matrix( x )     , weights = weight_lik )
      mod_NULL_id<- paste("glm( y ~ x", mod_NULL_id0, sep="")
    } else {
      mod_NULL   <- lm(y_org2 ~ as.matrix( x ) + as.matrix( xconst ), weights = weight_lik )
      mod_NULL_id<- paste("glm( y ~ x + xconst", mod_NULL_id0, sep="")
    }

    e_stat_N     <- data.frame( stat = c( logLik(mod_NULL), AIC(mod_NULL), BIC(mod_NULL) ) )
    lik_nam2_NULL<- paste("Gaussian ",lik_nam, " approximating the NULL model",sep="")
    rownames( e_stat_N ) <- c( lik_nam2_NULL, "AIC", "BIC")
    e_stat_NULL  <- list(NULL)
    e_stat_NULL[[1]]<- e_stat_N
    e_stat_NULL[[2]]<- mod_NULL_id

    r2_devrat<-dev_rat0
  }

  messs   <- 0
  if( sum( n_omit ) > 5 ) {
    if(y_type=="count" & r2_devrat <= 0 ){
      message( "Note: Singular fit. Simplify the model")
    } else {
      if( r2_devrat > -0.2 ){
        message( "Note: The model is nearly singular. Consider simplifying the model" )
      } else {
        message( "Note: Singular fit. Simplify the model")
      }
      messs <- 1
    }

  } else if( y_type=="continuous" & r2_devrat < -0.2 ){
    message("Note: Singular fit. Simplify the model")
    messs <- 1
  } else if( y_type=="count" & r2_devrat <= 0 ){
    message("Note: Singular fit. Simplify the model")
    messs <- 1
  }

  if( (messs == 0)&( loglik < logLik( mod_NULL )) ){
    message( "Note: The model is nearly singular. Consider simplifying the model")
  }

  other		<- list( res_int =res_int, r = r, sf_alpha = parR, x_id = x_id, nx=nx, nxf = nxf, xf_id = xf_id, df = df,null_dum3=null_dum3,
                  b_s = bb, b_covs = bb_cov, B_covs = b_cov2, sig = sig, sig_org=sig_org ,xg_levels = xg_levels, is_weight = !is.null( weight ),
                  eevSqrt=eevSqrt, evSqrts = evSqrts, evSqrts_n = evSqrts_n, evSqrts_c = evSqrts_c, model = "resf_vc", b_c = bf_s, b_covs_c = bf_covs,
                  Bias=Bias, nvc_x=nvc_x, nvc_xconst=nvc_xconst, nvc_num = nvc_num, sel_basis_c = sel_basis_c, sel_basis_n = sel_basis_n,
                  x = x, xconst = xconst, coords = meig$other$coords, dif=dif, y = y0, tr_num=tr_num, y_nonneg = y_nonneg, y_type = y_type,
                  method=method, y_added = y_added, jackup=jackup, np=np, offset = offset, e_NULL = e_stat_NULL,
                  w_scale = w_scale,tr_comp=tr_comp)

  result    <- list( b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc, B_vc_s = B_vc_s, B_vc_n = B_vc_n,
                     c = b_par, c_vc = bf_vc, cse_vc = bfse_vc, ct_vc = bft_vc, cp_vc = bfp_vc, b_g = bpar_g2,
                     s = s, s_c=s_nxconst, s_g = s_g, vc = vc, e = e_stat, pred = pred, pred_quantile=pq_dat,
                     tr_par=tr_par,tr_bpar=tr_bpar,tr_y=y,
                     resid = resid, pdf=y_probfun, skew_kurt = skew_kurt, other = other, call = match.call() )
  class( result) <- "resf_vc"
  return( result )
}

print.resf_vc <- function(x, ...)
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

  } else if( !is.null(x$c_vc) ){
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
  if( !is.null(x$s_g) ){
    cat("\nGroup effects:\n")
    print(x$s_g)
  }

  if( !is.null(x$skew_kurt)|!is.null(x$tr_bpar) ){
    cat("\n----Estimated probability distribution of y--------------\n")
    if( !is.null(x$skew_kurt) ) print(x$skew_kurt)
    if( !is.null(x$tr_bpar) & x$other$y_type =="continuous" ){
      cat( paste("(Box-Cox parameter: ", format(x$tr_bpar[1], digits=7),")\n",sep="") )
    } else if( x$other$y_type =="count" ){
      cat( paste("(dispersion parameter: ", format(x$other$sig_org, digits=7),")\n",sep="") )
    }
  }

  cat("\n----Error statistics-------------------------------------\n")
  print(x$e)

  loglik_NULL<- x$other$e_NULL[[1]][1,1]
  AIC_NULL   <- x$other$e_NULL[[1]][2,1]
  BIC_NULL   <- x$other$e_NULL[[1]][3,1]
  mod_NULL   <- x$other$e_NULL[[2]]
  if( x$other$y_type =="continuous" ){
    ml_name    <- ifelse( x$other$method=="reml", "(r)loglik: ", "loglik: " )
    cat( paste('\nNULL model: ', mod_NULL, sep="") )
    cat(paste("\n   ",ml_name,format(loglik_NULL,digits=7),sep=""))
    cat(paste(" ( AIC: ",
              format(AIC_NULL,digits=7), ",  BIC: ", format(BIC_NULL,digits=7)," )\n",sep=""))

  } else if( x$other$y_type =="count" ){
    ml_name0 <- ifelse( x$other$method=="reml", "(r)loglik", "loglik" )
    ml_name  <- paste("\n   Gaussian ",ml_name0, " approximating the model: ",sep="")
    cat( paste('\nNULL model: ', mod_NULL, sep="") )
    cat(paste(ml_name,format(loglik_NULL,digits=7), "\n",sep=""))
    cat(paste("   ( AIC: ",
              format(AIC_NULL,digits=7), ",  BIC: ", format(BIC_NULL,digits=7)," )\n",sep=""))
  }

  #if( (x$other$method=="reml")&(x$other$y_type=="continuous") ){
  #  cat('\nNote: The AIC and BIC values are based on the restricted likelihood.')
  #  cat('\n      Use method ="ml" for comparison of models with different fixed effects (x and xconst)\n')
  #}
  #if(x$other$y_type=="count"){
  #  cat('\nNote: Error statistics of the log-Gaussian model used to approximate')
  #  cat('\n      the Poisson model are available from other$e_stat_logG\n')
  #}
  invisible(x)
}

