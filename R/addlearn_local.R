addlearn_local<-function( mod, meig0 = NULL, x0 = NULL, xconst0=NULL, xgroup0=NULL,
                          cl_num=NULL, cl=NULL, parallel=FALSE, ncores=NULL ){

  if( !is.null( mod$other$coords_z ) ) stop(" Not yet supported for modeling considering coords_z")

  {
    if( mod$other$y_nonneg ) stop( "This function is yet not supported for non-Gaussian data" )
    y       <- mod$other$y
    dif     <- mod$other$dif
    x       <- mod$other$x
    weight  <- mod$other$weight
    xconst  <- mod$other$xconst
    coords  <- mod$other$coords
    method  <- mod$other$method
    penalty <- mod$other$penalty

    if( inherits(mod, "besf_vc") ){
      if( !is.null(meig0) ) message("Note: besf_vc is not supported for spatial prediction. x0 and meig0 are ignored" )
      meig0 = NULL; x0 = NULL; xconst0=NULL
    }

    resf_vc_internal<- function( y, x, xconst = NULL, xgroup = NULL, weight = NULL, offset = NULL,
                                 x_nvc = FALSE, xconst_nvc = FALSE,
                                 x_sel = TRUE, x_nvc_sel = TRUE, xconst_nvc_sel = TRUE, nvc_num = 5,
                                 meig, method = "reml", penalty = "bic", maxiter = 30, nongauss = NULL,
                                 mod_ref ){

      n_omit    <- 0
      par0      <- mod_ref$other$par0
      null_dum3 <- mod_ref$other$null_dum3
      null_dum3[mod_ref$other$idd %in% which(mod_ref$vc==0)]<-FALSE
      omit_list <- mod_ref$other$omit_list

      bb		    <- mod_ref$other$b_s
      bb_cov		<- mod_ref$other$b_covs
      evSqrts		<- mod_ref$other$evSqrts
      evSqrts_n <- mod_ref$other$evSqrts_n
      evSqrts_c <- mod_ref$other$evSqrts_c

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
          if( min( y ) < 0 ) stop("y must be non-negative when y_nonneg = TRUE")
          bc_value<- NULL
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
        if( is.null(weight) ){
          weight0<- 1
          weight <- y_org + 0.5
        } else {
          weight0<- weight
          weight <- weight0*(y_org+0.5)
        }

      } else if( y_type =="continuous" ){
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

      if( !is.null(weight) ){############## caution
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

      #if( n > 150000 ){
      #  message( paste( "Note: besf_vc function is available for large samples. see help(besf_vc)" ) )
      #}

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
      nev0_list<-list(NULL)
      for(j in 1:nsv){
        nev0_list[[j]] <-  max( min( round( n /nsv ) - 2, sum( omit_list[[j]] )) , 3)
      }
      #nev0	<- min( round( n /nsv ) - 2, length( meig$ev ))
      #nev0  <- max(nev0, 3)
      #meig$sf	<- meig$sf[ , 1 : nev0 ]
      #meig$ev	<- meig$ev[   1 : nev0 ]

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
        for(j in 1:nsv){
          if( nev0_list[[j]] > nev0g ){
            if( nev0g < 2 ){
              stop("Too many groups. Reduce variables in xgroup")
            } else {
              message( paste( "Note: ", nev0_list[[j]] - nev0g, " eigenvectors are omitted to stablize the estimates",sep=""))
              nev0_list[[j]]   <- nev0g
              omit_list[[j]][  1:nev0g ]<-TRUE
              omit_list[[j]][-(1:nev0g)]<-FALSE
              #meig$sf	<- meig$sf[ , 1 : nev0 ]
              #meig$ev	<- meig$ev[   1 : nev0 ]
            }
          }
        }
      }

      X2	    <- meig$sf[,omit_list[[1]]]
      #ev	    <- meig$ev[ omit_list[[1]]]
      id	    <- c( rep( 0, nsv ), rep( 0, nxf ), rep( 1, length( meig$ev[ omit_list[[1]]] ) ) )
      if( nsv >= 2 ){
        for( i in 2:nsv ){
          X2  <- cbind( X2, X1[, i-1 ] * meig$sf[,omit_list[[ i ]]] )
          id	<- c( id, rep( i, length( meig$ev[ omit_list[[ i ]]] )))
        }
      }

      if( is.null( xgroup ) == FALSE ){
        xg_id0	<- nsv + 1
        xg_idid2<- NULL
        xg_link_id <- list( NULL )
        xg_levels  <- list(NULL)
        for( ff in 1:ng){
          xg00	   <- factor( xgroup[ , ff ] )
          xg_levels[[ff]] <- levels( xg00 )

          xg0	     <- sparse.model.matrix( ~ 0 + xg00 )
          xg0s	   <- Matrix::colSums( xg0 )
          xg_idid	 <- max( which( xg0s == max( xg0s ) ) )
          Xg0	     <- as.matrix( xg0[ , -xg_idid ] )
          xg_id1	 <- rep( xg_id0, dim( Xg0 )[2] )
          if( ff == 1 ){
            Xg	   <- Xg0
            xg_id	 <- xg_id1
          } else {
            Xg	   <- cbind( Xg, Xg0 )
            xg_id	 <- c( xg_id, xg_id1 )
          }
          xg_id0	 <- xg_id0 + 1
          xg_idid2 <-c(xg_idid2, xg_idid)

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
      if( (( tr_num > 0 )|( y_nonneg == TRUE ))==FALSE ){
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

      #par0	  <- rep( 1, 2 * nsv )
      #par0[1] <- parVmax_sq / 3
      #par0_est<- c( 1, 1 )
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

      evSqrt	<- NULL
      par0_sq	<- par0 ^ 2
      for( i in 1:nsv ){
        evv_i   <- meig$ev[ omit_list[[ i ]] ]
        evv  	  <- evv_i ^ par0_sq[ nsv + i ] * sum( evv_i ) / sum( evv_i ^ par0_sq[ nsv + i ] )
        evSqrt	<- c( evSqrt, par0_sq[ i ] * sqrt( evv ) )
      }

      M0	<- M
      for( j in 1:nsv ){
        id_sub	<- ( id == j )
        diag( M0[ id_sub, id_sub ] )<-diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == j ] ^ 2
      }

      if( ng != 0 ){
        #for( i2 in 1:ng ){
        #  xgg	<- rep( 1, sum( id == ( nsv + i2 ) ) )
        #  evSqrt	<- c( evSqrt, par0_sq[ 2 * nsv + i2 ] * xgg )
        #}
        #for( j2 in 1:ng ){
        #  if( j2 != ( par0_sel - nsv ) ){
        #    id_sub	<- ( id == ( j2 + nsv ) )
        #    if( sum(id_sub)==1 ){
        #      M0[ id_sub, id_sub ] <- M0[ id_sub, id_sub ] + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
        #    } else {
        #      diag( M0[ id_sub, id_sub ] )<-
        #        diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2
        #    }
        #  } else {
        #    id_sub	<- ( id == ( j2 + nsv ) )
        #    g_add   <- (1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv) ] ^ 2)/100
        #
        #    if( sum(id_sub)==1 ){
        #      M0[ id_sub, id_sub ]        <- M0[ id_sub, id_sub ] + g_add
        #    } else {
        #      diag( M0[ id_sub, id_sub ] )<- diag( M0[ id_sub, id_sub ] ) + g_add
        #    }
        #  }
        #}
      }

      if( nnxf != 0 ){
        #for( i2 in 1:nnxf ){
        #  xgg	<- rep( 1, sum( id == ( nsv+ ng + i2 ) ) )
        #  evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng ) + i2 ] * xgg )
        #}

        #for( j2 in 1:nnxf ){
        #  id_sub	<- ( id == ( j2 + nsv + ng ) )
        #  diag( M0[ id_sub, id_sub ] )<-
        #    diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng) ] ^ 2
        #}
      }

      if( nnsv != 0 ){
        #for( i2 in 1:nnsv ){
        #  xgg	<- rep( 1, sum( id == ( nsv + ng + nnxf + i2 ) ) )
        #  evSqrt	<- c( evSqrt, par0_sq[ (2 * nsv + ng + nnxf ) + i2 ] * xgg )
        #}

        #for( j2 in 1:nnsv ){
        #  if( j2 != ( par0_sel - nsv - ng - nnxf ) ){
        #    id_sub	<- ( id == ( j2 + nsv + ng + nnxf) )
        #    diag( M0[ id_sub, id_sub ] )<-
        #      diag( M0[ id_sub, id_sub ] ) + 1/evSqrt[ id[ - ( 1:nx ) ] == (j2 + nsv + ng + nnxf) ] ^ 2
        #  }
        #}
      }

      null_dum4 <-c( null_dum3[1:nx], null_dum3[-(1:nx)][ unlist(omit_list) ]) ##################### added
      MM	<- M [ null_dum4, null_dum4 ]
      MM0	<- M0[ null_dum4, null_dum4 ]
      mm	<- m [ null_dum4 ]
      idd	<- id[ null_dum4 ]
      id_omit1  <- c( diff( id ), 1)
      id_omit2_0<- which( id_omit1 != 0 )
      id_omit2  <-id_omit2_0[2:(nsv + 1)]
      id_omit1[ id_omit2_0[(id_omit2_0 %in% id_omit2)==FALSE] ] <- 0

      weight_lik<-weight

      par2   	  <- par0 ^ 2
      evSqrt	  <- NULL
      evSqrt2	  <- list(NULL)
      for( i in 1:nsv ){
        evv_i   <- meig$ev[ omit_list[[ i ]] ]
        evv  	  <- evv_i ^ par2[ nsv + i ] * sum( evv_i ) / sum( evv_i ^ par2[ nsv + i ] )
        evSqrt  <- c( evSqrt, par2[ i ] * sqrt( evv ) )
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

      eevSqrt		<- evSqrt[ null_dum4[ -( 1:nx )] ]
      MM		    <- M [ null_dum4, null_dum4 ]
      MM[ -( 1:nx ), -( 1:nx ) ]	<- t( MM[ -( 1:nx ), -( 1:nx ) ] * eevSqrt ) * eevSqrt
      MM[ -( 1:nx ),    1:nx   ]	<-    MM[ -( 1:nx ),    1:nx   ] * eevSqrt
      MM[    1:nx  , -( 1:nx ) ]	<- t( MM[ -( 1:nx ),    1:nx   ] )
      diag( MM [ -( 1:nx ), -( 1:nx ) ] ) <- diag( MM[ -( 1:nx ), -( 1:nx ) ] ) + 1

      MMinv	  <- solve( MM, tol = 1e-30 )
      MM0	    <- M0[ null_dum4, null_dum4 ]
      idd	    <- id[ null_dum4 ]

      if( method == "reml" ){
        term1	<- determinant( MM )$modulus
      } else if( method == "ml" ){
        term1	<- determinant( as.matrix( MM[ -( 1:nx ), -( 1:nx ) ] ) )$modulus
      }

      mm		          <- m [ null_dum4 ]
      mm[ -( 1:nx ) ]	<- mm[ -( 1:nx ) ] * eevSqrt
      b               <- MMinv %*% mm
      b[ -( 1:nx ) ]	<- b[ -( 1:nx ) ] * eevSqrt

      X3        <- X
      X3[, null_dum4 ][,-( 1:nx )]<- t(t(X[, null_dum4 ][,-( 1:nx )])* eevSqrt)
      nxx       <- nx

      ############################################# Poisson estimator
      #ulim      <- NULL
      if( y_type == "count" ){
        #pred0_lm   <- X_org[ , null_dum3 ] %*% b
        #
        #if( tr_num > 0 ){
        #  z0_p     <- sal_k(par=tr_par,y=y0,k=tr_num,noconst_last=noconst_last,
        #                    bc_par=tr_bpar,jackup=jackup,y_added2=y_added2)#
        #  z_ms     <- z0_p$z_ms
        #  y_ms     <- z0_p$y_ms
        #  #ulim     <- max(z0_p$y)
        #  pred0_ex <- i_sal_k(par=tr_par,y=pred0_lm,k=tr_num,noconst_last=noconst_last,
        #                      bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)
        #} else {
        #  pred0_ex  <-pred0_lm
        #}

        #pred0_ex   <- exp(pred0_ex)
        #if( !is.null( offset ) ) pred0_ex <- offset*pred0_ex

        #########
        #weight00    <- c(  weight0 * pred0_ex )
        #w_scale_test<- n/sum(weight00)
        #weight_test <- c( weight00*w_scale_test )

        #y_p        <- sqrt( weight_test ) * (pred0_lm + (y_org - pred0_ex)/pred0_ex)
        #X_test     <- sqrt( weight_test ) * X_org
        #M_test     <- crossprod( X_test )
        #m_test	   <- crossprod( X_test , y_p )

        #MM_test		  <- M_test[ null_dum3, null_dum3 ]
        #MM_test[ -( 1:nx ), -( 1:nx ) ]	<- t( MM_test[ -( 1:nx ), -( 1:nx ) ] * eevSqrt ) * eevSqrt
        #MM_test[ -( 1:nx ),    1:nx   ]	<-    MM_test[ -( 1:nx ),    1:nx   ] * eevSqrt
        #MM_test[    1:nx  , -( 1:nx ) ]	<- t( MM_test[ -( 1:nx ),    1:nx   ] )
        #diag( MM_test [ -( 1:nx ), -( 1:nx ) ] ) <- diag( MM_test[ -( 1:nx ), -( 1:nx ) ] ) + 1
        #MMinv_test	  <- solve( MM_test, tol = 1e-30 )
        #mm_test		            <- m_test [ null_dum3 ]
        #mm_test[ -( 1:nx ) ]	<- mm_test[ -( 1:nx ) ] * eevSqrt
        #b_test                <- MMinv_test %*% mm_test
        #b_test[ -( 1:nx ) ]	  <- b_test[ -( 1:nx ) ] * eevSqrt

        #X3_test        <- X_test
        #X3_test[, null_dum3 ][,-( 1:nx )]<- t(t(X_test[, null_dum3 ][,-( 1:nx )])* eevSqrt)

        #pred0_lm_test  <- X_org[ , null_dum3 ] %*% b_test

        #if( tr_num > 0 ){
        #  pred0_ex_test <- i_sal_k(par=tr_par,y=pred0_lm_test,k=tr_num,noconst_last=noconst_last,
        #                           bc_par=tr_bpar,y_ms=y_ms,z_ms=z_ms,jackup=jackup)#y_added2=y_added2,ulim=ulim
        #} else {
        #  pred0_ex_test <- pred0_lm_test
        #}

        #pred0_ex_test   <- exp(pred0_ex_test)
        #if( !is.null( offset ) ) pred0_ex_test <- offset*pred0_ex_test

        ######## devance
        #nz_dum     <- (y_org > 0) & ( pred0_ex_test[,1] > 0) & ( pred0_ex[,1] > 0)
        #dev1        <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/pred0_ex[nz_dum])) - 2*sum(y_org - pred0_ex)
        #dev1_test   <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/pred0_ex_test[nz_dum])) - 2*sum(y_org - pred0_ex_test)
        #if(dev1 >= dev1_test){
        #  X    <- X_test
        #  X3   <- X3_test
        #  M    <- m_test
        #  m    <- m_test
        #  MM   <- MM_test
        #  mm   <- mm_test
        #  MMinv<- MMinv_test
        #  weight<-weight_test
        #  w_scale<-w_scale_test
        #  b    <- b_test
        #  pred0_ex<-pred0_ex_test
        #  dev1 <- dev1_test
        #}

        ######## null deviance
        #y_mean    <- mean( y_org[nz_dum] )
        #dev0      <- 2*sum(y_org[nz_dum]*log(y_org[nz_dum]/y_mean )) - 2*sum(y_org - y_mean )
        #dev_rat0  <- (sum( dev0 ) - sum( dev1 ))/sum( dev0 )
        #dev_rat   <- ifelse( dev_rat0 < 0, 0, dev_rat0)

        ######## mean squared error
        #rmse      <- sqrt( mean((y_org - pred0_ex)^2) )
      }

      if( y_type == "continuous"){
        if( !is.null( weight ) ){
          pred0	  <- X_org[ , null_dum4 ] %*% b
          resid	  <- y - pred0
          SSE     <- sum( resid ^ 2 )
          SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
          sig_org <- SSE /( n - nxx )

          wsq_y   <- sqrt(weight)*y
          resid_w <- wsq_y - X[ , null_dum4 ] %*% b
          SSE     <- sum( resid_w ^ 2 )
          SSY		  <- sum( ( wsq_y - mean( wsq_y ) ) ^ 2 )
          sig		  <- SSE /( n - nxx )
          b_cov		<- sig * MMinv
          b_cov2  <- b_cov[ b!=0, b!=0 ]
          pred0_se<- sqrt( colSums( t( X3[ , null_dum4 ][ ,b!=0 ] ) *
                                      ( b_cov2 %*% t( X3[ , null_dum4 ][ ,b!=0 ] ) ) )  + sig )
          pred0_se<- pred0_se / sqrt( weight )

        } else {
          pred0	  <- X[ , null_dum4 ] %*% b
          resid	  <- y - pred0
          SSE		  <- sum( resid ^ 2 )
          SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
          sig		  <- sig_org <- SSE /( n - nxx )
          b_cov		<- sig * MMinv
          b_cov2  <- b_cov[ b!=0, b!=0 ]
          pred0_se<- sqrt( colSums( t( X3[ , null_dum4 ][, b!=0 ] ) *
                                      ( b_cov2 %*% t( X3[ , null_dum4 ][, b!=0 ] ) ) ) + sig )
        }
      } else if( y_type == "count" ){
        #pred0	  <- X_org[ , null_dum3 ] %*% b
        #resid	  <- y - pred0
        #sig_org <- sum( ( y_org - pred0_ex )^2 / pred0_ex )/( n - nxx )

        #wsq_y   <- sqrt(weight)*y
        #resid_w <- wsq_y - X[ , null_dum3 ] %*% b
        #SSE     <- sum( resid_w ^ 2 )
        #SSY		  <- sum( ( wsq_y - mean( wsq_y ) ) ^ 2 )
        #sig     <- SSE /( n - nxx )

        #b_cov		<- sig * MMinv
        #b_cov2  <- b_cov[ b!=0, b!=0 ]

        ######## (Gaussian approx)
        #pred0_se<- sqrt( colSums( t( X3[ , null_dum3 ][ ,b!=0 ] ) *
        #                            ( b_cov2 %*% t( X3[ , null_dum3 ][ ,b!=0 ] ) ) ) + sig )
        #pred0_se<- pred0_se
      }

      bse		  <- sqrt( diag( b_cov ) )
      #nev		  <- length( ev )

      for(pp in 1:length(null_dum)){
        null_dum[pp]<- ifelse(sum(idd==pp)==0,1,0)
      }

      b_vc_s0	  <- NULL
      b_vc_n0	  <- NULL
      B_vc_s	  <- NULL
      B_vc_n  	<- NULL

      evSqrts		    <- evSqrts_n <- list(NULL)
      evSqrts[1:nsv]<- evSqrts_n[1:nsv] <- NA

      if( nnsv == 0 ){
        b_vc		  <- matrix( 0, nrow = n, ncol = nsv )
        bse_vc		<- matrix( 0, nrow = n, ncol = nsv )
        moran_vc  <- rep( 0, nsv )
        bb		    <- list(NULL)
        bb_cov		<- list(NULL)
        j0		    <- 1
        for( j in 1:nsv ){
          bid0		<- ifelse( j == 1, 1, nxf + j )
          if( null_dum[ j ] == 0 ){
            bid0_vc		    <- which( idd == j )
            bid		        <- c( bid0, bid0_vc )
            moran_vc[ j ] <- sum(b[ bid0_vc ]^2*meig$ev[ omit_list[[j]] ])/(meig$ev[1]*sum(b[ bid0_vc ]^2)) ##########
            b_vc[ , j ]	  <- b[ bid0 ] + meig$sf[,omit_list[[j]]] %*% b[ bid0_vc ]
            b_cov_sub	    <- b_cov[ bid, bid ]
            sf2		        <- t( t( meig$sf[,omit_list[[j]]] ) * evSqrt2[[ j ]] )
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
            #evSqrts[[ j ]]<- NA
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
        j_nvc     <- 1
        for( j in 1:nsv ){
          bid0		<- ifelse( j == 1, 1, nxf + j )
          bse_vc_n0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] )## added
          if( null_dum[ j ] == 0 ){
            bid0_vc		    <- which( idd == j )
            b_vc_s0[ , j ]<- meig$sf[,omit_list[[j]]] %*% b[ bid0_vc ]
            moran_vc[ j ] <- sum(b[ bid0_vc ]^2*meig$ev[ omit_list[[j]] ])/(meig$ev[1]*sum(b[ bid0_vc ]^2))

            x_sf_s        <- t( t( meig$sf[,omit_list[[j]]] ) * evSqrt2[[ j ]] )
            x_sf_ss       <- as.matrix( cbind( 1, x_sf_s ))              #### added
            b_cov_sub_s0  <- b_cov[ c(bid0, bid0_vc), c(bid0, bid0_vc) ] #### changed
            bse_vc_s0[ , j ]<- sqrt( colSums( t( x_sf_ss ) * ( b_cov_sub_s0 %*% t( x_sf_ss ) ) ) )

          } else {
            bid0_vc       <- NULL
            bse_vc_s0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] )
            x_sf_s        <- NULL
            moran_vc[ j ] <- NA
          }

          if( (( 1:nsv ) %in% ( nvc_x + 1 ))[ j ] ){
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
            }# else {
            #  evSqrts_n[[ j ]]<- NA
            #}
            if( null_dum[ j ] == 1 ) evSqrts[[ j ]]<- NA

            b_vc_s0[ , j ]<-b[ bid0 ] +b_vc_s0[ , j ]   #### added
            b_vc_n0[ , j ]<-b[ bid0 ] +b_vc_n0[ , j ]   #### added
          } else {
            b_vc[ , j ]	  <- b_vc_s0[ , j ]  <- b_vc_n0[ , j ]  <- b[ bid0 ] ##### added
            bse_vc[ , j ]	<- bse_vc_s0[ , j ]<- bse_vc_n0[ , j ]<- sqrt( b_cov[ bid0, bid0 ] ) #### added
            bb[[ j ]]	  <- b[ bid0 ]
            bb_cov[[ j ]]	<- b_cov[ bid0, bid0 ]
            #evSqrts  [[ j ]]	<- NA
            #evSqrts_n[[ j ]]	<- NA
          }
        }
      }

      if( nnxf != 0 ){
        bf_vc		  <- matrix( 0, nrow = n, ncol = nxf )
        bfse_vc		<- matrix( 0, nrow = n, ncol = nxf )
        bf_s		  <- list(NULL)
        bf_covs		<- list(NULL)
        evSqrts_c	<- list(NULL)
        evSqrts_c[1:nxf]<-NA
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
          if( null_dum[ nsv + ng + j ] == 1 ) evSqrts_c[[ j ]]<- NA  #added
        }
      } else {
        bf_vc		  <- NULL
        bfse_vc		<- NULL
        bf_s		  <- NULL
        bf_covs		<- NULL
        #evSqrts_c	<- NA
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
      Xm		  <- X[ , null_dum4 ]
      Xm[ , -( 1:nx ) ] <- t( t( Xm[ , -( 1:nx ) ] ) * eevSqrt )
      np		  <- nxx + nsv2 * 2 + ng + nnxf2 + nnsv2 + 1 + tr_npar

      b_sub   <- b[ -( 1:nx ) ][eevSqrt>0]/eevSqrt[eevSqrt>0]

      if( is.null(weight) ){
        dd      <- sum( resid^2 )+sum( (b_sub) ^ 2 )
      } else {
        dd      <- sum( weight*resid^2 )+sum( (b_sub) ^ 2 )
      }
      if( method == "reml" ){
        term2_0<- log( 2 * pi * dd / ( n - nx ) )
        term2	 <- ( n - nx )  + ( n - nx ) * term2_0
      } else if( method == "ml" ){
        term2_0<- log( 2 * pi * dd / n )
        term2	 <- n + n * term2_0
      }
      loglik		<- loglik0<- -0.5*(term1 + term2 -2*tr_comp)

      if( !is.null(weight) ){
        loglik<-loglik + 0.5*sum( log( weight_lik ))########## check!!
        Xm_org    <- X_org[ , null_dum4 ]
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
        r       <- b[( nx + 1 ):( nx + length( meig$ev[ omit_list[[j]] ] )) ]
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

      } else if( y_nonneg ==TRUE ){

      } else {
        #dif          <- 1
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
          #dif      <- dif*y_org
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

        if( y_nonneg ==FALSE & tr_num ==0 ){
          rownames( e_stat ) <- c( "resid_SE", "adjR2(cond)", lik_nam, "AIC", "BIC")
        } else {
          rownames( e_stat ) <- c( "resid_SE(trans)", "adjR2(cond)", lik_nam, "AIC", "BIC")
        }

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

      }

      messs   <- 0
      if( sum( n_omit ) > 5 ) {
        if(y_type=="count" & r2_devrat <= -0.1 ){
          message( "Note: Singular fit. Simplify the model")
        } else {
          if( r2_devrat <= -0.1 ){
            #  message( "Note: The model is nearly singular. Consider simplifying the model" )
            #} else {
            message( "Note: Singular fit. Simplify the model")
          }
          messs <- 1
        }

      } else if( y_type=="continuous" & r2_devrat < -0.1 ){
        message("Note: Singular fit. Simplify the model")
        messs <- 1
      } else if( y_type=="count" & r2_devrat <= -0.1 ){
        message("Note: Singular fit. Simplify the model")
        messs <- 1
      }

      #if( (messs == 0)&( loglik < logLik( mod_NULL )) ){
      #  message( "Note: The model is nearly singular. Consider simplifying the model")
      #}

      evSqrts_t <- evSqrts_tc <- evSqrts_t_int <- evSqrts_tc_int<- list(NULL)
      evSqrts_t[1:nsv] <- evSqrts_tc[1:nsv] <- evSqrts_t_int[1:nsv] <- evSqrts_tc_int[1:nsv]<- NA

      other		<- list(  r = r, sf_alpha = parR, x_id = x_id, nx=nx, nxf = nxf, xf_id = xf_id, df = df,null_dum3=null_dum4,
                       b_s = bb, b_covs = bb_cov, B_covs = b_cov2, sig = sig, sig_org=sig_org ,xg_levels = xg_levels, is_weight = !is.null( weight ),
                       eevSqrt=eevSqrt, evSqrts = evSqrts, evSqrts_n = evSqrts_n, evSqrts_c = evSqrts_c, model = "resf_vc", b_c = bf_s, b_covs_c = bf_covs,
                       Bias=Bias, nvc_x=nvc_x, nvc_xconst=nvc_xconst, nvc_num = nvc_num, sel_basis_c = sel_basis_c, sel_basis_n = sel_basis_n,
                       x = x, xconst = xconst, coords = meig$other$coords, dif=dif, y = y0, tr_num=tr_num, y_nonneg = y_nonneg, y_type = y_type,
                       method=method, y_added = y_added, jackup=jackup, np=np, offset = offset, e_NULL = e_stat_NULL,
                       w_scale = w_scale,tr_comp=tr_comp,par0=par0, b_sub=b_sub,omit_list=omit_list, nev0=length(meig$ev),
                       interact_sel_all=NULL,evSqrts_t = evSqrts_t, evSqrts_tc = evSqrts_tc, evSqrts_t_int = evSqrts_tc, evSqrts_tc_int = evSqrts_tc)

      result    <- list( b_vc = b_vc, bse_vc = bse_vc, t_vc = bt_vc, p_vc = bp_vc, B_vc_s = B_vc_s, B_vc_n = B_vc_n,
                         c = b_par, c_vc = bf_vc, cse_vc = bfse_vc, ct_vc = bft_vc, cp_vc = bfp_vc, b_g = bpar_g2,
                         s = s, s_c=s_nxconst, s_g = s_g, vc = vc, e = e_stat, pred = pred, pred_quantile=pq_dat,
                         tr_par=tr_par,tr_bpar=tr_bpar,tr_y=y,
                         resid = resid, pdf=y_probfun, skew_kurt = skew_kurt, other = other, call = match.call() )
      class( result) <- "resf_vc_internal"
      return( result )
    }


    kmeans2   <- function(coords, cl_num, min_size=200, max_size=900){

      coords_uni  <- unique(coords)
      kres        <- kmeans(coords_uni,cl_num)

      cl_tab <- table(kres$cluster)
      while(max(cl_tab) > max_size){
        cc        <-which.max(cl_tab)[1]
        sub_kmeans<-kmeans(coords_uni[kres$cluster==cc,],2)
        kres$cluster[kres$cluster==cc][sub_kmeans$cluster==1]<-cc
        kres$cluster[kres$cluster==cc][sub_kmeans$cluster==2]<-max(kres$cluster)+1
        kres$centers[cc,]<-sub_kmeans$centers[1,]
        kres$centers<-rbind(kres$centers,sub_kmeans$centers[2,])
        cl_tab <- table(kres$cluster)
      }

      cl_num <-length(cl_tab)
      cl_list<-1:cl_num
      cl_tab <-table(kres$cluster)
      while(min(cl_tab) < min_size){
        cc           <-cl_list[which.min(cl_tab)[1]]
        cc_exclude   <-cl_list!=cc
        cl_list      <-cl_list[cc_exclude]
        coords_sub   <-as.matrix(coords_uni[which(kres$cluster==cc),])
        kres$centers <-kres$centers[cc_exclude,]#-cc
        if( is.null(ncol(kres$centers)) ){
          cl_nn_id     <-get.knnx(cbind(kres$centers[1],kres$centers[2]),coords_sub,1)$nn.index
        } else {
          cl_nn_id     <-get.knnx(kres$centers,coords_sub,1)$nn.index
        }
        kres$cluster[which(kres$cluster==cc)] <- cl_list[cl_nn_id]
        cl_tab <-table(kres$cluster)
      }

      cl_id0  <- sort(unique(kres$cluster))
      for(cl_id in 1:length(cl_id0)) kres$cluster[kres$cluster==cl_id0[cl_id]]<-cl_id

      kres2           <- kres
      kres2$cluster   <- kres$cluster[get.knnx(coords_uni,coords,1)$nn.index]

      return(kres2)
    }

    samp_sel  <- function(coords,coords_sub1,cut,r,cl,k){
      id_sel1    <- (coords[,1]>min(coords_sub1[,1])-cut*r)&
        (coords[,1]<max(coords_sub1[,1])+cut*r)&
        (coords[,2]>min(coords_sub1[,2])-cut*r)&
        (coords[,2]<max(coords_sub1[,2])+cut*r)& (cl!=k)
      if(sum(id_sel1)==0){
        id_sel2    <-0
        coords_sub2<-NULL
        dist_sub2  <-NULL
      } else {
        coords_sub2<-coords[id_sel1,]
        if(sum(id_sel1)==1){
          dist_sub2  <- sqrt( (coords_sub1[,1]-coords_sub2[1])^2 +
                                (coords_sub1[,2]-coords_sub2[2])^2 )
        } else {
          dist_sub2  <- get.knnx(coords_sub1,coords_sub2,1)$nn.dist
        }

        id_sel2      <- dist_sub2 < cut*r
        max_id       <- nrow( coords_sub1 )
        if( sum( id_sel2 ) > max_id ){
          id_sel2[ order(dist_sub2)[-(1:max_id)] ] <- FALSE
        }
      }
      return(list(id_sel1=id_sel1,id_sel2=id_sel2, coords_sub2=coords_sub2,
                  dist_sub2=dist_sub2))
    }

    cut          <- 2.2
    n            <- length(y)
    coords       <- as.matrix(coords)
    coords0      <- NULL
    if( is.null( meig0 ) == FALSE ) coords0 <- as.matrix(meig0$other$coords0)

    if( is.null( x ) == FALSE ){
      x	<- as.matrix( x )
      if( is.numeric( x ) == FALSE ){
        mode( x ) <- "numeric"
      }
      x_id	<- apply( x, 2, sd ) != 0
      nx	  <- sum( x_id )
      if( nx == 0 ){
        x	  <- NULL
      } else {
        x	  <- as.matrix( cbind(1, x[ , x_id ] ))
      }
      xname	<- names( as.data.frame( x ) )[-1]
    } else {
      x	  <- as.matrix( rep(1, length(y) ))
      nx    <- 0
      x_id  <- NULL
      x_nvc <- FALSE
      xname <- FALSE
    }

    if( is.null( xconst ) == FALSE ){
      xconst	        <- as.matrix( xconst )
      if( is.numeric( xconst ) == FALSE ) mode( xconst )<- "numeric"
      nxf    <- ncol(xconst)
      xconst_name	    <- names( as.data.frame( xconst ) )
    } else {
      nxf         <- 0
      xf_id       <- NULL
      xconst_name	<- NULL
    }

    id           <- 1:n
    if( !is.null( weight )) weight<-weight/sum(weight)*n

    if( is.null(cl) ){
      if( is.null( cl_num ) ) cl_num <- round(n/600)
      cl_num     <- max(cl_num, 3)
      cl         <- kmeans2(coords, cl_num=cl_num)$cluster
      table(cl)
      if(length(unique(cl))==1){
        stop("Number of cluster becomes one. Use the resf_vc function for this small data")
      }
    }

    cl_num       <- max( cl )
    cl_num2      <- cl_num+1

    Id           <- list(NULL)
    Id0          <- list(NULL)
    R            <- list(NULL)
    Wei          <- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
    for(k in 1:cl_num){ # This part is cl_num
      coords_sub1<- coords[cl==k,]
      cdist1     <- rep(0, sum(cl==k))
      r          <- max( spantree( rdist( coords_sub1 ) )$dist )+1e-15#sd(rdist( coords_sub1 ))/20#
      samp_sel1  <- samp_sel(coords=coords,coords_sub1=coords_sub1,
                             cut=cut,r=r,cl=cl,k=k)
      id_sel1    <- samp_sel1$id_sel1
      id_sel2    <- samp_sel1$id_sel2
      Id0[[k]]   <- id[cl==k]
      Id[[k]]    <- c(id[cl==k],id[id_sel1][id_sel2])
      R[[k]]     <- r

      dist_sub2  <- samp_sel1$dist_sub2
      if( sum(id_sel2)>0 ){
        cdist_sub <- c(cdist1,dist_sub2[id_sel2])
      } else {
        cdist_sub <- cdist1
      }
      Wei[Id[[k]],k]   <- exp(-(cdist_sub/r)^2)#exp(-cdist_sub/r)
    }

    Wei [,cl_num2]     <-1
    Wei     <- Wei/rowSums(Wei,na.rm=TRUE)

    if(!is.null(coords0)){
      if( is.null( x0 ) == FALSE ){
        x0	<- as.matrix( x0 )
        if( is.numeric( x0 ) == FALSE ){
          mode( x0 ) <- "numeric"
        }
        if( nx == 0 ){
          x0	  <- NULL
        } else {
          if(nrow(coords0)==1){
            x0  <- as.matrix( cbind(1, t(x0[ , x_id ]) ))
          } else {
            x0  <- as.matrix( cbind(1, x0[ , x_id ] ))
          }
        }
      }

      if( is.null( xconst0 ) == FALSE ){
        xconst0	    <- as.matrix( xconst0 )
        if( is.numeric( xconst0 ) == FALSE ) mode( xconst0 )<- "numeric"
      }

      nn         <- get.knnx(coords,coords0,4)
      nn$nn.dist[nn$nn.dist< 1e-10]<-1e-10
      nn_w       <- 1/nn$nn.dist/rowSums(1/nn$nn.dist)
      n0         <- nrow(coords0)
      Wei0       <- Matrix(nrow = nrow(coords0), ncol = cl_num2, data = 0, sparse = TRUE)
      for(k in 1:cl_num2){ # This part is cl_num
        for(kk in 1:4){
          Wei0[,k]   <- Wei0[,k] + nn_w[,kk]*Wei[nn$nn.index[,kk],k]
        }
      }
      Wei0       <- Wei0/rowSums(Wei0,na.rm=TRUE)
      Pred0      <- Pred0_se <- Wei0
      B_vc0      <- Bse_vc0  <- list(NULL)
      for(kk in 1:ncol(mod$b_vc)){
        B_vc0[[kk]]<-Bse_vc0[[kk]]<- Matrix(nrow = n0, ncol = cl_num2, data = 0, sparse = TRUE)
      }
    } else {
      Pred0      <- Pred0_se <- Wei0       <- NULL
    }

    ######################### global model
    if( !is.null(weight) ){
      weight_sub3<- weight
    } else {
      weight_sub3<- 1
    }

    Mod0        <- mod
    x_sel_global<- which( Mod0$vc[rownames(Mod0$vc)=="Spatial",]==1 )
    weight_sc   <- weight_sub3*Mod0$other$w_scale*Wei[,cl_num2]
    n_sub       <- n - ifelse(method=="reml",Mod0$other$nx,0)####changed
    nw_sub      <- sum(Wei[,cl_num2])
    nw_sub2     <- sum(Wei[,cl_num2]) - ifelse(method=="reml",Mod0$other$nx,0)
    eWe_unsc    <- sum(Wei[,cl_num2]  * Mod0$resid^2)
    eWe_sc      <- sum(weight_sc* Mod0$resid^2)
    uu          <- sum(Mod0$other$b_sub^2)
    term2_unsc  <- (nw_sub2/2)*log(2*pi*(eWe_unsc+uu)/n_sub)
    term2_sc    <- (n_sub/2)*log(2*pi*(eWe_sc+uu)/n_sub)

    pen         <- -1/2*sum(log(weight_sc)) + term2_sc - term2_unsc + n_sub/2 - nw_sub/2
    loglik      <- Mod0$e$stat[3] + pen

    if( !is.null(coords0) ){
      mod$other$x<-x
      pmod0     <- predict0(mod,x0=x0,xconst0=xconst0,xgroup0=xgroup0,
                            meig0=meig0,compute_se = TRUE,weight0=Wei0[,cl_num2])
    } else {
      pmod0     <- NA
    }

    Mod_global  <- list(Mod0,loglik, pmod0)#meig_sub,
    ######################### setting for parallelization

    #if( inherits(mod, "besf_vc")&!is.null(mod$other$ncores) ){
    #    ncores <- makeCluster(mod$other$ncores,setup_strategy = "sequential")
    #    registerDoParallel(ncores)
    #    doop    <- `%dopar%`
    #} else if(parallel){
    if(parallel){
      if(is.null(ncores)) {
        ncores00<- detectCores()
        ncores0 <- ifelse(ncores00<=2,ncores00, ifelse(ncores00<=5, ncores00-1,ncores00-2))
        ncores <- makeCluster(ncores0,setup_strategy = "sequential")
      } else {
        ncores0<- ncores
        ncores <- makeCluster(ncores0,setup_strategy = "sequential")
      }
      message(paste0("Parallelized with ", ncores0," cores"))
      registerDoParallel(ncores)
      doop     <- `%dopar%`
    } else {
      doop     <- `%do%`
    }
  }

  ######################### local models
  print(paste0("-------- Aggregating ", cl_num, " local sub-models ---------"))
  Mod  <- doop(foreach(k = 1:cl_num, .combine = "list", .multicombine=TRUE,
                       .export=c("spantree","rdist","get.knnx","meigen","resf_vc","besf_vc",#
                                 "meigen","meigen0","predict0")),{
                                   #Mod<- list(NULL)
                                   #for(k in 1:cl_num){
                                   print( paste0( k,"/",cl_num ) )

                                   coords_sub1<- coords |> subset(cl==k)
                                   y_sub1     <- y |> subset(cl==k)
                                   if(!is.null(x)){
                                     x_sub1           <- x |> subset(cl==k)
                                     x_sub_sel_local  <- c(1, which(apply(as.matrix(x_sub1),2,sd)>0)) ### similar one required for xconst?
                                     x_sub_sel        <- intersect( x_sel_global, x_sub_sel_local )
                                     x_sub_sel_const  <- setdiff(x_sub_sel_local,x_sel_global)
                                     x_sub1_const     <- x_sub1[,x_sub_sel_const]
                                   } else {
                                     x_sub1     <- x_sub_sel<- x_sub1_const<-NULL
                                   }

                                   if(!is.null(xconst)){
                                     xconst_sub1      <- xconst |> subset(cl==k)
                                     xconst_sub_sel   <- which(apply(as.matrix(xconst_sub1),2,sd)>0) ### similar one required for xconst?

                                   } else {
                                     xconst_sub1 <- xconst_sub_sel <- NULL
                                   }
                                   if( !is.null(weight) ){
                                     weight_sub1<- weight |> subset(cl==k)
                                   } else {
                                     weight_sub1<- 1
                                   }

                                   samp_sel1  <- samp_sel(coords=coords,coords_sub1=coords_sub1,
                                                          cut=cut,r=R[[k]],cl=cl,k=k)
                                   id_sel1    <- samp_sel1$id_sel1
                                   id_sel2    <- samp_sel1$id_sel2
                                   coords_sub2<- samp_sel1$coords_sub2

                                   if(!is.null(coords_sub2[id_sel2,][1])){
                                     coords_sub12<-rbind(coords_sub1,coords_sub2[id_sel2,])
                                     suppressMessages(meig_sub0 <- meig_sub <- meigen(coords=coords_sub12))
                                     meig_sub0$sf<- meig_sub$sf[1:nrow(coords_sub1),]
                                     meig_sub0$other$sfsf<- meig_sub0$other$sfsf[1:nrow(coords_sub1)]
                                     meig_sub0$other$coordk <- meig_sub$other$coordk <-coords_sub12

                                   } else {
                                     suppressMessages(meig_sub0 <-meig_sub <- meigen(coords=coords_sub1))
                                     meig_sub0$other$coordk <- meig_sub$other$coordk <-coords_sub1
                                   }


                                   #y=y_sub1;x=x_sub1[,x_sub_sel]#
                                   #xconst=cbind(xconst_sub1[,xconst_sub_sel],x_sub1_const)
                                   #meig=meig_sub0;x_sel=TRUE;miniter=1;maxiter=2#x_sel
                                   #method=method; penalty = penalty
                                   #weight=weight_sub1*Wei[Id0[[k]],k]
                                   #xgroup=NULL;x_nvc=FALSE;xconst_nvc = FALSE
                                   #offset<-NULL;x_nvc_sel<-TRUE;xconst_nvc_sel = TRUE; nvc_num = 5
                                   #nongauss = NULL
                                   invisible(capture.output(mod_sub2<- resf_vc(y=y_sub1,x=x_sub1[,x_sub_sel],#
                                                                               xconst=cbind(xconst_sub1[,xconst_sub_sel],x_sub1_const),
                                                                               meig=meig_sub0,x_sel=TRUE,miniter=1,maxiter=2,#x_sel
                                                                               method=method, penalty = penalty,
                                                                               weight=weight_sub1*Wei[Id0[[k]],k],
                                                                               xgroup=NULL,x_nvc=FALSE,xconst_nvc = FALSE)))

                                   if( sum(id_sel1)!=0 | sum(id_sel2)!=0 ){
                                     y_sub2        <- y |> subset(id_sel1) |> subset( id_sel2 ) |> c(y_sub1, x=_)
                                     if( !is.null(x) ){
                                       x_sub2      <- x |> subset(id_sel1) |> subset( id_sel2 ) |> rbind(x_sub1, x=_)
                                       x_sub2_const<- x_sub2[,x_sub_sel_const]

                                     } else {
                                       x_sub2    <- x_sub2_const<- NULL
                                     }
                                     if(!is.null(xconst)){
                                       xconst_sub2<- xconst |> subset(id_sel1) |> subset( id_sel2 ) |> rbind(xconst_sub1, x=_)

                                     } else {
                                       xconst_sub2<- NULL
                                     }
                                     if( !is.null(weight) ){
                                       weight_sub2<- weight |> subset(id_sel1) |> subset( id_sel2 ) |> c(weight_sub1, x=_)
                                     } else {
                                       weight_sub2<- 1
                                     }

                                     meig_sub$sf<-meig_sub$sf[,1:mod_sub2$other$nev0]
                                     meig_sub$ev<-meig_sub$ev[ 1:mod_sub2$other$nev0]

                                     #y=y_sub2;x=x_sub2[,x_sub_sel]
                                     #xconst=cbind(xconst_sub2[,xconst_sub_sel],x_sub2_const)
                                     #meig=meig_sub;x_sel=TRUE;mod_ref=mod_sub2
                                     #weight=weight_sub2*Wei[Id[[k]],k];method=method;penalty = penalty
                                     #xgroup=NULL;x_nvc=FALSE;xconst_nvc = FALSE
                                     #offset<-NULL;x_nvc_sel<-TRUE;xconst_nvc_sel = TRUE; nvc_num = 5
                                     #maxiter = 30; nongauss = NULL
                                     Mod0       <- resf_vc_internal(y=y_sub2,x=x_sub2[,x_sub_sel],
                                                                    xconst=cbind(xconst_sub2[,xconst_sub_sel],x_sub2_const),
                                                                    meig=meig_sub,x_sel=TRUE,mod_ref=mod_sub2,
                                                                    weight=weight_sub2*Wei[Id[[k]],k],method=method,penalty = penalty,
                                                                    xgroup=NULL,x_nvc=FALSE,xconst_nvc = FALSE)
                                     #Mod0b       <- resf_vc(y=y_sub2,x=x_sub2[,x_sub_sel],
                                     #                              xconst=cbind(xconst_sub2[,xconst_sub_sel],x_sub2_const),
                                     #                              meig=meig_sub,x_sel=TRUE,
                                     #                              weight=weight_sub2*Wei[Id[[k]],k],method=method,penalty = penalty,
                                     #                              xgroup=NULL,x_nvc=FALSE,xconst_nvc = FALSE)

                                     if( !is.null(coords0) ){
                                       sel0          <- which( Wei0[,k]>0 )
                                       if(length(sel0)>0){
                                         x0_sub        <- xconst0_sub<-x0_sub1_const<-xgroup0_sub<-NULL
                                         if( !is.null(x0))      x0_sub      <- as.matrix(x0)[sel0,x_sub_sel]
                                         if( !is.null(xconst0 ))xconst0_sub <- as.matrix(xconst0)[sel0,xconst_sub_sel]
                                         if(ncol(x0)>length(x_sub_sel)) x0_sub1_const<-as.matrix(x0)[sel0,-x_sub_sel]
                                         xconst0_sub   <- cbind(xconst0_sub, x0_sub1_const)
                                         if(nrow(coords0)==1){
                                           meig0_sub     <- meigen0(meig_sub, coords0)
                                         } else {
                                           meig0_sub     <- meigen0(meig_sub, coords0[sel0,])
                                         }
                                         pmod0           <- predict0(Mod0,x0=x0_sub,xconst0=xconst0_sub,xgroup0=NULL,
                                                                     meig0=meig0_sub,compute_se = TRUE,weight0=Wei0[sel0,k])
                                         #mod<-Mod0; meig0=meig0_sub; x0 = x0_sub; xconst0 = xconst0_sub; xgroup0 = NULL
                                         #offset0 = NULL; weight0 = Wei0[sel0,k]; compute_se=TRUE
                                         #compute_quantile = FALSE
                                       } else {
                                         pmod0           <- NA
                                       }
                                     }

                                     Mod0$other$x_id<-rep(FALSE,ncol(x_sub2))
                                     Mod0$other$x_id[x_sub_sel]<-TRUE
                                     Mod0$other$x_id[1]<-FALSE

                                     #Mod0$other$x_id <- Mod0$other$x_id[-1]
                                     Mod0$other$xf_id<-Mod0$other$xf_id[xconst_sub_sel]
                                     Mod0$other$x    <-x_sub2
                                     Mod0$other$nx   <-ncol(x_sub2) + length(xconst_sub_sel)####### caution 230925
                                     Mod0$other$nxf  <-length(xconst_sub_sel)
                                     evSqrts_Mod0    <-Mod0$other$evSqrts
                                     b_covs_Mod0     <-Mod0$other$b_covs
                                     b_s_Mod0        <-Mod0$other$b_s
                                     Mod0$other$evSqrts<-list(NULL)
                                     Mod0$other$b_covs<-list(NULL)
                                     Mod0$other$b_s<-list(NULL)

                                     qqq_a <- qqq_b  <- 1
                                     for(ppp in 1:Mod0$other$nx){
                                       if(ppp %in% x_sub_sel){
                                         Mod0$other$evSqrts[[ppp]] <- evSqrts_Mod0[[qqq_a]]
                                         Mod0$other$b_covs[[ppp]]  <- b_covs_Mod0[[qqq_a]]
                                         Mod0$other$b_s[[ppp]]     <- b_s_Mod0[[qqq_a]]

                                         qqq_a <- qqq_a + 1
                                       } else {
                                         Mod0$other$evSqrts[[ppp]] <- 0
                                         Mod0$other$b_covs[[ppp]]  <- Mod0$c$SE[Mod0$other$nxf + qqq_b]^2
                                         Mod0$other$b_s[[ppp]]     <- Mod0$c$Estimate[Mod0$other$nxf + qqq_b]
                                         qqq_b <- qqq_b + 1
                                       }
                                     }

                                   } else {
                                     Mod0       <- mod_sub2
                                     y_sub2     <- y_sub1
                                     x_sub2     <- x_sub1
                                     xconst_sub2<- xconst_sub1
                                     weight_sub2<- weight_sub1
                                     Id         <- Id0
                                   }

                                   if( is.null( x_sub2 ) ){
                                     b_vc_rev  <-data.frame(Mod0$b_vc)
                                     bse_vc_rev<-data.frame(Mod0$bse_vc)

                                     X1	    <- as.matrix( rep(1,nrow(Mod0$b_vc) ) )
                                     names(b_vc_rev)<- names(bse_vc_rev)<-xname	<- "(Intercept)"
                                     Mod0$b_vc  <- b_vc_rev
                                     Mod0$bse_vc<- bse_vc_rev

                                     if( !is.null(coords0) & !is.na(pmod0[1]) ){
                                       b_vc_rev0  <-data.frame(pmod0$b_vc)
                                       bse_vc_rev0<-data.frame(pmod0$bse_vc)
                                       names(b_vc_rev0)<- names(bse_vc_rev0)<-xname	<- "(Intercept)"
                                       pmod0$b_vc  <- b_vc_rev0
                                       pmod0$bse_vc<- bse_vc_rev0
                                     }

                                   } else if( ncol( x_sub2 )!= ncol( Mod0$b_vc ) ){# reshape b_vc and bse_vc in Mod0
                                     b_vc_rev  <- bse_vc_rev<- data.frame(matrix(NA,nrow=nrow(x_sub2),ncol=ncol(x_sub2)))
                                     b_vc_rev[,x_sub_sel]   <- Mod0$b_vc
                                     bse_vc_rev[,x_sub_sel] <- Mod0$bse_vc

                                     if( !is.null(coords0) & !is.na(pmod0[1])){
                                       b_vc_rev0  <- bse_vc_rev0<- data.frame(matrix(NA,nrow=nrow(x0_sub),ncol=ncol(x_sub2)))
                                       b_vc_rev0[,x_sub_sel]   <- pmod0$b_vc
                                       bse_vc_rev0[,x_sub_sel] <- pmod0$bse_vc
                                     }

                                     if( length( xconst_sub_sel ) != nrow(Mod0$c) ){
                                       nxf_sel  <- length(xconst_sub_sel)
                                       for(p in 1:length( x_sub_sel_const ) ){
                                         b_vc_rev[,x_sub_sel_const[p]]  <- Mod0$c$Estimate[nxf_sel+p]
                                         bse_vc_rev[,x_sub_sel_const[p]]<- Mod0$c$SE[nxf_sel+p]
                                         if( !is.null( coords0 ) & !is.na(pmod0[1])){
                                           b_vc_rev0[,x_sub_sel_const[p]]  <- Mod0$c$Estimate[nxf_sel+p]
                                           bse_vc_rev0[,x_sub_sel_const[p]]<- Mod0$c$SE[nxf_sel+p]
                                         }
                                       }
                                       Mod0$c     <- Mod0$c[-(nxf_sel + 1:length( x_sub_sel_const )),]
                                     }

                                     names(b_vc_rev)  <- names(bse_vc_rev)<- c( "(Intercept)", xname )
                                     Mod0$b_vc        <- b_vc_rev
                                     Mod0$bse_vc      <- bse_vc_rev
                                     if( !is.null(coords0) & !is.na(pmod0[1])){
                                       names(b_vc_rev0)  <- names(bse_vc_rev0)<- c( "(Intercept)", xname )
                                       pmod0$b_vc        <- b_vc_rev0
                                       pmod0$bse_vc      <- bse_vc_rev0
                                     }

                                     s_rev            <- data.frame(matrix(c(0,NA),nrow=2,ncol=ncol(x_sub2)))
                                     s_rev[,x_sub_sel]<- Mod0$s
                                     rownames(s_rev)  <- c("random_SD","Moran.I/max(Moran.I)")
                                     names(s_rev)     <- c( "(Intercept)", xname )
                                     Mod0$s           <- s_rev
                                   }

                                   if( !is.null( xconst_sub2 ) ){
                                     if( ncol( xconst_sub2 )!= nrow( Mod0$c ) ){# reshape b_vc and bse_vc in Mod0
                                       c_rev          <- data.frame(matrix(NA,nrow=ncol(xconst_sub2),ncol=4))
                                       c_rev[xconst_sub_sel,]<-Mod0$c
                                       rownames(c_rev)<- xconst_name
                                       names(c_rev)   <- c("Estimate","SD","z_value","p_value")
                                       Mod0$c         <- c_rev
                                     }
                                   }

                                   weight_sc   <- weight_sub2*Wei[Id[[k]],k]*Mod0$other$w_scale
                                   n_sub       <- length(y_sub2)- ifelse(method=="reml",Mod0$other$nx,0)####changed
                                   nw_sub      <- sum(Wei[Id[[k]],k])#
                                   nw_sub2     <- sum(Wei[Id[[k]],k]) - ifelse(method=="reml",Mod0$other$nx,0)
                                   eWe_unsc    <- sum(Wei[Id[[k]],k]*Mod0$resid^2)
                                   eWe_sc      <- sum(weight_sc*Mod0$resid^2)
                                   uu          <- sum(Mod0$other$b_sub^2)
                                   term2_unsc  <- (nw_sub2/2)*log(2*pi*(eWe_unsc+uu)/n_sub)
                                   term2_sc    <- (n_sub/2)*log(2*pi*(eWe_sc+uu)/n_sub)

                                   pen         <- -1/2*sum(log(weight_sc)) + term2_sc - term2_unsc + n_sub/2 - nw_sub/2
                                   loglik      <- Mod0$e$stat[3] + pen
                                   #Mod[[i]]<-list(Mod0,loglik,pmod0)
                                   #}

                                   list(Mod0,loglik,pmod0)
                                 })
  if(parallel) stopCluster( ncores )

  Mod[[cl_num2]] <- Mod_global

  ######################## prediction
  Pred         <- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
  Pred_se      <- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
  for(k in 1:cl_num){ ### this part is cl_num
    Pred   [Id[[k]],k]<-Mod[[k]][1][[1]]$pred$pred
    Pred_se[Id[[k]],k]<-Mod[[k]][1][[1]]$pred$pred_se
  }

  b_g                 <- Mod[[cl_num2]][1][[1]]$b_g
  if( inherits( Mod[[cl_num2]][1][[1]], "besf_vc" )){
    Pred[,cl_num2]    <- Mod[[cl_num2]][1][[1]]$pred
    Pred_se[,cl_num2] <- Mod[[cl_num2]][1][[1]]$e[1,1]
  } else {
    Pred[,cl_num2]      <- Mod[[cl_num2]][1][[1]]$pred$pred
    Pred_se[,cl_num2]   <- Mod[[cl_num2]][1][[1]]$pred$pred_se
  }

  pred_se2<- sqrt(1/rowSums(Wei/Pred_se^2,na.rm=TRUE))
  pred2   <- pred_se2^2*rowSums(Wei*Pred/Pred_se^2,na.rm=TRUE)
  pred3   <- data.frame(pred=pred2,pred_se=pred_se2)
  resid   <- y - pred3$pred

  ######################## error statistics
  loglik <-0
  np     <-0
  for(k in 1:cl_num2){
    loglik<-loglik + Mod[[k]][[2]]
    np     <-np + Mod[[k]][[1]]$other$np
  }
  AIC		  <- -2 * loglik + np * 2
  BIC		  <- -2 * loglik + np * log( n )
  SSE     <- sum( resid ^ 2 )
  SSY		  <- sum( ( y - mean( y ) ) ^ 2 )
  r2_0    <-1 - SSE / SSY

  #r2_0		<- cor(y,pred3$pred)^2#1 - SSE / SSY
  r2		  <- 1 - ( 1- r2_0 ) * ( n - 1 ) / ( n - np - 1 )
  sig_org <- sum( resid ^2 )/( n -nx - nxf )
  e_stat  <- data.frame( stat = c( sqrt( sig_org ), r2, loglik, AIC, BIC ) )

  if( method == "reml" ) lik_nam	<- "rlogLik"
  if( method == "ml"   ) lik_nam	<- "logLik"
  rownames( e_stat ) <- c( "resid_SD", "adjR2(cond)", lik_nam, "AIC", "BIC")

  ######################## prediction at missing sites
  pred3_0<-BETA0<-BETA0_se<-BETA0_z<-BETA0_p<-NULL
  if( !is.null(coords0) ){
    for(k in 1:ncol(Wei)){
      sel0        <- which( Wei0[,k]>0 )
      if( !is.na(Mod[[k]][[3]][1]) ){
        Pred0[sel0,k]   <-Mod[[k]][[3]]$pred$pred
        Pred0_se[sel0,k]<-Mod[[k]][[3]]$pred$pred_se
        for(kk in 1:(nx+1)){
          B_vc0[[kk]][sel0,k]  <- Mod[[k]][[3]]$b_vc[,kk]
          Bse_vc0[[kk]][sel0,k]<- Mod[[k]][[3]]$bse_vc[,kk]
        }
      }
    }

    pred_se2_0<- sqrt(1/rowSums(Wei0/Pred0_se^2,na.rm=TRUE))
    pred2_0   <- pred_se2_0^2*rowSums(Wei0*Pred0/Pred0_se^2,na.rm=TRUE)
    pred3_0   <- data.frame(pred=pred2_0,pred_se=pred_se2_0)

    BETA0     <- BETA0_se <- matrix(NA,nrow=n0,ncol=(nx + 1))
    for(kk in 1:(nx + 1)){
      BETA0  [ ,kk ] <- pred_se2_0^2*rowSums(Wei0*B_vc0[[kk]]/Pred0_se^2,na.rm=TRUE)
      BETA0_se[ ,kk ] <- sqrt(1/rowSums(Wei0/Bse_vc0[[kk]]^2,na.rm=TRUE))
    }

    BETA0       <- as.data.frame(BETA0)
    BETA0_se    <- as.data.frame(BETA0_se)
    BETA0_z     <- as.data.frame(BETA0/BETA0_se)
    names(BETA0)<- names(BETA0_se)<-names(BETA0_z)<-names(mod$b_vc)
    BETA0_p	    <- 2 - 2 * pt( abs( as.matrix( BETA0_z ) ), df = Inf )
  }

  ######################## SVCs
  Beta    <- Beta_se<-matrix(NA,nrow = n, ncol =ncol(Mod[[1]][1][[1]]$b_vc) )
  for(kk in 1:ncol(Mod[[1]][1][[1]]$b_vc) ){
    Beta0   <- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
    Beta0_se<- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
    for(k in 1:cl_num){
      Beta0[Id[[k]],k]   <-Mod[[k]][1][[1]]$b_vc[,kk]
      Beta0_se[Id[[k]],k]<-Mod[[k]][1][[1]]$bse_vc[,kk]
    }
    Beta0[,cl_num2]       <- Mod[[cl_num2]][1][[1]]$b_vc[,kk]
    Beta0_se[,cl_num2]    <- Mod[[cl_num2]][1][[1]]$bse_vc[,kk]

    Wei2 <- Wei
    Wei2[is.na(Beta0)|Beta0==0]<-0
    Wei2 <-Wei2/rowSums(Wei2)
    Beta[,kk]   <- pred_se2^2*rowSums(Wei2*Beta0/Pred_se^2,na.rm=TRUE)
    Beta_se[,kk]<- sqrt(1/rowSums(Wei2/Beta0_se^2,na.rm=TRUE))

    NA_row            <- rowSums(Beta0,na.rm=TRUE)==0
    Beta[NA_row,kk]   <- NA
    Beta_se[NA_row,kk]<- NA
  }

  b_const    <- setdiff(1:ncol(Beta),x_sel_global)
  x_const_avc<- TRUE
  if( x_const_avc==TRUE & length(b_const) > 0 ){
    for(bb in b_const){
      Beta[,bb]   <- mean(Beta[,bb])
      Beta_se[,bb]<- sqrt(mean(Beta_se[,bb]^2))
    }
  }

  Beta       <- as.data.frame(Beta)
  Beta_se    <- as.data.frame(Beta_se)
  Beta_z     <- as.data.frame(Beta/Beta_se)
  Beta_p	   <- 2 - 2 * pt( abs( as.matrix( Beta_z ) ), df = Inf )
  Beta_p_adj <- data.frame( apply(Beta_p, 2, function(x) p.adjust(x, method="BY" ) ))
  names(Beta)<- names(Beta_se)<-names(Beta_z)<-names(Beta_p)<-names(Beta_p_adj)<-names( Mod[[1]][1][[1]]$b_vc )

  ######################## Constant coefficients
  if( !is.null(xconst) ){
    C      <- C_se<-matrix(0,nrow = n, ncol =ncol(xconst))#nrow(Mod[[k]][1][[1]]$c
    for(kk in 1:ncol(xconst)){#nrow(Mod[[k]][1][[1]]$c)
      C0   <- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
      C0_se<- Matrix(nrow = n, ncol = cl_num2, data = 0, sparse = TRUE)
      for(k in 1:cl_num){
        C0[Id[[k]],k]   <-Mod[[k]][1][[1]]$c[kk,1]
        C0_se[Id[[k]],k]<-Mod[[k]][1][[1]]$c[kk,2]
      }
      C0[,cl_num2]     <- Mod[[cl_num2]][1][[1]]$c[kk,1]
      C0_se[,cl_num2]  <- Mod[[cl_num2]][1][[1]]$c[kk,2]

      Wei2              <- Wei
      Wei2[is.na(C0)|C0==0]<- 0
      Wei2              <- Wei2/rowSums(Wei2)
      C[,kk]            <- pred_se2^2*rowSums(Wei2*C0/Pred_se^2,na.rm=TRUE)
      C_se[,kk]         <- sqrt(1/rowSums(Wei2/C0_se^2,na.rm=TRUE))
    }

    c_est <- colMeans(C)
    c_se  <- apply(C_se,2,function(x) sqrt(mean(x^2)))
    c_z   <- c_est/c_se
    c_p   <- pnorm(abs(c_z), lower.tail = FALSE)
    cest  <- data.frame(Estimate=c_est,SE=c_se,z_value=c_z,p_value=c_p)
    row.names(cest)<- names( as.data.frame( xconst ) )
  } else {
    C     <- C_se <-cest <-NULL
  }

  ######################## Variance parameters
  s           <- NULL
  m           <- NULL
  for(k in 1:cl_num){
    s0        <- data.frame(Mod[[k]][[1]]$s[1,])
    s         <- rbind(s,s0)

    m0        <- data.frame(Mod[[k]][[1]]$s[2,])
    m         <- rbind(m,m0)
  }
  #xname       <- names( data.frame(Mod[[k]][[1]]$b_vc ))[-1]
  names(s)    <- names(m)<-names( Beta ) #c( "(Intercept)", xname )
  rownames(s) <- paste0("random_SD (cluster",1:cl_num,")")
  rownames(m) <- paste0("Moran.I/max(Moran.I)(cluster",1:cl_num,")")
  m[is.nan(as.matrix(m))]<-NA
  nan_id      <- is.nan(as.matrix(m))
  sm          <- list(s,m)

  sm_global   <- s_g <-NULL
  sm_global <- data.frame(Mod[[cl_num2]][[1]]$s)
  names(sm_global)<- names( Beta ) #c( "(Intercept)", xname)
  rownames(sm_global)[1:2]<-c("random_SD (global)","Moran.I/max(Moran.I) (global)")
  s_g       <- Mod[[cl_num2]][[1]]$s_g

  ######################## NULL model
  if( is.null( weight ) ){
    mod_NULL_id0<-" )"
  } else {
    if( length(weight)==1 ) weight <-rep( weight, n )
    mod_NULL_id0<-", weights = weight )"
  }

  resf_flag  <- ifelse( is.null(x), 1, 0)
  if( is.null( x ) & is.null( xconst ) ){
    mod_NULL   <- lm(y ~ 1, weights = weight )
    mod_NULL_id<- paste("lm( y ~ 1", mod_NULL_id0, sep="")
  } else if( is.null( x ) ){
    lm_dat     <- data.frame(y, xconst)
    mod_NULL   <- lm(y ~ ., data=lm_dat, weights = weight )
    if( resf_flag ){
      mod_NULL_id<- paste("lm( y ~ x", mod_NULL_id0, sep="")
    } else {
      mod_NULL_id<- paste("lm( y ~ xconst", mod_NULL_id0, sep="")
    }
  } else if( is.null( xconst ) ){
    lm_dat     <- data.frame(y, x)
    mod_NULL   <- lm(y ~ ., data=lm_dat, weights = weight )
    mod_NULL_id<- paste("lm( y ~ x", mod_NULL_id0, sep="")
  } else {
    lm_dat     <- data.frame(y, x, xconst)
    mod_NULL   <- lm(y ~ ., data = lm_dat, weights = weight )
    mod_NULL_id<- paste("lm( y ~ x + xconst", mod_NULL_id0, sep="")
  }

  e_stat_N     <- data.frame( stat = c( logLik(mod_NULL), AIC(mod_NULL), BIC(mod_NULL) ) )
  rownames( e_stat_N ) <- c( lik_nam, "AIC", "BIC")
  e_stat_NULL  <- list(NULL)
  e_stat_NULL[[1]]<- e_stat_N
  e_stat_NULL[[2]]<- mod_NULL_id

  Mod_local        <-list( NULL )
  for(i in 1:cl_num){
    Mod_local[[i]]        <- Mod[[i]][1][[1]]
  }

  Mod_global          <- NULL
  Mod_global        <- Mod[[cl_num2]][1][[1]]

  other   <- list(C=C,C_se,R=R,Wei=Wei,Wei0 = Wei0, Beta0=Beta0,Beta0_se=Beta0_se,
                  method = method, e_NULL = e_stat_NULL,coords=coords, Pred_se=Pred_se, Pred0=Pred0, Pred0_se=Pred0_se,
                  Mod_local=Mod_local, Mod_global=Mod_global, x=x, y=y, dif=dif)

  result  <- list(b_vc=Beta,bse_vc=Beta_se, z_vc=Beta_z, p_vc=Beta_p_adj, p_vc_naive=Beta_p,c=cest,b_g=b_g,
                  s=sm, s_global=sm_global, s_g = s_g, e = e_stat, pred=pred3,resid=resid,cl=cl,
                  pred0 = pred3_0, b_vc0 = BETA0, bse_vc0 = BETA0_se, z_vc0 = BETA0_z, p_vc0 = BETA0_p,
                  other=other, call = match.call() )
  class( result ) <- "addlearn_local"
  return( result )
}

print.addlearn_local <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n----Spatially varying coefficients on x (summary)----\n")

  cat("\nCoefficient estimates:\n")
  print( summary( x$b_vc ) )
  cat("\nStatistical significance (adjusted by Benjamini-Yekutieli method):\n")
  p01<-apply(x$p_vc,2,function(x) sum(x<0.01,na.rm=TRUE))
  p05<-apply(x$p_vc,2,function(x) sum(x<0.05,na.rm=TRUE)) - p01
  p10<-apply(x$p_vc,2,function(x) sum(x<0.10,na.rm=TRUE)) - p01 - p05
  p90<-length(x$p_vc[,1]) - p01 - p05 - p10
  pv <-data.frame(rbind( p90, p10, p05, p01))
  names(pv)[1]  <- "Intercept"
  row.names(pv) <- c("Not significant", "Significant (10% level)",
                     "Significant ( 5% level)","Significant ( 1% level)")
  print(pv)
  if( !is.null(x$c) ){
    cat("\n----Constant coefficients on xconst----------------------------\n")
    print( x$c )
  }

  cat("\n----Variance parameters----------------------------------\n")
  av_s<-data.frame(rbind(colMeans(x$s[[1]],na.rm=TRUE),colMeans(x$s[[2]],na.rm=TRUE)))
  av_s[is.na(av_s)]<-NA
  row.names(av_s)<-c("random_SD","Moran.I/max(Moran.I)")
  names(av_s)[1] <-"(Intercept)"
  cat("\nSpatial effects (Local sub-models; Average):\n")
  print( av_s )

  if(!is.null(x$s_global)){
    cat("\nSpatial effects (Global sub-model):\n")
    row.names(x$s_global)<-c("random_SD","Moran.I/max(Moran.I)")
    print( x$s_global )
  }

  if( !is.null(x$s_g) ){
    cat("\nGroup effects:\n")
    print(x$s_g)
  }

  cat("\n----Error statistics-------------------------------------\n")
  print(x$e)

  loglik_NULL<- x$other$e_NULL[[1]][1,1]
  AIC_NULL   <- x$other$e_NULL[[1]][2,1]
  BIC_NULL   <- x$other$e_NULL[[1]][3,1]
  mod_NULL   <- x$other$e_NULL[[2]]
  ml_name    <- ifelse( x$other$method=="reml", "(r)loglik: ", "loglik: " )
  cat( paste('\nNULL model: ', mod_NULL, sep="") )
  cat( paste("\n   ",ml_name,format(loglik_NULL,digits=7),sep=""))
  cat( paste(" ( AIC: ",
             format(AIC_NULL,digits=7), ",  BIC: ", format(BIC_NULL,digits=7)," )\n",sep=""))

  invisible(x)
}
