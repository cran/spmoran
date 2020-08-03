
plot_s   <- function( mod, xnum = 0, btype = "snvc", xtype = "x", pmax = NULL, ncol = 8, col = NULL,
                      inv =FALSE, brks = "regular", cex = 1, nmax = 20000 ){

  pmax0      <- pmax
  xnum       <- xnum + 1
  if( (class(mod) == "resf_vc")|(class(mod) == "besf_vc") ){

    if( xtype == "xconst" ){
      xnum   <- xnum - 1
      dd     <- data.frame(coords = mod$other$coords,b_vc = mod$c_vc, p_vc = mod$cp_vc )
      xname    <-names( mod$c_vc )
      nx       <-length( xname )
      names(dd)<-c("px","py",xname,paste(xname,"_p",sep=""))
    } else {
      if( btype == "snvc" ){
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$b_vc, p_vc = mod$p_vc )
      } else if( btype == "svc" ) {
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$B_vc_s[[1]], p_vc = mod$B_vc_s[[4]] )
      } else if( btype == "nvc" ) {
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$B_vc_n[[1]], p_vc = mod$B_vc_n[[4]] )

      }
      xname    <-names( mod$b_vc )
      xname[ 1 ]<-"Spatially.dependent.intercept"
      nx       <-length( xname )
      names(dd)<-c("px","py",xname,paste(xname,"_p",sep=""))
    }

  } else if( (class(mod) == "esf")|(class(mod) == "resf")|(class(mod) == "besf") ){
    if( ( (xnum >1)|(btype=="nvc") )&( !is.null( mod$c_vc[1,1] ) ) ){
      dd     <- data.frame(coords = mod$other$coords,b_vc = mod$c_vc, p_vc = mod$cp_vc )
      xnum   <- xnum - 1
      xname    <-names( mod$c_vc )
      nx       <-length( xname )
      names(dd)<-c("px","py",xname,paste(xname,"_p",sep=""))
    } else {#
      dd       <-data.frame(coords = mod$other$coords,sf = mod$sf )
      names(dd)<-c("px","py","Spatially.depepdent.component" )
      xnum     <-1
      nx       <-0
    }
  }

  if( dim( dd )[ 1 ] > nmax ){
    samp     <- sample(dim( dd )[ 1 ], nmax)
    dd       <- dd[samp,]
  }

  xlim   <-range( dd$px )
  ylim   <-range( dd$py )
  xlim[1]<-xlim[1] - abs(xlim[2]-xlim[1])/20
  xlim[2]<-xlim[2] + abs(xlim[2]-xlim[1])/20
  ylim[1]<-ylim[1] - abs(ylim[2]-ylim[1])/20
  ylim[2]<-ylim[2] + abs(ylim[2]-ylim[1])/20
  coordinates(dd) <- c("px","py")

  if( brks == "quantile" ){
    cols <- quantile(dd@data[, xnum ],probs=seq(0,1,len=ncol+1))
  } else if( brks == "regular" ){
    cols <- seq(min(dd@data[, xnum ]),max(dd@data[, xnum ]),len=ncol+1)
  }

  if( is.null(col)){
    palette <- get_col_regions()
  } else {
    palette <- brewer.pal(n = ncol, name = col)
  }

  if(inv) palette <- rev(palette)
  if( is.null(pmax) ) pmax <- Inf
  if( xnum == 1 ) pmax <- Inf

  SE_rat      <- sd(dd@data[,xnum]) / abs( mean(dd@data[,xnum]) )
  if( SE_rat  < 0.0001 ){
    #coll      <- palette[ round( length( palette )/2 ) ]
    #pp   <-spplot( dd[dd@data[,nx + xnum] < pmax,], xnum, colorkey=TRUE, cex=cex, col.regions = coll, #palette,
    #        col = "transparent",xlim=xlim, ylim=ylim,main=names(dd)[xnum])
    #message( "Note: The coefficients are uniform over space" )
    if( SE_rat == 0){
      message("Note: Coefficients are not plotted because they are constant")
    } else {
      message("Note: Coefficients are not plotted because they are approximately constant")
    }
    #pp <-NULL

  } else {
    test<-try(pp  <-spplot( dd[dd@data[,nx + xnum] < pmax,], xnum, colorkey=TRUE, cex=cex, col.regions = palette,
            cuts = cols, col = "transparent",xlim=xlim, ylim=ylim,main=names(dd)[xnum]),silent=TRUE)
    if( class(test)=="try-error" ){
      names(dd)[xnum]<-paste( "x", xnum, sep = "" )
      pp  <-spplot( dd[dd@data[,nx + xnum] < pmax,], xnum, colorkey=TRUE, cex=cex, col.regions = palette,
                              cuts = cols, col = "transparent",xlim=xlim, ylim=ylim,main=names(dd)[xnum])
      message( "Note: Covariate name is changed to aviod error caused by special characters" )
    }

    if( sum(dd@data[,nx + xnum] < pmax)==0 ){
      message( "Note: Estimates are not plotted because all the p-values are above pmax" )
    }

    if( ( xnum ==1 )&( !is.null( pmax0 ) ) ){
      message( "Note: pmax is not supported for spatially varying intercept" )
    }
    pp
  }
}
