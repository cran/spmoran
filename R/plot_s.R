
plot_s   <- function( mod, xnum = 0, btype = "all", xtype = "x", pmax = NULL, ncol = 8, col = NULL,
                      inv =FALSE, brks = "regular", cex = 1, pch = 20, nmax = 20000,
                      coords_z1_lim=NULL, coords_z2_lim = NULL){

  pmax0        <- pmax
  xnum         <- xnum + 1
  #if( (class(mod) == "resf_vc")|(class(mod) == "besf_vc") ){

  if( inherits(mod, "addlearn_local") ){
    dd     <- data.frame(coords = mod$other$coords,b_vc = mod$b_vc, p_vc = mod$p_vc )
    xname     <- names( mod$b_vc )
    xname[ 1 ]<- "Spatially.dependent.intercept"
    nx        <- length( xname )
    names(dd) <- c("px","py",xname,paste(xname,"_p",sep=""))
  } else if( inherits(mod, "resf_vc")|inherits(mod, "besf_vc") ){

    if( xtype == "xconst" ){
      xnum     <- xnum - 1
      dd       <- data.frame(coords = mod$other$coords,b_vc = mod$c_vc, p_vc = mod$cp_vc )
      xname    <- names( mod$c_vc )
      nx       <- length( xname )
      names(dd)<- c("px","py",xname,paste(xname,"_p",sep=""))
    } else {
      if( btype == "all" ){
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$b_vc, p_vc = mod$p_vc )
      } else if( btype == "svc" ) {
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$B_vc_s[[1]], p_vc = mod$B_vc_s[[4]] )
      } else if( btype == "nvc" ) {
        dd     <- data.frame(coords = mod$other$coords,b_vc = mod$B_vc_n[[1]], p_vc = mod$B_vc_n[[4]] )

      }
      xname     <- names( mod$b_vc )
      xname[ 1 ]<- "Spatially.dependent.intercept"
      nx        <- length( xname )
      names(dd) <- c("px","py",xname,paste(xname,"_p",sep=""))
    }

  } else if( inherits(mod, "esf" )|inherits(mod, "resf" )|inherits( mod, "besf" ) ){
    if( ( (xnum >1)|(btype=="nvc") )&( !is.null( mod$c_vc[1,1] ) ) ){
      dd       <- data.frame(coords = mod$other$coords,b_vc = mod$c_vc, p_vc = mod$cp_vc )
      xnum     <- xnum - 1
      xname    <- names( mod$c_vc )
      nx       <- length( xname )
      names(dd)<- c("px","py",xname,paste(xname,"_p",sep=""))
    } else {
      dd       <- data.frame(coords = mod$other$coords,sf = mod$sf )
      names(dd)<- c("px","py","Spatially.depepdent.component" )
      xnum     <- 1
      nx       <- 0
    }
  }

  if( !is.null(mod$other$coords_z) ){
    nz         <- dim(mod$other$coords_z)[2]

    if( !is.null(coords_z1_lim) ){
      if( length(coords_z1_lim)==1 ){
        sel_coords_z1 <- mod$other$coords_z[,1]==coords_z1_lim
      } else if(length(coords_z1_lim)==2){
        sel_coords_z1 <- mod$other$coords_z[,1]>=min(coords_z1_lim) & mod$other$coords_z[,1]<=max(coords_z1_lim)
      }
      if(sum(sel_coords_z1)==0) stop("No sample in the value range: coords_z1_lim")
      dd<-dd[sel_coords_z1,]
    }

    if( !is.null(coords_z2_lim) & nz==2 ){

      coords_z2<-mod$other$coords_z[,2]
      if( !is.null(coords_z1_lim) ) coords_z2<-coords_z2[mod$other$coords_z[,1]==coords_z1_lim]

      if( length(coords_z2_lim)==1 ){
        sel_coords_z2 <- coords_z2==coords_z2_lim
      } else if(length(coords_z2_lim)==2){
        sel_coords_z2 <- coords_z2>=min(coords_z2_lim) & coords_z2<=max(coords_z2_lim)
      }
      if(sum(sel_coords_z2)==0) stop("No sample in the value range: coords_z1_lim and coords_z2_lim")
      dd<-dd[sel_coords_z2,]
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
  #coordinates(dd) <- c("px","py")
  dd_s   <- st_as_sf(dd, coords=c("px","py") )

  if( brks == "quantile" ){
    #cols <- quantile(dd@data[, xnum ],probs=seq(0,1,len=ncol+1))
    cols  <- quantile(dd[, xnum+2 ],probs=seq(0,1,len=ncol+1))
  } else if( brks == "regular" ){
    #cols <- seq(min(dd@data[, xnum ]),max(dd@data[, xnum ]),len=ncol+1)
    cols  <- seq(min(dd[, xnum+2 ]),max(dd[, xnum+2 ]),len=ncol+1)
  }

  if( is.null(col)){
    palette <- NULL
  } else {
    palette <- brewer.pal(n = ncol, name = col)
  }

  if(inv) palette <- rev(palette)
  if( is.null(pmax) ) pmax <- Inf
  if( xnum == 1 ) pmax <- Inf

  SE_rat      <- sd(dd[,xnum+2]) / abs( mean(dd[,xnum+2]) )
  if( SE_rat  < 0.0001 ){
    if( SE_rat == 0){
      message("Note: Coefficients are not plotted because they are constant")
    } else {
      message("Note: Coefficients are not plotted because they are approximately constant")
    }

  } else {

    test<-try(pp  <-plot( dd_s[dd[,nx + xnum+2] < pmax, xnum], cex=cex, pal = palette, key.pos = 4, pch = pch,
                          breaks = cols, xlim=xlim, ylim=ylim, axes=TRUE), silent=TRUE )
    #test<-try(pp  <-spplot( dd[dd[,nx + xnum+2] < pmax,], xnum, colorkey=TRUE, cex=cex, col.regions = palette,
    #                        cuts = cols, col = "transparent",xlim=xlim, ylim=ylim,main=names(dd)[xnum]),silent=TRUE)
    if( inherits( test, "try-error" ) ){#class(test)=="try-error"
      names(dd)[xnum+2]<-paste( "x", xnum+2, sep = "" )
      #pp  <-spplot( dd[dd[,nx + xnum+2] < pmax,], xnum, colorkey=TRUE, cex=cex, col.regions = palette,
      #              cuts = cols, col = "transparent",xlim=xlim, ylim=ylim,main=names(dd)[xnum])
      pp  <-plot( dd_s[dd[,nx + xnum+2] < pmax, xnum], cex=cex, pal = palette, key.pos = 4, pch = pch,
                  breaks = cols, xlim=xlim, ylim=ylim, axes=TRUE)
      message( "Note: Covariate name is changed to aviod error due to special characters" )
    }

    if( sum(dd[,nx + xnum + 2] < pmax)==0 ){
      message( "Note: Estimates are not plotted because all the p-values are above pmax" )
    }

    if( ( xnum ==1 )&( !is.null( pmax0 ) ) ){
      message( "Note: pmax is not supported for spatially varying intercept" )
    }
    pp
  }
  invisible()
}
