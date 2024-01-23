
plot_n <- function( mod, xnum = 1, xtype = "x", cex.lab = 20,
                    cex.axis = 15, lwd = 1.5, ylim = NULL, nmax = 20000 ){

  #if( (class(mod) == "resf_vc")|(class(mod) == "besf_vc") ){
  if( inherits( mod, "resf_vc")|inherits( mod, "besf_vc" ) ){

    if( ( xtype == "x" )&!is.null(mod$c_vc)&is.null(mod$B_vc_n)){
      xtype <- "xconst"
      message('Note: xtype = "x" is assumed because no varying coefficients were estimated on xconst')
    }

    if( !is.null( mod$other$xconst )&( xtype == "xconst" )&is.null(mod$c_vc) ){
      stop( "Varying coefficients on xconst were not estimated. Specify xconst_nvc=TRUE when using resf_vc/besf_vc" )
    }

    if( xtype =="x" ){
      bvc    <- mod$B_vc_n
      xnum_b <- xnum+1
      xnum_x <- xnum
      Covariate<-as.matrix( mod$other$x)[, xnum_x]

    } else if( xtype == "xconst" ){
      bvc    <- list( mod$c_vc, mod$cse_vc )
      xnum_b <- xnum
      xnum_x <- xnum
      Covariate<-as.matrix( mod$other$xconst )[, xnum_x]

    }
  #} else if( (class(mod) == "resf")|(class(mod) == "besf") ){
  } else if( inherits( mod, "resf")|inherits( mod, "besf") ){
    bvc   <- list( mod$c_vc, mod$cse_vc )
    xnum_b <- xnum
    xnum_x <- xnum+1
    Covariate<-mod$other$x[, xnum_x]
  }

  if( !is.na(bvc[[1]][1,1]) ){
    xname    <-names(bvc[[1]])[ xnum_b ]
    Coefficients<-bvc[[1]][ ,xnum_b ]
    bse      <-bvc[[2]][ ,xnum_b ]
    CI_upper<-Coefficients + 1.96*bse
    CI_lower<-Coefficients - 1.96*bse
    d       <-data.frame(Covariate, Coefficients,CI_upper, CI_lower)
    names(d)[1]<-xname

    if(dim(d)[1] > nmax){
      samp <- c( which.max(d$Covariate), which.min(d$Covariate), sample(dim(d)[1], nmax))
      d    <- d[samp,]
    }

    p   <-ggplot(d, aes_string(x=xname, y="Coefficients")) +
        geom_ribbon( aes_string( ymin = "CI_lower", ymax = "CI_upper" ), alpha = 0.2 )

    p   <- p + geom_line(size = lwd ) +
        geom_abline(intercept=0,slope=0,size=1,linetype="dashed") +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.title.x=element_text(size=cex.lab) ,axis.title.y=element_text(size=cex.lab)) +
        theme(axis.text.x =element_text(size=cex.axis ),axis.text.y=element_text(size=cex.axis))
  } else {
    stop( "Coefficients are missing. Check xtype and xnum" )
  }

  if( !is.null( ylim ) ) p    <-p + ylim(ylim)
  plot(p)
}
