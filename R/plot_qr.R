
plot_qr <- function( mod, pnum = 1, par = "b", cex.main = 22, cex.lab = 20, cex.axis = 18, lwd = 1.5 ){
  if( is.null( mod$B ) ) {
    ci <- FALSE
  } else {
    ci <- TRUE
  }
  xx    <- eval(parse(text=paste("mod$", par, "[ pnum, ]", sep="")))
  xname <- rownames(xx)

  if( ci == FALSE ){
    d          <- data.frame(mod$tau, t( xx ) )
    names( d ) <- c( "Quantile", "Estimates" )
    p          <-ggplot(d, aes_string(x="Quantile", y="Estimates"))
  } else if( ci == TRUE ) {
    d        <- NULL
    for( k in 1:length( mod$tau )){
      xx_b0<- paste( "mod$", toupper( par ), "[[ k ]][ pnum, ]", sep="" )
      xx_b <- eval( parse( text = xx_b0 ) )
      d    <- rbind( d, xx_b )
    }
    d$Quantile <- mod$tau
    p          <-ggplot(d, aes_string(x="Quantile", y="Estimates")) +
      geom_ribbon( aes_string( ymin = "CI_lower", ymax = "CI_upper" ), alpha = 0.2 )
  }
    p      <- p + geom_line(size = lwd ) +
    geom_abline(intercept=0,slope=0,size=1,linetype="dashed") +
    labs(title = xname ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(title       =element_text(size=cex.main)) +
    theme(axis.title.x=element_text(size=cex.lab) ,axis.title.y=element_text(size=cex.lab)) +
    theme(axis.text.x =element_text(size=cex.axis ),axis.text.y=element_text(size=cex.axis))
    if( ( pnum ==2 ) & ( par == "s" ) ){
      p    <- p + ylim(0,1)
    }
  plot(p)
}
