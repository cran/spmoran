\name{plot_qr}
\alias{plot_qr}
\title{Plot quantile regression coefficients estimated from SF-UQR}
\usage{
plot_qr( mod, pnum = 1, par = "b", cex.main = 20, cex.lab = 18, cex.axis = 15, lwd = 1.5 )
}
\description{
This function plots regression coefficients estimated from the spatial filter unconditional quantile regression (SF-UQR) model.
}
\arguments{
\item{mod}{Outpot from the \code{\link{resf_qr}} function}
\item{pnum}{A number specifying the parameter being plotted. If par = "b", the coefficients on the pnum-th explanatory variable are plotted (intercepts are plotted if pnum = 1). If par = "s" and pnum = 1, the estimated standard errors for the reidual spatial process are plotted. If par = "s" and pnum = 2, the Moran's I values of the residual spatial process are plotted. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{par}{If it is "b", regression coefficeints are plotted. If it is "s", shrinkage (variance) parameters for the residual spatial process are plotted. Default is "b"}
\item{cex.main}{Graphical parameter specifying the size of the main title}
\item{cex.lab}{Graphical parameter specifying the size of the x and y axis labels}
\item{cex.axis}{Graphical parameter specifying the size of the tick label numbers}
\item{lwd}{Graphical parameters specifying the width of the line drawing the coefficient estimates}
}
\note{
See \code{\link{par}} for the graphical parameters
}
\seealso{
\code{\link{resf_qr}}
}
