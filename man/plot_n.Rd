\name{plot_n}
\alias{plot_n}
\title{Plot non-spatially varying coefficients (NVCs)}
\usage{
plot_n( mod, xnum = 1, xtype = "x", cex.lab = 20,
        cex.axis = 15, lwd = 1.5, ylim = NULL, nmax = 20000 )
}
\description{
This function plots non-spatially varying coefficients (NVCs; coefficients varying with respect to explanatory variable value) and their 95 percent confidence intervals
}
\arguments{
\item{mod}{Outpot from \code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, or  \code{\link{besf_vc}} function}
\item{xnum}{The NVC on the xnum-th explanatory variable is plotted. Default is 1}
\item{xtype}{Effective for \code{\link{resf_vc}} and \code{\link{besf_vc}}.  If "x", the num-th NVC in the spatially and non-spatially varying coefficients on x is plotted. If "xconst", the num-th NVC on xconst is plotted. Default is "x"}
\item{cex.lab}{The size of the x and y axis labels}
\item{cex.axis}{The size of the tick label numbers}
\item{lwd}{The width of the line drawing the coefficient estimates}
\item{ylim}{The limints of the y-axis}
\item{nmax}{If sample size exceeds nmax, nmax samples are randomly selected and plotted. Default is 20,000}
}
\seealso{
\code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, \code{\link{besf_vc}}
}
