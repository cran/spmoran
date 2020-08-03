\name{plot_s}
\alias{plot_s}
\title{Mapping spatially (and non-spatially) varying coefficients (SVCs or SNVC)}
\usage{
plot_s( mod, xnum = 0, btype = "snvc", xtype = "x", pmax = NULL, ncol = 8,
        col = NULL, inv =FALSE, brks = "regular", cex = 1, nmax = 20000)
}
\description{
This function plots spatially and non-spatially varying coefficients (SNVC) or spatially varying coefficients (SVC). Note that SNVC = SVC + NVC (NVC is a coefficient varying with respect to explanatory variable value)
}
\arguments{
\item{mod}{Outpot from \code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, or  \code{\link{besf_vc}} function}
\item{xnum}{For \code{\link{resf_vc}} and \code{\link{besf_vc}}, xnum-th S(N)VC on x is plotted. If num = 0, spatially varying intercept is plotted. For \code{\link{resf}} and \code{\link{besf}}, estimated spatially dependent component in the residuals is plotted irrespective of the xnum value. Default is 0}
\item{btype}{Effective for \code{\link{resf_vc}} and \code{\link{besf_vc}}. If "snvc" (default), SNVC (= SVC + NVC) is plotted. If "svc" , SVC is plotted. If "nvc", NVC is plotted}
\item{xtype}{If "x" (default), coefficients on x is plotted. If "xconst", those on xconst is plotted}
\item{pmax}{The maximum p-value for the S(N)VC to be displayed. For example, if pmax = 0.05, only coefficients that are statistically significant at the 5 percent level are plotted. If NULL, all the coefficients are plotted. Default is NULL}
\item{ncol}{Number of colors in the color palette. Default is 8}
\item{col}{Color palette used for the mapping. If NULL, the blue-pink-yellow color scheme is used. Palettes in the RColorBrewer package are available. Default is NULL}
\item{inv}{If TRUE, the color palett is inverted. Default is FALSE}
\item{brks}{If "regular", color is changed at regular intervals. If "quantile", color is changed for each quantile}
\item{cex}{Size of the dots representing sample sites}
\item{nmax}{If sample size exceeds nmax, nmax samples are randomly selected and plotted. Default is 20,000}
}
\seealso{
\code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, \code{\link{besf_vc}}
}
