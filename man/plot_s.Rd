\name{plot_s}
\alias{plot_s}
\title{Mapping spatially and spatio-temporally varying coefficients}
\usage{
plot_s( mod, xnum = 0, btype = "all", xtype = "x", pmax = NULL, ncol = 8,
        col = NULL, inv =FALSE, brks = "regular", cex = 1, pch = 20,
        nmax = 20000, coords_z1_lim=NULL, coords_z2_lim = NULL)
}
\description{
This function plots spatially varying coefficients (SVC) and spatio-temporally varying coefficients (STVC) with/without coefficient varying with respect to the value of the explanatory variable (NVC). Namely, the full varying coefficient equals STVC + NVC.
}
\arguments{
\item{mod}{Outpot from \code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, \code{\link{besf_vc}}, or \code{\link{addlearn_local}} function}
\item{xnum}{For \code{\link{resf_vc}}, \code{\link{besf_vc}}, and  \code{\link{addlearn_local}}, xnum-th SVC/STVC is plotted. If num = 0, varying intercept is plotted. For \code{\link{resf}} and \code{\link{besf}}, estimated spatially dependent residual process is plotted irrespective of the xnum value. Default is 0}
\item{btype}{Effective if x_nvc =TRUE in \code{\link{resf_vc}} and \code{\link{besf_vc}}. If "all" (default), the estiamted varying coefficients (S(T)VC + NVC) are plotted as they are. If "svc", S(T)VC is plotted. If "nvc", NVC is plotted.}
\item{xtype}{If "x" (default), coefficients on x is plotted. If "xconst", those on xconst is plotted}
\item{pmax}{The maximum p-value for the varying coefficients to be displayed. For example, if pmax = 0.05, the only coefficients that are statistically significant at the 5 percent level are plotted. If NULL, all the coefficients are plotted. Default is NULL}
\item{ncol}{Number of colors in the color palette. Default is 8}
\item{col}{Color palette used for the mapping. If NULL, the blue-pink-yellow color scheme is used. Palettes in the RColorBrewer package are available. Default is NULL}
\item{inv}{If TRUE, the color palett is inverted. Default is FALSE}
\item{brks}{If "regular", color is changed at regular intervals. If "quantile", color is changed for each quantile}
\item{cex}{Size of the dots representing sample sites}
\item{pch}{A number indicating the symbol to use}
\item{nmax}{If sample size exceeds nmax, nmax samples are randomly selected and plotted. Default is 20,000}
\item{coords_z1_lim}{Value range for coords_z[,1] in the \code{\link{meigen}}/\code{\link{meigen_f}} function (vector). If is has two elements, the samples whose coords_z[,1] values are in between these values are plotted. If it is a scalar, samples satisfying coords_z[,1]==coords_z1_lim is plotted}
\item{coords_z2_lim}{Value range for coords_z[,2] (vector).}
}
\seealso{
\code{\link{resf}}, \code{\link{besf}}, \code{\link{resf_vc}}, \code{\link{besf_vc}}, \code{\link{addlearn_local}}
}
