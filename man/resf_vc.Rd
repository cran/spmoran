\name{resf_vc}
\alias{resf_vc}
\title{Spatially varying coefficient modeling with/without automatic coefficient selection}
\usage{
resf_vc(y, x=NULL, xgroup = NULL, xconst=NULL, meig, method="reml",
        penalty="bic", allsvc=FALSE, maxiter=30, sizelimit=2000 )
}
\description{
This function estimates a spatially varying coefficient model based on the random effects eigenvector spatial filtering (RE-ESF) approach. Spatially varying coefficients are selected to stablize the estimates.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables whose coefficients are allowed to vary across geographical space (N x K). Default is NULL}
\item{xgroup}{Matrix of group indexes. The indeces may be group IDs (numbers) or group names (N x K_{group}). Default is NULL}
\item{xconst}{Matrix of explanatory variables whose coefficients are assumed constant across space (N x K_{const}). Default is NULL}
\item{meig}{Moran's eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select varying coefficients and stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion  (2K) (see Muller et al., 2013). Default is "bic"}
\item{allsvc}{If it is TRUE, SVCs are selected using the penalty. If it is FALSE, all the coefficients on the explanatory variables in x are SVCs. Default is FALSE}
\item{maxiter}{Maximum number of iterations. Default is 30}
\item{sizelimit}{Maximum size of matrix being inverted. Roughly speaking, this value contraints K + KL, where K and L are the numbers of the explanatory variables and eigen-pairs, respectively. If K + KL exceeds sizelimit, L is reduced to fullfill the contraint. Default is 2000}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on xconst, their standard errors, t-values, and p-values (K_{const} x 4)}
\item{b_g}{List of K_{group} matrices with columns for the estimated group effects, their standard errors, and t-values}
\item{s}{Matrix of estimated variance parameters (2 x K). The (1, k)-th element denotes the standard error of the k-th SVCs, while the (2, k)-th element denotes the Moran's I value of the SVCs. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{s_g}{Vector of estimated standard errors of the group effects}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{b_vc}{Matrix of estimated spatially varying coefficients (SVCs) on x (N x K)}
\item{bse_vc}{Matrix of estimated standard errors for the SVCs (N x k)}
\item{t_vc}{Matrix of estimated t-values for the SVCs (N x K)}
\item{p_vc}{Matrix of estimated p-values for the SVCs (N x K)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{vc}{Vector indicating whether spatial variations are found or not from each coefficients (K x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Muller, S., Scealy, J.L., and Welsh, A.H. (2013) Model selection in linear mixed models. Statistical Science, 28 (2), 136-167.

Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.

Murakami, D., and Griffith, D.A. (2019). Spatially varying coefficient modeling for large datasets: Eliminating N from spatial regressions. Spatial Statistics, 30, 39-64.

Griffith, D. A. (2003). Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{resf}},\code{\link{besf_vc}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV"]
x       <- boston.c[,c("ZN", "LSTAT")]
xconst  <- boston.c[,c("CRIM", "NOX", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "RM")]
xgroup  <- boston.c[,"TOWN"]
coords  <- boston.c[,c("LAT","LON")]
meig 	  <- meigen(coords=coords)
# meig_f 	<- meigen_f(coords=coords)  ## for fast computation


######## Without the penalty-based SVC selection (just like the usual GWR)
res	    <- resf_vc(y=y,x=x,xconst=xconst,meig=meig,allsvc=TRUE)

######## With the penalty-based SVC selection
# res2	<- resf_vc(y=y,x=x,xconst=xconst,meig=meig)

######## SVC modeling with group effects
# res3	<- resf_vc(y=y,x=x,xconst=xconst,xgroup=xgroup,meig=meig)

res$b
res$s
res$e

res$vc
res$b_vc[1:10,]
res$bse_vc[1:10,]
res$t_vc[1:10,]
res$p_vc[1:10,]

###Plot the first SVC (INDUS)
#
#require(ggplot2)
#ggplot(mapping = aes(x = coords$LON, y = coords$LAT)) +
#  geom_point(aes(colour = res$b_vc[,2])) +
#  scale_color_gradientn(colours=rev(rainbow(4)))
}
