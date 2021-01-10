\name{lslm}
\alias{lslm}
\title{Low rank spatial lag model (LSLM) estimation}
\usage{
lslm( y, x, weig, method = "reml", boot = FALSE, iter = 200 )
}
\description{
This function estimates the low rank spatial lag model.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K)}
\item{weig}{eigenvectors and eigenvalues of a spatial weight matrix. Output from \code{\link{weigen}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{boot}{If it is TRUE, confidence intervals for the spatial dependence parameters (s), the mean direct effects (de), and the mean indirect effects (ie), are estimated through a parametric bootstrapping. Default is FALSE}
\item{iter}{The number of bootstrap replicates. Default is 200}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{s}{Vector of estimated shrinkage parameters (2 x 1). The first and the second elements denote the estimated rho parameter (sp_rho) quantfying the scale of spatial dependence, and the standard error of the spatial dependent component (sp_SE), respectively. If boot = TRUE, their 95 percent confidence intervals and the resulting p-values are also provided}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{de}{Matrix with columns for the estimated mean direct effects on x. If boot = TRUE, their 95 percent confidence intervals and the resulting p-values are also provided}
\item{ie}{Matrix with columns for the estimated mean indirect effects on x. If boot = TRUE, their 95 percent confidence intervals and the resulting p-values are also provided}
\item{r}{Vector of estimated random coefficients on the spatial eigenvectors (L x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outputs, which are internally used}
}
\references{
Murakami, D., Seya, H. and Griffith, D.A. (2018) Low rank spatial econometric models. Arxiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{weigen}}, \code{\link{lsem}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords  <- boston.c[,c("LON", "LAT")]
weig 	  <- weigen(coords)
res	    <- lslm(y=y,x=x,weig=weig)
## res  <- lslm(y=y,x=x,weig=weig, boot=TRUE)
res

}

