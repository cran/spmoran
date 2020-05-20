\name{lsem}
\alias{lsem}
\title{Low rank spatial error model (LSEM) estimation}
\usage{
lsem( y, x, weig, method = "reml" )
}
\description{
This function estimates the low rank spatial error model.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K)}
\item{weig}{eigenvectors and eigenvalues of a spatial weight matrix. Output from \code{\link{weigen}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements denote the estimated rho parameter (sp_lambda) quantfying the scale of spatial dependent process, and the standard error of the process (sp_SE), respectively.}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{r}{Vector of estimated random coefficients on the spatial eigenvectors (L x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
}
\references{
Murakami, D., Seya, H. and Griffith, D.A. (2018) Low rank spatial econometric models. Arxiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords<- boston.c[,c("LAT","LON")]
weig 	<- weigen( coords )
res	  <- lsem(y=y,x=x,weig=weig)
res
}

