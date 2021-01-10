\name{esf}
\alias{esf}
\title{Spatial regression with eigenvector spatial filtering}
\usage{
esf( y, x = NULL, vif = NULL, meig, fn = "r2" )
}
\description{
This function estimates the linear eigenvector spatial filtering (ESF) model. The eigenvectors are selected by a forward stepwise method.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{vif}{Maximum acceptable value of the variance inflation factor (VIF) (scalar). For example, if vif = 10, eigenvectors are selected so that the maximum VIF value among explanatory variables and eigenvectors is equal to or less than 10. Default is NULL}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{fn}{Objective function for the stepwise eigenvector selection. The adjusted R2 ("r2"), AIC ("aic"), or BIC ("bic") are available. Alternatively, all the eigenvectors in meig are use if fn = "all". This is acceptable for large samples (see Murakami and Griffith, 2019). Default is "r2"}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{s}{Vector of statistics for the estimated spatial component (2 x 1). The first element is the standard error and the second element is the Moran's I value of the estimated spatially dependent component. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{r}{Matrix with columns for the estimated coefficients on Moran's eigenvectors, their standard errors, t-values, and p-values (L x 4)}
\item{vif}{Vector of variance inflation factors of the explanatory variables (N x 1)}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted R2 (adjR2), log-likelihood (logLik), AIC, and BIC}
\item{sf}{Vector of estimated spatial dependent component (E\eqn{\gamma}) (N x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outputs, which are internally used}
}
\references{
Griffith, D. A. (2003). Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Tiefelsdorf, M., and Griffith, D. A. (2007). Semiparametric filtering of spatial autocorrelation: the eigenvector approach. Environment and Planning A, 39 (5), 1193-1221.

Murakami, D. and Griffith, D.A. (2019) Eigenvector spatial filtering for large data sets: fixed and random effects approaches. Geographical Analysis, 51 (1), 23-49.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{resf}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE")]
coords  <- boston.c[,c("LON", "LAT")]

#########Distance-based ESF
meig 	<- meigen(coords=coords)
esfD	<- esf(y=y,x=x,meig=meig, vif=5)
esfD

#########Fast approximation
meig_f<- meigen_f(coords=coords)
esfD	<- esf(y=y,x=x,meig=meig_f, vif=10, fn="all")
esfD

############################Not run
#########Topoligy-based ESF (it is commonly used in regional science)
#
#cknn	<- knearneigh(coordinates(coords), k=4) #4-nearest neighbors
#cmat	<- nb2mat(knn2nb(cknn), style="B")
#meig <- meigen(cmat=cmat, threshold=0.25)
#esfT	<- esf(y=y,x=x,meig=meig)
#esfT
}
