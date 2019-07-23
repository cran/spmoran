\name{besf}
\alias{besf}
\title{Memory-free implementation of the RE-ESF model for very large datasets}
\usage{
besf( y, x = NULL, coords,  method = "reml", covmodel="exp",
      enum = 200, bsize = 3000, cl=NULL)
}
\description{
This function estimates the RE-ESF model through a parallel computation whose memory consumption is independent of the sample size.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{covmodel}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel}
\item{enum}{Number of eigenvectors and eigenvalues to be extracted (scalar). Default is 200}
\item{bsize}{Block/badge size. bsize x bsize elements are iteratively processed during the parallelized computation. Default is 3000}
\item{cl}{Number of cores used for the parallel computation. If cl=NULL, which is the default, the number of available cores is detected and used}
}

\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements denote the standard error and the Moran's I value of the estimated spatially dependent component, respectively. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{r}{Vector of estimated random coefficients on Moran's eigenvectors (L x 1)}
\item{sf}{Vector of estimated spatial dependent component (N x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Griffith, D. A. (2003). Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

Murakami, D. and Griffith, D.A. (2019) A memory-free spatial additive mixed modeling for big spatial data. ArXiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{besf}}, \code{\link{meigen}}, \code{\link{meigen_f}}
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
xgroup  <- boston.c[,"TOWN"]
coords  <- boston.c[,c("LAT","LON")]
# res	  <- besf(y = y, x = x, coords=coords)
# res$b
# res$s
# res$e

}

