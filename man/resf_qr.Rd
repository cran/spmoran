\name{resf_qr}
\alias{resf_qr}
\title{Spatial filter unconditional quantile regression}
\usage{
resf_qr( y, x = NULL, meig, tau = NULL, boot = TRUE, iter = 200, cl=NULL )
}
\description{
This function estimates the spatial filter unconditional quantile regression (SF-UQR) model.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{tau}{The quantile(s) to be modeled. It must be a number (or a vector of numbers) strictly between 0 and 1. By default, tau = c(0.1, 0.2, ..., 0.9)}
\item{boot}{If it is TRUE, confidence intervals of regression coefficients are estimated by a semiparametric bootstrapping. Default is TRUE}
\item{iter}{The number of bootstrap replications. Default is 200}
\item{cl}{Number of cores used for the parallel computation. If cl=NULL, which is the default, the number of available cores is detected and used}
}
\value{
\item{b}{Matrix of estimated regression coefficients (K x Q), where Q is the number of quantiles (i.e., the length of tau)}
\item{r}{Matrix of estimated random coefficients on Moran eigenvectors (L x Q)}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements denote the standard error and the Moran's I value of the estimated spatially dependent component, respectively. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{e}{Vector whose elements are residual standard error (resid_SE) and adjusted quasi conditional R2 (quasi_adjR2(cond))}
\item{B}{Q matrices (K x 4) summarizing bootstrapped estimates for the regression coefficients. Columns of these matrices consist of the estimated coefficients, the lower and upper bounds for the 95 percent confidencial intervals, and p-values. It is returned if boot = TRUE}
\item{S}{Q matrices (2 x 3) summarizing bootstrapped estimates for the variance parameters. Columns of these matrices consist of the estimated parameters, the lower and upper bounds for the 95 percent confidencial intervals. It is returned if boot = TRUE}
\item{B0}{List of Q matrices (K x iter) summarizing bootstrapped coefficients. The q-th matrix consists of the coefficients on the q-th quantile. Effective if boot = TRUE}
\item{S0}{List of Q matrices (2 x iter) summarizing bootstrapped variance parameters. The q-th matrix consists of the parameters on the q-th quantile. Effective if boot = TRUE}
}
\seealso{
\code{\link{plot_qr}}
}
\references{
Murakami, D. and Seya, H. (2017) Spatially filtered unconditional quantile regression. ArXiv.
}
\author{
Daisuke Murakami
}
\examples{
require(spdep)
data(boston)
y	<- boston.c[, "CMEDV" ]
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
coords  <- boston.c[,c("LAT","LON")]
meig 	  <- meigen(coords=coords)
res	    <- resf_qr(y=y,x=x,meig=meig, boot=FALSE)
res
plot_qr(res,1)     # Intercept
plot_qr(res,2)     # Coefficient on CRIM
plot_qr(res,1,"s") # spcomp_SE
plot_qr(res,2,"s") # spcomp_Moran.I/max(Moran.I)

###Not run
#res	<- resf_qr(y=y,x=x,meig=meig, boot=TRUE)
#res
#plot_qr(res,1)     # Intercept + 95 percent confidence interval (CI)
#plot_qr(res,2)     # Coefficient on CRIM + 95 percent CI
#plot_qr(res,1,"s") # spcomp_SE + 95 percent CI
#plot_qr(res,2,"s") # spcomp_Moran.I/max(Moran.I) + 95 percent CI
}

