\name{besf}
\alias{besf}
\title{Spatial regression with RE-ESF for very large samples}
\usage{
besf( y, x = NULL, nvc = FALSE, nvc_sel = TRUE, coords, s_id = NULL,
      covmodel="exp", enum = 200, method = "reml", penalty = "bic", nvc_num = 5,
      maxiter = 30, bsize = 4000, cl = NULL )
}
\description{
Memory-free implementation of RE-ESF-based spatial regression for very large samples. This model estimates residual spatial dependence, constant coefficients, and non-spatially varying coefficients (NVC; coefficients varying with respect to explanatory variable value).
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K)}
\item{nvc}{If TRUE, NVCs are assumed on x. Otherwise, constant coefficients are assumed. Default is FALSE}
\item{nvc_sel}{If TRUE, type of coefficients (NVC or constant) is selected through a BIC (default) or AIC minimization. If FALSE, NVCs are assumed across x. Alternatively, nvc_sel can be given by column number(s) of x. For example, if nvc_sel = 2, the coefficient on the second explanatory variable in x is NVC and the other coefficients are constants. The Default is TRUE}
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{s_id}{Optional. ID specifying groups modeling spatially dependent process (N x 1). If it is specified, group-level spatial process is estimated. It is useful. e.g., for multilevel modeling (s_id is given by the group ID) and panel data modeling (s_id is given by individual location id). Default is NULL}
\item{covmodel}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel}
\item{enum}{Number of Moran eigenvectors to be used for spatial process modeling (scalar). Default is 200}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select type of coefficients (NVC or constant) to stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion  (2K) (see Muller et al., 2013). Default is "bic"}
\item{nvc_num}{Number of basis functions used to model NVC. An intercept and nvc_num natural spline basis functions are used to model each NVC. Default is 5}
\item{maxiter}{Maximum number of iterations. Default is 30}
\item{bsize}{Block/badge size. bsize x bsize elements are iteratively processed during the parallelized computation. Default is 4000}
\item{cl}{Number of cores used for the parallel computation. If cl = NULL, the number of available cores is detected. Default is NULL}
}
\value{
\item{b}{Matrix with columns for the estimated coefficients on x, their standard errors, z-values, and p-values (K x 4). Effective if nvc =FALSE}
\item{c_vc}{Matrix of estimated NVCs on x (N x K). Effective if nvc =TRUE}
\item{cse_vc}{Matrix of standard errors for the NVCs on x (N x K). Effective if nvc =TRUE}
\item{ct_vc}{Matrix of t-values for the NVCs on x (N x K). Effective if nvc =TRUE}
\item{cp_vc}{Matrix of p-values for the NVCs on x (N x K). Effective if nvc =TRUE}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements denote the standard error and the Moran's I value of the estimated spatially dependent component, respectively. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{vc}{List indicating whether NVC are removed or not during the BIC/AIC minimization. 1 indicates not removed whreas 0 indicates removed}
\item{r}{Vector of estimated random coefficients on Moran's eigenvectors (L x 1)}
\item{sf}{Vector of estimated spatial dependent component (N x 1)}
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outputs, which are internally used}
}
\references{
Griffith, D. A. (2003). Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

Murakami, D. and Griffith, D.A. (2019) A memory-free spatial additive mixed modeling for big spatial data. Japan Journal of Statistics and Data Science. DOI:10.1007/s42081-019-00063-x.
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
x	<- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                       "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
xgroup  <- boston.c[,"TOWN"]
coords  <- boston.c[,c("LON", "LAT")]

######## Regression considering spatially dependent residuals
#res	  <- besf(y = y, x = x, coords=coords)
#res

######## Regression considering spatially dependent residuals and NVC
######## (coefficients or NVC is selected)
#res2  <- besf(y = y, x = x, coords=coords, nvc = TRUE)

######## Regression considering spatially dependent residuals and NVC
######## (all the coefficients are NVCs)
#res3  <- besf(y = y, x = x, coords=coords, nvc = TRUE, nvc_sel=FALSE)

}

