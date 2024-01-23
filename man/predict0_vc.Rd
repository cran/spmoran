\name{predict0_vc}
\alias{predict0_vc}
\title{Spatial predictions for explained variables and spatially varying coefficients}
\usage{
predict0_vc( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL,
             offset0 = NULL, weight0 = NULL, compute_se=FALSE, compute_quantile = FALSE )
}
\description{
This function predicts explained variables and spatially and non-spatially varying coefficients. IMPORTANT: This function will be retired in the next update. Use the predict0 function (auguments are the same).
}
\arguments{
\item{mod}{Output from \code{\link{resf_vc}} or  \code{\link{besf_vc}}}
\item{meig0}{Moran eigenvectors at prediction sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at prediction sites whose coefficients are allowed to vary across geographical space (N_0 x K). Default is NULL}
\item{xgroup0}{Matrix of group indeces at prediction sites that may be group IDs (integers) or group names (N_0 x K_g). Default is NULL}
\item{xconst0}{Matrix of explanatory variables at prediction sites whose coefficients are assumed constant (or NVC) across space (N_0 x K_const). Default is NULL}
\item{offset0}{Vector of offset variables at prediction sites (N x 1). Available if y is count (see \code{\link{nongauss_y}}). Default is NULL}
\item{weight0}{Vector of weights for predicted sites (N_0 x 1). Required if compute_se = TRUE or compute_quantile = TRUE}
\item{compute_se}{If TRUE, predictive standard error is evaulated. It is currently supported only for continuous variables. If nongauss is specified in mod, standard error for the transformed y is evaluated. Default is FALSE}
\item{compute_quantile}{If TRUE, Matrix of the quantiles for the predicted values (N x 15) is evaulated. Default is FALSE}
}
\value{
\item{pred}{Matrix with the first column for the predicted values (pred). The second and the third columns are the predicted trend component (i.e., component explained by x0 and xconst0) (xb) and the residual spatial process (sf_residual). If xgroup0 is specified, the fourth column is the predicted group effects (group) If tr_num > 0 or tr_nonneg ==TRUE (i.e., y is transformed) in \code{\link{resf_vc}}, another column including the predicted values in the transformed/normalized scale (pred_trans) is inserted into the second column}
\item{b_vc}{Matrix of estimated spatially (and non-spatially) varying coefficients (S(N)VCs) on x0 (N_0 x K)}
\item{bse_vc}{Matrix of estimated standard errors for the S(N)VCs (N_0 x K)}
\item{t_vc}{Matrix of estimated t-values for the S(N)VCs (N_0 x K)}
\item{p_vc}{Matrix of estimated p-values for the S(N)VCs (N_0 x K)}
\item{c_vc}{Matrix of estimated non-spatially varying coefficients (NVCs) on xconst0 (N_0 x K)}
\item{cse_vc}{Matrix of estimated standard errors for the NVCs (N_0 x K)}
\item{ct_vc}{Matrix of estimated t-values for the NVCs (N_0 x K)}
\item{cp_vc}{Matrix of estimated p-values for the NVCs (N_0 x K)}
}

\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.

Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.
}
\seealso{
\code{\link{meigen0}}, \code{\link{predict0}}
}

\examples{
require(spdep)
data(boston)
samp    <- sample( dim( boston.c )[ 1 ], 300)

d       <- boston.c[ samp, ]    ## Data at observed sites
y	      <- d[, "CMEDV"]
x       <- d[,c("ZN", "LSTAT")]
xconst  <- d[,c("CRIM", "NOX", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "RM")]
coords  <- d[,c("LON", "LAT")]

d0      <- boston.c[-samp, ]    ## Data at unobserved sites
y0      <- d0[, "CMEDV"]
x0      <- d0[,c("ZN", "LSTAT")]
xconst0 <- d0[,c("CRIM", "NOX", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "RM")]
coords0 <- d0[,c("LON", "LAT")]

############ Model estimation
meig    <- meigen( coords )
mod	    <- resf_vc(y=y, meig=meig, x=x, xconst=xconst )

############ Spatial prediction of y and spatially varying coefficients
meig0   <- meigen0( meig=meig, coords0=coords0)
pred0   <- predict0_vc( mod = mod, meig0=meig0, x0 = x0, xconst0=xconst0 )

pred0$pred[1:5,]  # Predicted explained variables
pred0$b_vc[1:5,]  # Predicted SVCs
pred0$bse_vc[1:5,]# Predicted standard errors of the SVCs
pred0$z_vc[1:5,]  # Predicted z-values of the SVCs
pred0$p_vc[1:5,]  # Predicted p-values of the SVCs

plot(y0,pred0$pred[,1]);abline(0,1)

}
