\name{predict0}
\alias{predict0}
\title{Spatial predictions}
\usage{
predict0( mod, meig0, x0 = NULL, xconst0 = NULL, xgroup0 = NULL, offset0 = NULL,
weight0 = NULL, compute_se=FALSE, compute_quantile = FALSE )
}
\description{
It is a function for spatial prediction using the model estimated from \code{\link{esf}}, \code{\link{resf}}, or \code{\link{resf_vc}} function.
}
\arguments{
\item{mod}{Output from \code{\link{esf}} \code{\link{resf}}}, or \code{\link{resf_vc}}
\item{meig0}{Moran eigenvectors at prediction sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at prediction sites (N_0 x K). Each column of x0 must correspond to those in x in the input model (mod). Default is NULL}
\item{xconst0}{Effective for \code{\link{resf_vc}}. Matrix of explanatory variables at prediction sites whose coefficients are assumed constant across space (N_0 x K_const). Each column of xconst0 must correspond to those in xconst in the input model. Default is NULL}
\item{xgroup0}{Matrix/vector of group IDs at prediction sites that may be integer or name by group (N_0 x K_g). Default is NULL}
\item{offset0}{Vector of offset variables at prediction sites (N_0 x 1). Effective if y is count (see \code{\link{nongauss_y}}). Default is NULL}
\item{weight0}{Vector of weights for prediction sites (N_0 x 1). Required if compute_se = TRUE or compute_quantile = TRUE, and weight in the input model is not NULL}
\item{compute_se}{If TRUE, predictive standard error is evaulated. It is currently supported only for continuous variables. If nongauss is specified in the input model (mod), standard error for the transformed y is evaluated. Default is FALSE}
\item{compute_quantile}{If TRUE, Matrix of the quantiles for the predicted values (N x 15) is evaulated. It is currently supported only for continuous variables. Default is FALSE}
}
\value{
\item{pred}{Matrix with the first column for the predicted values (pred). The second and the third columns are the predicted trend component (xb) and the residual spatial process (sf_residual). If xgroup0 is specified, the fourth column is the predicted group effects (group). If tr_num > 0 or tr_nonneg ==TRUE (i.e., y is transformed) in the input model \code{\link{resf}}, another column including the predicted values in the transformed/normalized scale (pred_trans) is added. In addition, if compute_quantile =TRUE, predictive standard errors (pred_se) is evaluated and added as another column}
\item{pred_quantile}{Effective if compute_quantile = TRUE. Matrix of the quantiles for the predicted values (N x 15). It is useful for evaluating uncertainty in the predictive values}
\item{b_vc}{Matrix of estimated spatially (and non-spatially) varying coefficients (S(N)VCs) on x0 (N_0 x K)}
\item{bse_vc}{Matrix of estimated standard errors for the S(N)VCs (N_0 x K)}
\item{t_vc}{Matrix of estimated t-values for the S(N)VCs (N_0 x K)}
\item{p_vc}{Matrix of estimated p-values for the S(N)VCs (N_0 x K)}
\item{c_vc}{Matrix of estimated non-spatially varying coefficients (NVCs) on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
\item{cse_vc}{Matrix of standard errors for the NVCs on x0 (N x K).Effective if nvc =TRUE in \code{\link{resf}}}
\item{ct_vc}{Matrix of t-values for the NVCs on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
\item{cp_vc}{Matrix of p-values for the NVCs on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
}
\seealso{
\code{\link{meigen0}}
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

meig 	  <- meigen( coords = coords )
meig0 	<- meigen0( meig = meig, coords0 = coords0 )

############ Spatial prediction ############
#### model with residual spatial dependence
mod	    <- resf(y=y, x=x, meig=meig)
pred0   <- predict0( mod = mod, x0 = x0, meig0 = meig0 )
pred0$pred[1:5,]  # Predicted values

#### model with spatially varying coefficients (SVCs)
mod	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig )
pred0   <- predict0( mod = mod, x0 = x0, xconst0=xconst0, meig0 = meig0 )
pred0$pred[1:5,]  # Predicted values
pred0$b_vc[1:5,]  # SVCs
pred0$bse_vc[1:5,]# standard errors of the SVCs
pred0$t_vc[1:5,]  # t-values of the SNVCs
pred0$p_vc[1:5,]  # p-values of the SNVCs

plot(y0,pred0$pred[,1]);abline(0,1)


}
