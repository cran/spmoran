\name{predict0}
\alias{predict0}
\title{Spatial predictions}
\usage{
predict0( mod, meig0, x0 = NULL, xgroup0 = NULL, offset0 = NULL,
weight0 = NULL, compute_quantile = FALSE )
}
\description{
This function predicts explained variables using eigenvector spatial filtering (ESF) or random effects ESF. The Nystrom extension is used to perform a prediction minimizing the expected prediction error
}
\arguments{
\item{mod}{Output from \code{\link{esf}} or \code{\link{resf}}}
\item{meig0}{Moran eigenvectors at predicted sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at predicted sites (N_0 x K). Default is NULL}
\item{xgroup0}{Matrix of group IDs that may be group IDs (integers) or group names (N_0 x K_{group}). Default is NULL}
\item{offset0}{Vector of offset variables at predicted sites (N_0 x 1). Effective if y is count (see \code{\link{nongauss_y}}). Default is NULL}
\item{weight0}{Vector of weights for predicted sites (N_0 x 1). Required if compute_quantile = TRUE}
\item{compute_quantile}{If TRUE, Matrix of the quantiles for the predicted values (N x 15) is evaulated. It is currently supported only for continuous variables. Default is FALSE}
}
\value{
\item{pred}{Matrix with the first column for the predicted values (pred). The second and the third columns are the predicted trend component (xb) and the residual spatial process (sf_residual). If xgroup0 is specified, the fourth column is the predicted group effects (group). If tr_num > 0 or tr_nonneg ==TRUE (i.e., y is transformed) in \code{\link{resf}}, another column including the predicted values in the transformed/normalized scale (pred_trans) is inserted as the second column. In addition, if compute_quantile =TRUE, predictive standard errors (pred_se) is evaluated and inserted as another column}
\item{pred_quantile}{Effective if compute_quantile = TRUE. Matrix of the quantiles for the predicted values (N x 15). It is useful to evaluate uncertainty in the predictive value}
\item{c_vc}{Matrix of estimated non-spatially varying coefficients (NVCs) on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
\item{cse_vc}{Matrix of standard errors for the NVCs on x0 (N x K).Effective if nvc =TRUE in \code{\link{resf}}}
\item{ct_vc}{Matrix of t-values for the NVCs on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
\item{cp_vc}{Matrix of p-values for the NVCs on x0 (N x K). Effective if nvc =TRUE in \code{\link{resf}}}
}
\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.
}
\seealso{
\code{\link{meigen0}}, \code{\link{predict0_vc}}
}
\examples{
require(spdep)
data(boston)
samp    <- sample( dim( boston.c )[ 1 ], 400)

d       <- boston.c[ samp, ]  ## Data at observed sites
y	      <- d[, "CMEDV"]
x       <- d[,c("ZN","INDUS", "NOX","RM", "AGE", "DIS")]
coords  <- d[,c("LON", "LAT")]

d0      <- boston.c[-samp, ]  ## Data at unobserved sites
y0	    <- d0[, "CMEDV"]
x0      <- d0[,c("ZN","INDUS", "NOX","RM", "AGE", "DIS")]
coords0 <- d0[,c("LON", "LAT")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod	    <- resf(y=y, x=x, meig=meig)
## or
# mod   <- esf(y=y,x=x,meig=meig)

############ Spatial prediction
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0( mod = mod, x0 = x0, meig0 = meig0 )
pred0$pred[1:10,]

######################## If NVCs are assumed
#mod2	  <- resf(y=y, x=x, meig=meig, nvc=TRUE)
#pred02  <- predict0( mod = mod2, x0 = x0, meig0 = meig0 )
#pred02$pred[1:10,]  # Predicted explained variables
#pred02$c_vc[1:10,]  # Predicted NVCs

}
