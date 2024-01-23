\name{addlearn_local}
\alias{addlearn_local}
\title{Additional learning of local processes and prediction for large samples}

\usage{
addlearn_local( mod, meig0 = NULL, x0 = NULL, xconst0=NULL, xgroup0=NULL,
           cl_num=NULL, cl=NULL, parallel=FALSE, ncores=NULL )
}
\description{
This function performs an additional learning of local variations in spatially varying coefficients. While the SVC model implemented in \code{\link{resf_vc}} or \code{\link{besf_vc}} can be less accurate for large samples (e.g., n > 5,000) due to a degeneracy/over-smoothing problem, this additional learning mitigates this problem by synthesizing/averaging the model with local SVC models. The resulting spatial prediction implemented in this function is expected to be more accurate than the resf_vc function.
}
\arguments{
\item{mod}{Outpot from \code{\link{resf_vc}} or  \code{\link{besf_vc}} function}
\item{meig0}{Moran eigenvectors at prediction sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at prediction sites whose coefficients are allowed to vary across geographical space (N_0 x K). Default is NULL}
\item{xconst0}{Matrix of explanatory variables at prediction sites whose coefficients are assumed constant (or NVC) across space (N_0 x K_const). Default is NULL}
\item{xgroup0}{Matrix of group indeces at prediction sites that may be group IDs (integers) or group names (N_0 x K_g). Default is NULL}
\item{cl_num}{Number of local sub-models being aggregated/averaged. If NULL, the number is determined so that the number of samples per sub-model equals approximately 600. Default is NULL}
\item{cl}{Vector of cluster ID for each sample (N x 1). If specified, the local sub-models are given by this ID. If NULL, k-means clustering based on spatial coordinates is performed to obtain spatial clusters each of which contain approximately 600 samples. Default is NULL}
\item{parallel}{If TRUE, the model is estimated through parallel computation. The default is FALSE for \code{\link{resf_vc}} while TRUE for \code{\link{besf_vc}}}
\item{ncores}{Number of cores used for the parallel computation. If ncores = NULL and parallel = TRUE, the number of available cores is detected. Default is NULL}
}

\value{
\item{b_vc}{Matrix of estimated spatially varying coefficients (SVCs) on x (N x K)}
\item{bse_vc}{Matrix of standard errors for the SVCs on x (N x k)}
\item{z_vc}{Matrix of z-values for the SVCs on x (N x K)}
\item{p_vc}{Matrix of p-values for the SVCs on x (N x K)}
\item{c}{Matrix with columns for the estimated coefficients on xconst, their standard errors, z-values, and p-values (K_c x 4)}
\item{b_g}{List of K_g matrices with columns for the estimated group effects, their standard deviations, and t-values}
\item{s}{List of 2 elements summarizing variance parameters characterizing SVCs of each local sub-model. The first element contains standard deviations of each SVCs while the second elementcontains their Moran's I values that are scaled to take a value between 0 (no spatial dependence) and 1 (strongest positive spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{s_global}{The same variance parameters for the globa sub-model}
\item{s_g}{Vector of standard deviations of the group effects}
\item{e}{Error statistics. It includes residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC)}
\item{pred}{Matrix of predicted values for y (pred) and their standard errors (pred_se) (N x 2)}
\item{resid}{Vector of residuals (N x 1)}
\item{cl}{Vector of cluster ID being used (N x 1)}
\item{pred0}{Matrix of predicted values for y (pred) and their standard errors (pred_se) at prediction sites (N_0 x 2)}
\item{b_vc0}{Matrix of estimated spatially varying coefficients (SVCs) at prediction sites (N_0 x K)}
\item{bse_vc0}{Matrix of standard errors for the SVCs at prediction sites (N_0 x k)}
\item{z_vc0}{Matrix of z-values for the SVCs at prediction sites (N x K)}
\item{p_vc0}{Matrix of p-values for the SVCs at prediction sites (N x K)}
\item{other}{List of other outputs, which are internally used}
}
\references{
Murakami, D., Sugasawa, S., T., Seya, H., and Griffith, D.A. (2024) Sub-model aggregation-based scalable eigenvector spatial filtering: application to spatially varying coefficient modeling. Arxiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{resf_vc}}, \code{\link{besf_vc}}
}
\examples{
require(spdep)
data(house)
dat0    <- data.frame(house@coords,house@data)
dat     <- dat0[dat0$yrbuilt>=1980,]

###### purpose 1: improve SVC modeling accuracy ######
###### (i.e., addressing the over-smoothing problem) #
y	      <- log(dat[,"price"])
x       <- dat[,c("age","rooms")]
xconst  <- dat[,c("lotsize","s1994","s1995","s1996","s1997","s1998")]
coords  <- dat[ ,c("long","lat")]
meig    <- meigen_f( coords )

## Not run: Remove # and run
# res0  <- resf_vc(y = y,x = x, xconst = xconst, meig = meig)
# res   <- addlearn_local(res0) # It adjusts SVCs to model local patterns
# res

####### parallel version for very large samples (e.g., n >100,000)
# bes0  <- besf_vc(y = y,x = x, xconst = xconst, coords=coords)
# bes	  <- addlearn_local( bes0 )


####### purpose 2: improve predictive accuracy ########

#samp    <- sample( dim( dat )[ 1 ], 2500)
#d       <- dat[ samp, ]    ## Data at observed sites
#y	     <- log(d[,"price"])
#x       <- d[,c("age","rooms")]
#xconst  <- d[,c("lotsize","s1994","s1995","s1996","s1997","s1998")]
#coords  <- d[ ,c("long","lat")]

#d0      <- dat[-samp, ]    ## Data at observed sites
#y0	     <- log(d0[,"price"])
#x0      <- d0[,c("age","rooms")]
#xconst0 <- d0[,c("lotsize","s1994","s1995","s1996","s1997","s1998")]
#coords0 <- d0[ ,c("long","lat")]

#meig    <- meigen_f( coords )
#res0    <- resf_vc(y = y,x = x, xconst = xconst, meig = meig)
#meig0   <- meigen0( meig=meig, coords0=coords0 )
#res     <- addlearn_local(res0, meig0=meig0, x0=x0, xconst0=xconst0) #
#pred    <- res$pred0       ## Predictive values

}
