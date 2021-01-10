\name{predict0_vc}
\alias{predict0_vc}
\title{Prediction of explained variables and spatially varying coefficients}
\usage{
predict0_vc( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL )
}
\description{
This function predicts explained variables and spatially and non-spatially varying coefficients. The Nystrom extension is used to perform a prediction minimizing the expected prediction error
}
\arguments{
\item{mod}{Output from \code{\link{resf_vc}} or  \code{\link{besf_vc}}}
\item{meig0}{Moran eigenvectors at predicted sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at predicted sites whose coefficients are allowed to vary across geographical space (N_0 x K). Default is NULL}
\item{xgroup0}{Matrix of group indeces that may be group IDs (integers) or group names (N_0 x K_{group}). Default is NULL}
\item{xconst0}{Matrix of explanatory variables at predicted sites whose coefficients are assumed constant (or NVC) across space (N_0 x K_const). Default is NULL}
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
x0      <- d0[,c("ZN", "LSTAT")]
xconst0 <- d0[,c("CRIM", "NOX", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "RM")]
coords0 <- d0[,c("LON", "LAT")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig )

############ Spatial prediction of y and spatially varying coefficients
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0_vc( mod = mod, x0 = x0, xconst0=xconst0, meig0 = meig0 )

pred0$pred[1:10,]  # Predicted explained variables
pred0$b_vc[1:10,]  # Predicted SVCs
pred0$bse_vc[1:10,]# Predicted standard errors of the SVCs
pred0$t_vc[1:10,]  # Predicted t-values of the SNVCs
pred0$p_vc[1:10,]  # Predicted p-values of the SNVCs

############ or spatial prediction of spatially varying coefficients only
# pred00  <- predict0_vc( mod = mod, meig0 = meig0 )
# pred00$b_vc[1:10,]
# pred00$bse_vc[1:10,]
# pred00$t_vc[1:10,]
# pred00$p_vc[1:10,]

######################## If SNVCs are assumed on x
# mod2	   <- resf_vc(y=y, x=x, xconst=xconst, meig=meig, x_nvc=TRUE,xconst_nvc=TRUE )
# pred02  <- predict0_vc( mod = mod2, x0 = x0, xconst0=xconst0 ,meig0 = meig0 )
# pred02$pred[1:10,]  # Predicted explained variables
# pred02$b_vc[1:10,]  # Predicted SNVCs
# pred02$bse_vc[1:10,]# Predicted standard errors of the SNVCs
# pred02$t_vc[1:10,]  # Predicted t-values of the SNVCs
# pred02$p_vc[1:10,]  # Predicted p-values of the SNVCs

}
