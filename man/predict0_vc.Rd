\name{predict0_vc}
\alias{predict0_vc}
\title{Prediction of explained variables and spatially varying coefficients}
\usage{
predict0_vc( mod, meig0, x0 = NULL, xgroup0 = NULL, xconst0 = NULL )
}
\description{
This function predicts explained variables and spatially varying coefficients using a Moran's eigenvector-based spatially varying coefficient model. The Nystrom extension is used to minimize the expected prediction error
}
\arguments{
\item{mod}{spatially varying coefficient model estimates. Output from \code{\link{resf_vc}}}
\item{meig0}{Moran's eigenvectors at predicted sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at predicted sites whose coefficients are allowed to vary across geographical space (N_0 x K). Default is NULL}
\item{xgroup0}{Matrix of group indeces that may be group IDs (integers) or group names (N_0 x K_{group}). Default is NULL}
\item{xconst0}{Matrix of explanatory variables at predicted sites whose coefficients are assumed constant across space (N_0 x K_const). Default is NULL}
}
\value{
\item{pred}{Matrix with the first column for the predicted values. The second and the third columns are the trend component (i.e., component explained by x0 and xconst0) and the spatial component in the predicted values (N x 3). The output has anoteher column including the estimated group effects if group effects are considered}
\item{b_vc}{Matrix of estimated spatially varying coefficients (SVCs) on x0 (N_0 x K)}
\item{bse_vc}{Matrix of estimated standard errors for the SVCs (N_0 x K)}
\item{t_vc}{Matrix of estimated t-values for the SVCs (N_0 x K)}
\item{p_vc}{Matrix of estimated p-values for the SVCs (N_0 x K)}

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
coords  <- d[,c("LAT","LON")]

d0      <- boston.c[-samp, ]    ## Data at unobserved sites
x0      <- d0[,c("ZN", "LSTAT")]
xconst0 <- d0[,c("CRIM", "NOX", "AGE", "DIS", "RAD", "TAX", "PTRATIO", "B", "RM")]
coords0 <- d0[,c("LAT","LON")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod	    <- resf_vc(y=y, x=x, xconst=xconst, meig=meig )

############ Spatial prediction of y and spatially varying coefficients
meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0_vc( mod = mod, x0 = x0, xconst0=xconst0, meig0 = meig0 )

pred0$pred[1:10,]
pred0$b_vc[1:10,]
pred0$bse_vc[1:10,]
pred0$t_vc[1:10,]
pred0$p_vc[1:10,]

############ or spatial prediction of spatially varying coefficients only
# pred00  <- predict0_vc( mod = mod, meig0 = meig0 )
# pred00$b_vc[1:10,]
# pred00$bse_vc[1:10,]
# pred00$t_vc[1:10,]
# pred00$p_vc[1:10,]
}
