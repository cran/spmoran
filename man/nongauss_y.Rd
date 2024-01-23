\name{nongauss_y}
\alias{nongauss_y}
\title{Parameter setup for modeling non-Gaussian continuous data and count data}
\usage{
nongauss_y( y_type = "continuous", y_nonneg = FALSE, tr_num = 0 )
}
\description{
Parameter setup for modeling non-Gaussian continuous data and count data. The SAL transformation (see details) is used to model a wide variety of non-Gaussian data without explicitly assuming data distribution (see Murakami et al., 2021 for further detail). In addition, Box-Cox transformation is used for non-negative continuous variables while another transformation approximating overdispersed Poisson distribution is used for count variables. The output from this function is used as an input of the resf and resf_vc functions. For further details about its implementation and case study examples, see Murakami (2021).

}
\arguments{
\item{y_type}{Type of explained variables y. "continuous" for continuous variables and "count" for count variables}
\item{y_nonneg}{Effective if y_type = "continuous". TRUE if y cannot take negative value. If y_nonneg = TRUE and tr_num = 0, the Box-Cox transformation is applied to y. If y_nonneg = TRUE and tr_num > 0, the Box-Cox transformation is applied first to roughly Gaussianize y. Then, the SAL transformation is iterated tr_num times to improve the modeling accuracy. Default is FALSE}
\item{tr_num}{Number of the SAL transformations (SinhArcsinh and Affine, where the use of "L" stems from the "Linear") applied to Gaussianize y. Default is 0}
}
\details{
If tr_num >0, the SAL transformation is iterated tr_num times to Gaussianize y. The SAL transformation is defined as SAL(y)=a+b*sinh(c*arcsinh(y)-d) where a,b,c,d are parameters. Based on Rios and Tobar (2019), the iteration of the SAL transformation approximates a wide variety of non-Gaussian distributions without explicitly assuming data distribution. The resf and resf_vc functions return tr_par, which is a list whose k-th element includes the a,b,c,d parameters used for the k-th SAL transformation.

In addition, for non-negative y (y_nonneg = TRUE), the Box-Cox transformation is applied prior to the iterative SAL transformation. tr_num and y_nonneg can be selected by comparing the BIC (or AIC) values across models. This compositionally-warped spatial regression approach is detailed in Murakami et al. (2021).

For count data (y_type = "count"), an overdispersed Poisson distribution (Gaussian approximation) is assumed. If tr_num > 0, the distribution is adjusted to fit the data (y) through the iterative SAL transformations. y_nonneg is ignored if y_type = "count".
}
\value{
\item{nongauss}{List of parameters for modeling non-Gaussian data}
}

\references{
Rios, G. and Tobar, F. (2019) Compositionally-warped Gaussian processes. Neural Networks, 118, 235-246.

Murakami, D. (2021) Transformation-based generalized spatial regression using the spmoran package: Case study examples, ArXiv.

Murakami, D., Kajita, M., Kajita, S. and Matsui, T. (2021) Compositionally-warped additive mixed modeling for a wide variety of non-Gaussian data. Spatial Statistics, 43, 100520.

Murakami, D., & Matsui, T. (2021). Improved log-Gaussian approximation for over-dispersed Poisson regression: application to spatial analysis of COVID-19. ArXiv, 2104.13588.
}

\seealso{
\code{\link{resf}}, \code{\link{resf_vc}}
}
\examples{

###### Regression for non-negative data (BC trans.)
ng1  <-nongauss_y( y_nonneg = TRUE )
ng1

###### General non-Gaussian regression for continuous data (two SAL trans.)
ng2  <-nongauss_y( tr_num = 2 )
ng2

###### General non-Gaussian regression for non-negative continuous data
ng3  <-nongauss_y( y_nonneg = TRUE, tr_num = 5 )
ng3

###### Over-dispersed Poisson regression for count data
ng4  <-nongauss_y( y_type = "count" )
ng4

###### A general non-Gaussian regression for count data
ng5  <-nongauss_y( y_type = "count", tr_num = 5 )
ng5

############################## Fitting example
require(spdep);require(Matrix)
data(boston)
y	    <- boston.c[, "CMEDV" ]
x	    <- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                     "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
xgroup<- boston.c[,"TOWN"]
coords<- boston.c[,c("LON","LAT")]
meig 	<- meigen(coords=coords)
res	  <- resf(y = y, x = x, meig = meig,nongauss=ng2)
res                    # Estimation results

plot(res$pdf,type="l") # Estimated probability density function
res$skew_kurt          # Skew and kurtosis of the estimated PDF
res$pred_quantile[1:2,]# predicted value by quantile
coef_marginal(res)     # Estimated marginal effects (dy/dx)
}

