\name{resf}
\alias{resf}
\title{spatial and spatio-temporal regression models}
\usage{
resf( y, x = NULL, xgroup = NULL, weight = NULL, offset = NULL,
      nvc = FALSE, nvc_sel = TRUE, nvc_num = 5, meig,
      method = "reml", penalty = "bic", nongauss = NULL )
}
\description{
This model estimates regression coefficients, coefficients varying depending on x (non-spatially varying coefficients; NVC), and group effects, considering residual spatial/spatio-temporal dependence. The random-effects eigenvector spatial filtering, which is an approximate Gaussian process approach, is used for modeling the residual dependence. If nonugauss is specified, non-Gaussian explained variables are Gaussianized using a compositional warping function (see \code{\link{nongauss_y}}). This augument allows the resf function to be applied to non-Gaussian explained variables, including count data.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{xgroup}{Matrix of group IDs. The IDs may be group numbers or group names (N x K_g). Default is NULL}
\item{weight}{Vector of weights for samples (N x 1). If non-NULL, the adjusted R-squared value is evaluated for weighted explained variables. Default is NULL}
\item{offset}{Vector of offset variables (N x 1). Available if y is count (y_type = "count" is specified in the \code{\link{nongauss_y}} function). Default is NULL}
\item{nvc}{If TRUE, a non-linear function of x (NVC; a spline function) is used as a varying coefficient. If FALSE, constant coefficients are assumed. Default is FALSE}
\item{nvc_sel}{If TRUE, type of each coefficient (NVC or constant) is selected through a BIC minimization. If FALSE, NVCs are assumed across x. Alternatively, nvc_sel can be given by column number(s) of x. For example, if nvc_sel = 2, the coefficient on the second explanatory variable is NVC and the other coefficients are constants. Default is TRUE}
\item{nvc_num}{Number of natural spline basis functions to be used to model NVC. Default is 5}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select type of coefficients (NVC or constant) to stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion (2K). Default is "bic"}
\item{nongauss}{Parameter setup for modeling non-Gaussian continuous data or count data. Output from \code{\link{nongauss_y}}}
}
\details{
For modeling non-Gaussian data including count data, see \code{\link{nongauss_y}}.
}
\value{
\item{b}{Matrix with columns for the estimated constant coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{b_g}{List of K_g matrices with columns for the estimated group effects, their standard errors, and t-values}
\item{c_vc}{Matrix of estimated NVCs on x (N x K). Effective if nvc = TRUE}
\item{cse_vc}{Matrix of standard errors for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{ct_vc}{Matrix of t-values for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{cp_vc}{Matrix of p-values for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements are the standard deviation and the Moran's I value of the estimated spatially (and temporally) dependent process, respectively. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{s_c}{Vector of standard deviations of the NVCs on xconst}
\item{s_g}{Vector of estimated standard deviations of the group effects}
\item{e}{Error statistics. When y_type="continuous", it includes residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). rlogLik is replaced with log-likelihood (logLik) if method = "ml". resid_SE is replaced with the residual standard error for the transformed y (resid_SE_trans) if nongauss is specified. When y_type="count", the error statistics contains root mean squared error (RMSE), Gaussian likelihood approximating the model, AIC and BIC based on the likelihood, and the proportion of the null deviance explained by the model (deviance explained (\%)). deviance explained, which is also used in the mgcv package, corresponds to the adjusted R2 in case of the linear regression}
\item{vc}{List indicating whether NVC are removed or not during the BIC minimization. 1 indicates not removed whreas 0 indicates removed}
\item{r}{Vector of estimated random coefficients on the Moran's eigenvectors (L x 1)}
\item{sf}{Vector of estimated spatial dependent component (N x 1)}
\item{pred}{Matrix of predicted values for y (pred) and their standard errors (pred_se) (N x 2). If y is transformed by specifying \code{\link{nongauss_y}}, the predicted values in the transformed/normalized scale are added as another column named pred_trans}
\item{pred_quantile}{Matrix of the quantiles for the predicted values (N x 15). It is useful to evaluate uncertainty in the predictive value}
\item{tr_par}{List of the parameter estimates for the tr_num SAL transformations. The k-th element of the list includes the four parameters for the k-th SAL transformation (see \code{\link{nongauss_y}})}
\item{tr_bpar}{The estimated parameter in the Box-Cox transformation}
\item{tr_y}{Vector of the transformed explaied variables}
\item{resid}{Vector of residuals (N x 1)}
\item{pdf}{Matrix whose first column consists of evenly spaced values within the value range of y and the second column consists of the estimated value of the probability density function for y if y_type in \code{\link{nongauss_y}} is "continuous" and probability mass function (PMF) if y_type = "count". If offset is specified (and y_type = "count"), the PMF given median offset value is evaluated}
\item{skew_kurt}{Skewness and kurtosis of the estimated probability density/mass function of y}
\item{other}{List of other outputs, which are internally used}
}
\references{
Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

Murakami, D., and Griffith, D.A. (2020) Balancing spatial and non-spatial variations in varying coefficient modeling: a remedy for spurious correlation. Geographical Analysis, DOI: 10.1111/gean.12310.

Murakami, D., Kajita, M., Kajita, S. and Matsui, T. (2021) Compositionally-warped additive mixed modeling for a wide variety of non-Gaussian data. Spatial Statistics, 43, 100520.

Murakami, D., Shirota, S., Kajita, S., and Kajita, S. (2024) Fast spatio-temporally varying coefficient modeling with reluctant interaction selection. ArXiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}, \code{\link{coef_marginal}}, \code{\link{besf}}
}
\examples{

#####################################################
############ Spatial regression modeling ############
#####################################################
require(spdep);require(Matrix)
data(boston)
y	    <- boston.c[, "CMEDV" ]
x	    <- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                     "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
xgroup<- boston.c[,"TOWN"]
coords<- boston.c[,c("LON","LAT")]
meig 	<- meigen(coords=coords)
# meig<- meigen_f(coords=coords)  ## for large samples

#####################################################
####### Gaussian regression #########################
res	  <- resf(y = y, x = x, meig = meig)
res
plot_s(res)    ## spatially dependent component (intercept)

#### Group-wise random intercepts
#res2 <- resf(y = y, x = x, meig = meig, xgroup = xgroup)

#### Group-level spatial dependence (s_id) + random intercepts (xgroup)

#meig_g<- meigen(coords=coords, s_id = xgroup)
#res3  <- resf(y = y, x = x, meig = meig_g, xgroup = xgroup)

#### Coefficients varying depending on x

#res4  <- resf(y = y, x = x, meig = meig, nvc = TRUE)
#res4

#plot_s(res4)   # spatially dependent component (intercept)
#plot_s(res4,5) # spatial plot of the 5-th NVC
#plot_s(res4,6) # spatial plot of the 6-th NVC
#plot_s(res4,13)# spatial plot of the 13-th NVC

#plot_n(res4,5) # 1D plot of the 5-th NVC
#plot_n(res4,6) # 1D plot of the 6-th NVC
#plot_n(res4,13)# 1D plot of the 13-th NVC

#####################################################
###### Non-Gaussian regression ######################

#### Model for non-Gaussian continuous data
# - Probability distribution is estimated from data

#ng5    <- nongauss_y( tr_num = 2 )# 2 SAL transformations to Gaussianize y
#res5	  <- resf(y = y, x = x, meig = meig, nongauss = ng5)
#res5              ## tr_num may be selected by comparing BIC

#plot(res5$pdf,type="l") # Estimated probability density function
#res5$skew_kurt          # Skew and kurtosis of the estimated PDF
#res5$pred_quantile[1:2,]# predicted value by quantile
#coef_marginal(res5)     # Estimated marginal effects (dy/dx)


#### Model for non-Gaussian and non-negative continuous data
# - Probability distribution is estimated from data

#ng6    <- nongauss_y( tr_num = 2, y_nonneg = TRUE )
#res6	  <- resf(y = y, x = x, meig = meig, nongauss = ng6 )
#coef_marginal(res6)

#### Overdispersed Poisson model for count data
# - y: count data

#ng7    <- nongauss_y( y_type = "count" )
#res7	  <- resf(y = y, x = x, meig = meig, nongauss = ng7 )

#### General model for count data
# - y: count data
# - Probability distribution is estimated from data

#ng8    <- nongauss_y( y_type = "count", tr_num = 2 )
#res8	  <- resf(y = y, x = x, meig = meig, nongauss = ng8 )


#####################################################
################ STVC modeling ######################
#####################################################
# See \url{https://github.com/dmuraka/spmoran}

#require(spData)
#data(house)
#dat0    <- st_as_sf(house)
#dat     <- data.frame(st_coordinates(dat0), dat0)
#y	     <- log(dat[,"price"])
#x       <- dat[,c("lotsize","TLA", "rooms","beds")]

#byear   <- house$yrbuilt
#syear   <- as.numeric(as.character(house$syear))#factor -> numeric
#coords_z<- cbind(byear,syear)
#meig 	 <- meigen_f(coords=coords, coords_z=cbind(byear,syear),interact=TRUE)
#res9	   <- resf(y=y,x=x,meig=meig )
#res9
}

