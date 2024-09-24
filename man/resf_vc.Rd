\name{resf_vc}
\alias{resf_vc}
\title{spatial and spatio-temporal regression models with varying coefficients}

\usage{
resf_vc(y, x, xconst = NULL, xgroup = NULL, weight = NULL, offset = NULL,
        x_nvc = FALSE, xconst_nvc = FALSE, x_sel = TRUE, x_nvc_sel = TRUE,
        xconst_nvc_sel = TRUE,
        nvc_num = 5, meig, method = "reml", penalty = "bic", nongauss = NULL,
        miniter=NULL, maxiter = 30, tol = 1e-30 )
}
\description{
This function estimates spatially varying coefficients (SVC) or spatio-temporally varying coefficients (STVC), group effects, considering residual spatial/spatio-temporal dependence. A non-linear function of x (NVC) can be added on each SVC/STVC mainly to stablize the estimation (see Murakami and Griffith, 2020). Approximate Gaussian processes based on Moran eigenvectors are used for modeling the spatio-temporal processes. Type of coefficients (constant or varying) is selected through a BIC minimization. If nonugauss is specified, non-Gaussian explained variables are Gaussianized using a compositional warping function (see \code{\link{nongauss_y}}). This augument allows the resf function to be applied to non-Gaussian explained variables, including count data.

Note that, for very large samples, this function can overlook small-scale spatial variations. \code{\link{addlearn_local}} applies an model aggregation/averaging technique to address this problem.
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables assuming SVC/STVC (N x K)}
\item{xconst}{Matrix of explanatory variables assuming constant coefficients (N x K_c). Default is NULL}
\item{xgroup}{Matrix of group IDs for modeling group-wise random effects. The IDs may be group numbers or group names (N x K_g). Default is NULL}
\item{weight}{Vector of weights for samples (N x 1). If non-NULL, the adjusted R-squared value is evaluated for weighted explained variables. Default is NULL}
\item{offset}{Vector of offset variables (N x 1). Available if y is count (y_type = "count" is specified in the \code{\link{nongauss_y}} function). Default is NULL}
\item{x_nvc}{If TRUE, a non-linear function of x (NVC) is added on each varying coefficient on x to stablize the estimate. Default is FALSE}
\item{xconst_nvc}{If TRUE, NVCs is added on each constant coefficient on xconst model estimate non-linear influence from xconst}
\item{x_sel}{If TRUE, type of coefficient on x (STVC, SVC, or constant) is selected through a BIC minimization. If FALSE, S(T)VCs are assumed across x. Alternatively, x_sel can be given by column number(s) of x. For example, if x_sel = 2, the coefficient on the second explanatory variable is S(T)VC and the other coefficients are constants. The Default is TRUE}
\item{x_nvc_sel}{If TRUE, with/without NVC on x is selected. If FALSE, NVCs are assumed across x. Alternatively, x_nvc_sel can be given by column number(s) of x. For example, if x_nvc_sel = 2, the coefficient on the second explanatory variable is NVC and the other coefficients are constants. The Default is TRUE}
\item{xconst_nvc_sel}{If TRUE, with/without NVC on xconst is selected. If FALSE, NVCs are assumed across xconst. Alternatively, xconst_nvc_sel can be given by column number(s) of xconst. For example, if xconst_nvc_sel = 2, the coefficient on xconst[,2] becomes constant + NVC while the other coefficients become constants.The Default is TRUE}
\item{nvc_num}{Number of natural spline basis functions to be used in NVC. Default is 5}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty for model estimation and selection. "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion  (2K). Default is "bic"}
\item{nongauss}{Parameter setup for modeling non-Gaussian continuous and count data. Output from \code{\link{nongauss_y}}}
\item{miniter}{Minimum number of iterations. Default is NULL}
\item{maxiter}{Maximum number of iterations. Default is 30}
\item{tol}{The tolerance for matrix inversion. Some errors regarding singular fit can be avoided by reducing the value, but the output can be unstable. Default is 1e-30}
}
\details{
For modeling non-Gaussian data including count data, see \code{\link{nongauss_y}}.
}
\value{
\item{b_vc}{Matrix of estimated spatially/spatio-temporally varying coefficients (S(T)VC + NVC) on x (N x K)}
\item{bse_vc}{Matrix of standard errors for the varying coefficients on x (N x k)}
\item{t_vc}{Matrix of t-values for the coefficients on x (N x K)}
\item{p_vc}{Matrix of p-values for the coefficients on x (N x K)}
\item{B_vc_s}{List of the estimated S(T)VCs in b_vc (= S(T)VC + NVC). The elements are the S(T)VCs (N x K), the standard errors (N x K), t-values (N x K), and p-values (N x K), respectively}
\item{B_vc_n}{List of the estimated NVCs in b_vc (= S(T)VC + NVC). The elements are the NVCs (N x K), the standard errors (N x K), t-values (N x K), and p-values (N x K), respectively}
\item{c}{Matrix with columns for the estimated coefficients on xconst, their standard errors, t-values, and p-values (K_c x 4). Effective if xconst_nvc = FALSE}
\item{c_vc}{Matrix of estimated NVCs on xconst (N x K_c). Effective if xconst_nvc = TRUE}
\item{cse_vc}{Matrix of standard errors for the NVCs on xconst (N x k_c). Effective if xconst_nvc = TRUE}
\item{ct_vc}{Matrix of t-values for the NVCs on xconst (N x K_c). Effective if xconst_nvc = TRUE}
\item{cp_vc}{Matrix of p-values for the NVCs on xconst (N x K_c). Effective if xconst_nvc = TRUE}
\item{b_g}{List of K_g matrices with columns for the estimated group effects, their standard errors, and t-values}
\item{s}{List of the variance parameters for the varying coefficient on x. The first element is a 2 x K matrix summarizing variance parameters for S(T)VC. The (1, k)-th element is the standard deviation of the k-th SVC, while the (2, k)-th element is the Moran's I value that  is scaled to take a value between 0 (no spatial dependence) and 1 (strongest spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked. The second element of s is the vector of standard deviations of the NVCs}
\item{s_c}{Vector of standard deviations of the NVCs on xconst}
\item{s_g}{Vector of standard deviations of the group effects}
\item{vc}{List indicating whether S(T)VC/NVC are removed or not during the BIC minimization. 1 indicates not removed (replaced with constant) whreas 0 indicates removed}
\item{e}{Error statistics. When y_type="continuous", it includes residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). rlogLik is replaced with log-likelihood (logLik) if method = "ml". resid_SE is replaced with the residual standard error for the transformed y (resid_SE_trans) if nongauss is specified. When y_type="count", the error statistics includes root mean squared error (RMSE), Gaussian likelihood approximating the model, AIC and BIC based on the likelihood, and the proportion of the null deviance explained by the model (deviance explained (\%)). deviance explained, which is also used in the mgcv package, corresponds to the adjusted R2 in case of the linear regression}
\item{pred}{Matrix of predicted values for y (pred) and their standard errors (pred_se) (N x 2). If y is transformed by specifying \code{\link{nongauss_y}}, the predicted values in the transformed/normalized scale are added as another column named pred_trans}
\item{pred_quantile}{Matrix of the quantiles for the predicted values (N x 15). It is useful to evaluate uncertainty in the predictive value}
\item{tr_par}{List of the parameter estimates for the tr_num SAL transformations. The k-th element of the list includes the four parameters for the k-th SAL transformation (see \code{\link{nongauss_y}})}
\item{tr_bpar}{The estimated parameter in the Box-Cox transformation}
\item{tr_y}{Vector of the transformed explaied variables}
\item{resid}{Vector of residuals (N x 1)}
\item{pdf}{Matrix whose first column consists of evenly spaced values within the value range of y and the second column consists of the estimated value of the probability density function for y if y_type in \code{\link{nongauss_y}} is "continuous" and probability mass function if y_type = "count".  If offset is specified (and y_type = "count"), the PMF given median offset value is evaluated}
\item{skew_kurt}{Skewness and kurtosis of the estimated probability density/mass function of y}
\item{other}{List of other outputs, which are internally used}
}
\references{
Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.

Murakami, D., Kajita, M., Kajita, S. and Matsui, T. (2021) Compositionally-warped additive mixed modeling for a wide variety of non-Gaussian data. Spatial Statistics, 43, 100520.

Murakami, D., and Griffith, D.A. (2021) Balancing spatial and non-spatial variations in varying coefficient modeling: a remedy for spurious correlation. Geographical Analysis, DOI: 10.1111/gean.12310.

Murakami, D., Shirota, S., Kajita, S., and Kajita, S. (2024) Fast spatio-temporally varying coefficient modeling with reluctant interaction selection. ArXiv.

Griffith, D. A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}, \code{\link{coef_marginal}}, \code{\link{besf_vc}}, \code{\link{addlearn_local}}
}
\examples{

#####################################################
################ SVC modeling #######################
#####################################################
require(spdep)
data(boston)
y	      <- boston.c[, "CMEDV"]
x       <- boston.c[,c("CRIM", "AGE")]
xconst  <- boston.c[,c("ZN","DIS","RAD","NOX",  "TAX","RM", "PTRATIO", "B")]
xgroup  <- boston.c[,"TOWN"]
coords  <- boston.c[,c("LON", "LAT")]
meig 	  <- meigen(coords=coords)
# meig	<- meigen_f(coords=coords)  ## for large samples

#####################################################
####### Gaussian regression with SVC ################

res	    <- resf_vc(y=y,x=x,xconst=xconst,meig=meig )
res
plot_s(res,0) # Spatially varying intercept
plot_s(res,1) # 1st SVC (Not shown because the SVC is estimated constant)
plot_s(res,2) # 2nd SVC

#### For large samples (e.g., n > 5,000), the following
#### additional learning often improves the modeling accuracy
# res_adj<- addlearn_local(res)
# res_adj
# plot_s(res_adj,0)
# plot_s(res_adj,1)
# plot_s(res_adj,2)

#### Group-level SVC (s_id) + random intercepts (xgroup)
# meig_g <- meigen(coords, s_id=xgroup)
# res2	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig_g,xgroup=xgroup)


#####################################################
####### Gaussian regression with SVC + NVC ##########

# res3	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE)
# plot_s(res3,0)            # Spatially varying intercept
# plot_s(res3,1)            # Spatial plot of the varying coefficient (SVC + NVC) on x[,1]
# plot_s(res3,1,btype="svc")# Spatial plot of SVC in the coefficient
# plot_s(res3,1,btype="nvc")# Spatial plot of NVC in the coefficient
# plot_n(res3,1)            # 1D plot of the NVC


#####################################################
######## Non-Gaussian regression with SVC ###########

#### Model for non-Gaussian continuous data
# - Probability distribution is estimated from data

# ng4     <- nongauss_y( tr_num = 2 )# 2 SAL transformations to Gaussianize y
# res4	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng4 )
# res4                    # tr_num may be selected by comparing BIC

# coef_marginal_vc(res4)  # marginal effects from x (dy/dx)
# plot(res4$pdf,type="l") # Estimated probability density function
# res4$skew_kurt          # Skew and kurtosis of the estimated PDF
# res4$pred_quantile[1:2,]# predicted value by quantile


#### Model for non-Gaussian and non-negative continuous data
# - Probability distribution is estimated from data

## 2 SAL trans. + 1 Box-Cox trans. to Gaussianize y
# ng5     <- nongauss_y( tr_num = 2, y_nonneg = TRUE )
# res5	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng5 )
# coef_marginal_vc(res5)

#### Overdispersed Poisson model for count data
# - y: count data

#ng6    <- nongauss_y( y_type = "count" )
#res6	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng6 )

#### General model for count data
# - y: count data
# - Probability distribution is estimated from data

#ng7    <- nongauss_y( y_type = "count", tr_num = 2 )
#res7	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng7 )


#####################################################
################ STVC modeling ######################
#####################################################
# See \url{https://github.com/dmuraka/spmoran}

#require(spData)
#data(house)
#dat0    <- st_as_sf(house)
#dat     <- data.frame(st_coordinates(dat0), dat0)
#y	     <- log(dat[,"price"])
#x       <- dat[,c("lotsize","TLA")]
#xconst  <- dat[,c("rooms","beds")]

#byear   <- house$yrbuilt
#syear   <- as.numeric(as.character(house$syear))#factor -> numeric
#coords_z<- cbind(byear,syear)
#meig 	 <- meigen_f(coords=coords, coords_z=cbind(byear,syear),interact=TRUE)
#res8	   <- resf_vc(y=y,x=x,xconst=xconst,meig=meig )
#res8

## Varying intercept for byear <=1950 and syear==1998
#plot_s(res8,0, coords_z1_lim=c(-Inf, 1950),coords_z2_lim=1998)

## 1st STVCs which are significant at the 5 percent level, for byear <= 1950
#plot_s(res8,1, coords_z1_lim=c(-Inf, 1950), pmax=0.05)

## 2nd STVC for byear >= 1951
#plot_s(res8,2, coords_z1_lim=c(1951,Inf))
}
