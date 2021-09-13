\name{resf_vc}
\alias{resf_vc}
\title{Gaussian and non-Gaussian spatial regression models with varying coefficients}

\usage{
resf_vc(y, x, xconst = NULL, xgroup = NULL, weight = NULL, offset = NULL,
        x_nvc = FALSE, xconst_nvc = FALSE, x_sel = TRUE, x_nvc_sel = TRUE,
        xconst_nvc_sel = TRUE, nvc_num = 5, meig, method = "reml",
        penalty = "bic", maxiter = 30, nongauss = NULL,
        tr_nonneg = NULL, tr_num = NULL )
}
\description{
This model estimates regression coefficients, spatially varying coefficients (SVCs), non-spatially varying coefficients (NVC; coefficients varying with respect to explanatory variable value), SNVC (= SVC + NVC), group effects, and residual spatial dependence. The random-effects eigenvector spatial filtering, which is an approximate Gaussian process approach, is used for modeling the spatial process in coefficients and residuals. While the resf_vc function estimates a SVC model by default, the type of coefficients (constant, SVC, NVC, or SNVC) can be selected through a BIC/AIC minimization. This function is available for modeling Gaussian and non-Gaussian data including continuous and count data (see \code{\link{nongauss_y}}).

Note that SNVCs can be mapped just like SVCs. SNVC model is more robust against spurious correlation (multicollinearity) and stable than SVC models (see Murakami and Griffith, 2020).
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables with spatially varying coefficients (SVC) (N x K)}
\item{xconst}{Matrix of explanatory variables with constant coefficients (N x K_{c}). Default is NULL}
\item{xgroup}{Matrix of group IDs. The IDs may be group numbers or group names (N x K_{g}). Default is NULL}
\item{weight}{Vector of weights for samples (N x 1). When non-NULL, the adjusted R-squared value is evaluated for weighted explained variables. Default is NULL}
\item{offset}{Vector of offset variables (N x 1). Available if y is count (y_type = "count" is specified in the \code{\link{nongauss_y}} function). Default is NULL}
\item{x_nvc}{If TRUE, SNVCs are assumed on x. Otherwise, SVCs are assumed. Default is FALSE}
\item{xconst_nvc}{If TRUE, NVCs are assumed on xconst. Otherwise, constant coefficients are assumed. Default is FALSE}
\item{x_sel}{If TRUE, type of coefficient (SVC or constant) on x is selected through a BIC (default) or AIC minimization. If FALSE, SVCs are assumed across x. Alternatively, x_sel can be given by column number(s) of x. For example, if x_sel = 2, the coefficient on the second explanatory variable in x is SVC and the other coefficients are constants. The Default is TRUE}
\item{x_nvc_sel}{If TRUE, type of coefficient (NVC or constant) on x is selected through the BIC (default) or AIC minimization. If FALSE, NVCs are assumed across x. Alternatively, x_nvc_sel can be given by column number(s) of x. For example, if x_nvc_sel = 2, the coefficient on the second explanatory variable in x is NVC and the other coefficients are constants. The Default is TRUE}
\item{xconst_nvc_sel}{If TRUE, type of coefficient (NVC or constant) on xconst is selected through the BIC (default) or AIC minimization. If FALSE, NVCs are assumed across xconst. Alternatively, xconst_nvc_sel can be given by column number(s) of xconst. For example, if xconst_nvc_sel = 2, the coefficient on the second explanatory variable in xconst is NVC and the other coefficients are constants.The Default is TRUE}
\item{nvc_num}{Number of basis functions used to model NVC. An intercept and nvc_num natural spline basis functions are used to model each NVC. Default is 5}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select varying coefficients and stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion  (2K). Default is "bic"}
\item{maxiter}{Maximum number of iterations. Default is 30}
\item{nongauss}{Parameter setup for modeling non-Gaussian continuous and count data. Output from \code{\link{nongauss_y}}}
\item{tr_nonneg}{Deprecated. Use the nongauss function}
\item{tr_num}{Deprecated. Use the nongauss function}
}
\details{
This function estimates Gaussian and non-Gaussian spatial model for continuous and count data. For non-Gaussian modeling, see \code{\link{nongauss_y}}.
}
\value{
\item{b_vc}{Matrix of estimated spatially and non-spatially varying coefficients (SNVC = SVC + NVC) on x (N x K)}
\item{bse_vc}{Matrix of standard errors for the SNVCs on x (N x k)}
\item{t_vc}{Matrix of t-values for the SNVCs on x (N x K)}
\item{p_vc}{Matrix of p-values for the SNVCs on x (N x K)}
\item{B_vc_s}{List summarizing estimated SVCs (in SNVC) on x. The four elements are the SVCs (N x K), the standard errors (N x K), t-values (N x K), and p-values (N x K), respectively}
\item{B_vc_n}{List summarizing estimated NVCs (in SNVC) on x. The four elements are the NVCs (N x K), the standard errors (N x K), t-values (N x K), and p-values (N x K), respectively}
\item{c}{Matrix with columns for the estimated coefficients on xconst, their standard errors, t-values, and p-values (K_{c} x 4). Effective if xconst_nvc = FALSE}
\item{c_vc}{Matrix of estimated NVCs on xconst (N x K_{c}). Effective if xconst_nvc = TRUE}
\item{cse_vc}{Matrix of standard errors for the NVCs on xconst (N x k_{c}). Effective if xconst_nvc = TRUE}
\item{ct_vc}{Matrix of t-values for the NVCs on xconst (N x K_{c}). Effective if xconst_nvc = TRUE}
\item{cp_vc}{Matrix of p-values for the NVCs on xconst (N x K_{c}). Effective if xconst_nvc = TRUE}
\item{b_g}{List of K_{g} matrices with columns for the estimated group effects, their standard errors, and t-values}
\item{s}{List of variance parameters in the SNVC (SVC + NVC) on x. The first element is a 2 x K matrix summarizing variance parameters for SVC. The (1, k)-th element is the standard error of the k-th SVC, while the (2, k)-th element is the Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (strongest spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked. The second element of s is the vector of standard errors of the NVCs}
\item{s_c}{Vector of standard errors of the NVCs on xconst}
\item{s_g}{Vector of standard errors of the group effects}
\item{vc}{List indicating whether SVC/NVC are removed or not during the BIC/AIC minimization. 1 indicates not removed (replaced with constant) whreas 0 indicates removed}
\item{e}{Error statistics. If y_type="continuous", it includes residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC) while rlogLik is replaced with log-likelihood (logLik) if method = "ml". If y_type="count", it includes deviance explained, Gaussian likelihood approximating the model, (Gaussian) AIC, and BIC}
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

Murakami, D., and Griffith, D.A. (2020) Balancing spatial and non-spatial variations in varying coefficient modeling: a remedy for spurious correlation. ArXiv.

Griffith, D. A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}, \code{\link{coef_marginal}}, \code{\link{besf_vc}}
}
\examples{
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
############## Gaussian SVC models ##################
#####################################################

#### SVC or constant coefficients on x ##############

res	    <- resf_vc(y=y,x=x,xconst=xconst,meig=meig )
res
plot_s(res,0) # Spatially varying intercept
plot_s(res,1) # 1st SVC (Not shown because the SVC is estimated constant)
plot_s(res,2) # 2nd SVC

#### SVC on x #######################################

#res2	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_sel = FALSE )

#### Group-level SVC or constant coefficients on x ##
#### Group-wise random intercepts ###################

#meig_g <- meigen(coords, s_id=xgroup)
#res3	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig_g,xgroup=xgroup)

#####################################################
############## Gaussian SNVC models #################
#####################################################

#### SNVC, SVC, NVC, or constant coefficients on x ###

#res4	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE)

#### SNVC, SVC, NVC, or constant coefficients on x ###
#### NVC or Constant coefficients on xconst ##########

#res5	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE, xconst_nvc=TRUE)
#plot_s(res5,0)            # Spatially varying intercept
#plot_s(res5,2)            # Spatial plot of the SNVC (SVC + NVC) on x[,2]
#plot_s(res5,2,btype="svc")# Spatial plot of SVC in the SNVC
#plot_s(res5,2,btype="nvc")# Spatial plot of NVC in the SNVC
#plot_n(res5,2)            # 1D plot of the NVC

#plot_s(res5,6,xtype="xconst")# Spatial plot of the NVC on xconst[,6]
#plot_n(res5,6,xtype="xconst")# 1D plot of the NVC on xconst[,6]


#####################################################
############## Non-Gaussian SVC models ##############
#####################################################

#### Generalized model for continuous data ##########
# - Probability distribution is estimated from data

#ng6    <- nongauss_y( tr_num = 2 )# 2 SAL transformations to Gaussianize y
#res6	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng6 )
#res6                   # tr_num may be selected by comparing BIC (or AIC)

#coef_marginal_vc(res6) # marginal effects from x (dy/dx)
#plot(res6$pdf,type="l") # Estimated probability density function
#res6$skew_kurt          # Skew and kurtosis of the estimated PDF
#res6$pred_quantile[1:2,]# predicted value by quantile


#### Generalized model for non-negative continuous data
# - Probability distribution is estimated from data

#ng7    <- nongauss_y( tr_num = 2, y_nonneg = TRUE )
#res7	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng7 )
#coef_marginal_vc(res7)

#### Overdispersed Poisson model for count data #####
# - y is assumed as a count data

#ng8    <- nongauss_y( y_type = "count" )
#res8	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng8 )

#### Generalized model for count data ###############
# - y is assumed as a count data
# - Probability distribution is estimated from data

#ng9    <- nongauss_y( y_type = "count", tr_num = 2 )
#res9	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, nongauss = ng9 )


}
