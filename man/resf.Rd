\name{resf}
\alias{resf}
\title{Spatial regression for Gaussian and non-Gaussian continuous data}
\usage{
resf( y, x = NULL, xgroup = NULL, weight = NULL, nvc = FALSE, nvc_sel = TRUE,
      nvc_num = 5, meig, method = "reml", penalty = "bic",
      tr_num = 0, tr_nonneg = FALSE )
}
\description{
This model estimates residual spatial dependence, constant coefficients, non-spatially varying coefficients (NVC; coefficients varying depending on x), and group effects. The random-effects eigenvector spatial filtering (RE-ESF), which is a low rank Gaussian process approach, is used for the spatial dependence modeling. Compositionally-warping function and/or the Box-Cox transformation function is available for non-Gaussian continuous data (see Murakami et al., 2020 for further detail).
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables (N x K). Default is NULL}
\item{xgroup}{Matrix of group IDs. The IDs may be group numbers or group names (N x K_{group}). Default is NULL}
\item{weight}{Vector of weights for samples (N x 1). When non-NULL, the adjusted R-squared value is evaluated for weighted explained variables. Default is NULL}
\item{nvc}{If TRUE, non-spatiallly varying coefficients (NVCs; coefficients varying with respect to explanatory variable value) are asumed. If FALSE, constant coefficients are assumed. Default is FALSE}
\item{nvc_sel}{If TRUE, type of each coefficient (NVC or constant) is selected through a BIC (default) or AIC minimization. If FALSE, NVCs are assumed across x. Alternatively, nvc_sel can be given by column number(s) of x. For example, if nvc_sel = 2, the coefficient on the second explanatory variable is NVC and the other coefficients are constants. Default is TRUE}
\item{nvc_num}{Number of basis functions used to model NVC. An intercept and nvc_num natural spline basis functions are used to model each NVC. Default is 5}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select type of coefficients (NVC or constant) to stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion (2K). Default is "bic"}
\item{tr_num}{Number of the SAL transformations ( SinhArcsinh and Affine, where the use of "L" stems from the "Linear") used to transform non-Gaussian explained variables to Gaussian variables. Default is 0}
\item{tr_nonneg}{If TRUE, the Box-Cox transfromation used to transform positive non-Gaussian explained variables to Gaussian variables. If tr_num > 0 and tr_nonneg == TRUE, the Box-Cox transformation is applied first. Then, th SAL transformation is applied tr_num times. Default is FALSE}
}
\details{
If tr_num >0, the resf function iterates the SAL transformation tr_num times to transform the explained variables to Gaussian variables. The SAL transformation is defined as SAL(y)=a+b*sinh(c*arcsinh(y)-d) where {a,b,c,d} are parameters. Based on Rois and Tober (2019), an iteration of the SAL transformation approximates a wide variety of non-Gaussian distributions without explicitly assuming data distribution. As a result, our spatial regression approach is applicable to a wide variety of non-Gaussian continuous data too. For non-negative explained variables, a Box-Cox transformation is available prior to the SAL transformations by specifying tr_nonneg >0. tr_num and tr_nonneg can be selected by comparing the BIC (or AIC) value across models. This compositionally-warped spatial regression approach is detailed in Murakami et al. (2020).
}
\value{
\item{b}{Matrix with columns for the estimated constant coefficients on x, their standard errors, t-values, and p-values (K x 4)}
\item{b_g}{List of K_{group} matrices with columns for the estimated group effects, their standard errors, and t-values}
\item{c_vc}{Matrix of estimated NVCs on x (N x K). Effective if nvc = TRUE}
\item{cse_vc}{Matrix of standard errors for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{ct_vc}{Matrix of t-values for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{cp_vc}{Matrix of p-values for the NVCs on x (N x K). Effective if nvc = TRUE}
\item{s}{Vector of estimated variance parameters (2 x 1). The first and the second elements are the standard error and the Moran's I value of the estimated spatially dependent process, respectively. The Moran's I value is scaled to take a value between 0 (no spatial dependence) and 1 (the maximum possible spatial dependence). Based on Griffith (2003), the scaled Moran'I value is interpretable as follows: 0.25-0.50:weak; 0.50-0.70:moderate; 0.70-0.90:strong; 0.90-1.00:marked}
\item{s_c}{Vector of standard errors of the NVCs on xconst}
\item{s_g}{Vector of estimated standard errors of the group effects}
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
\item{vc}{List indicating whether NVC are removed or not during the BIC/AIC minimization. 1 indicates not removed whreas 0 indicates removed}
\item{r}{Vector of estimated random coefficients on Moran's eigenvectors (L x 1)}
\item{sf}{Vector of estimated spatial dependent component (N x 1)}
\item{pred}{Vector of predicted values (N x 1). If tr_num > 0 or tr_nonneg == TRUE (i.e., y is trnsformed), another column including the predicted values in the transformed/normalized scale (pred_trans) is inserted into the second column}
\item{tr_par}{List of the estimated parameters in the tr_num SAL transformations}
\item{tr_bpar}{The estimated parameter in the Box-Cox transformation}
\item{tr_y}{Vector of the transformed explaied variables}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outputs, which are internally used}
}
\references{
Murakami, D., Kajita, M., Kajita, S. and Matsui, T. (2020) Compositionally-warped additive mixed modeling for a wide variety of non-Gaussian data, Arxiv.

Rios, G. and Tobar, F. (2019) Compositionally-warped Gaussian processes. Neural Networks, 118, 235-246.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

Griffith, D. A. (2003). Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}, \code{\link{besf}}
}
\examples{
require(spdep);require(Matrix)
data(boston)
y	    <- boston.c[, "CMEDV" ]
x	    <- boston.c[,c("CRIM","ZN","INDUS", "CHAS", "NOX","RM", "AGE",
                     "DIS" ,"RAD", "TAX", "PTRATIO", "B", "LSTAT")]
xgroup<- boston.c[,"TOWN"]
coords<- boston.c[,c("LON","LAT")]
meig 	<- meigen(coords=coords)
# meig<- meigen_f(coords=coords)  ## for large samples

######## Regression considering residual spatially dependence
res	  <- resf(y = y, x = x, meig = meig)
res
plot_s(res)    ## spatially dependent component (intercept)

######## Compositionally-warped spatial regression (2 SAL transformations)
res2	  <- resf(y = y, x = x, meig = meig, tr_num = 2)
res2               ## tr_num and tr_nonneg can be selected by comparing BIC (or AIC)
coef_marginal(res2)## marginal effects of x. The median might be useful as a summary statistic

######## Compositionally-warped spatial regression (2 SAL trans. + Box-Cox trans.)
res3	  <- resf(y = y, x = x, meig = meig, tr_num = 2, tr_nonneg=TRUE)
res3            ## tr_num and tr_nonneg can be selected by comparing BIC (or AIC)
coef_marginal(res3)

######## Regression considering residual spatially dependence and NVC
######## (constant coefficients or NVC is selected)
#res4 <- resf(y = y, x = x, meig = meig, nvc = TRUE)
#res4          ## Note: Coefficients on 5,6,and 13-th covariates
               ## are estimated non-spatially varying (NVC) depending on x

#plot_n(res4,5) ## 1D plot of the 5-th NVC
#plot_n(res4,6) ## 1D plot of the 6-th NVC
#plot_n(res4,13)## 1D plot of the 13-th NVC

#plot_s(res4)   ## spatially dependent component (intercept)
#plot_s(res4,5) ## spatial plot of the 5-th NVC
#plot_s(res4,6) ## spatial plot of the 6-th NVC
#plot_s(res4,13)## spatial plot of the 13-th NVC

######## Compositionally-warped spatial regression with NVC (2 SAL trans. + Box-Cox trans.)
######## (constant coefficients or NVC is selected)
#res5 <- resf(y = y, x = x, meig = meig, nvc = TRUE, tr_num = 2, tr_nonneg=TRUE)

######## Regression considering residual spatially dependence and NVC
######## (all the coefficients are NVCs)
#res6 <- resf(y = y, x = x, meig = meig, nvc = TRUE, nvc_sel=FALSE)

######## Regression considering residual spatially dependence and group effects
#res7 <- resf(y = y, x = x, meig = meig, xgroup = xgroup)

######## Regression considering group-level spatially dependence and group effects
#meig_g<- meigen(coords=coords, s_id = xgroup)
#res8 <- resf(y = y, x = x, meig = meig_g, xgroup = xgroup)

}

