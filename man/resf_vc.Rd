\name{resf_vc}
\alias{resf_vc}
\title{Spatially and non-spatially varying coefficient (SNVC) modeling for Gaussian and non-Gaussian continuous data}

\usage{
resf_vc(y, x, xconst = NULL, xgroup = NULL, weight = NULL,
        x_nvc = FALSE, xconst_nvc = FALSE, x_sel = TRUE, x_nvc_sel = TRUE,
        xconst_nvc_sel = TRUE, nvc_num = 5, meig, method = "reml",
        penalty = "bic", maxiter = 30, tr_num = 0, tr_nonneg = FALSE )
}
\description{
The model estimates residual spatial dependence, constant coefficients, spatially varying coefficients (SVCs), non-spatially varying coefficients (NVC; coefficients varying with respect to explanatory variable value), SNVC (= SVC + NVC), and group effects. The random-effects eigenvector spatial filtering (RE-ESF), which is a low rank Gaussian process approach, is used for the spatial process modeling. While the resf_vc function estimates a SVC model by default, the type of coefficients (constant, SVC, NVC, or SNVC) can be selected through a BIC/AIC minimization. In addition, to transform non-Gaussian explained variables to Gaussian variables, the compositionally-warping function and/or the Box-Cox transformation is available (see Murakami et al., 2020).

Note: SNVCs can be mapped just like SVCs. SNVC model is more robust against spurious correlation (multicollinearity) and stable than SVC models (see Murakami and Griffith, 2020).
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables with spatially varying coefficients (SVC) (N x K)}
\item{xconst}{Matrix of explanatory variables with constant coefficients (N x K_{c}). Default is NULL}
\item{xgroup}{Matrix of group IDs. The IDs may be group numbers or group names (N x K_{g}). Default is NULL}
\item{weight}{Vector of weights for samples (N x 1). When non-NULL, the adjusted R-squared value is evaluated for weighted explained variables. Default is NULL}
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
\item{tr_num}{Number of the SAL transformations ( SinhArcsinh and Affine, where the use of "L" stems from the "Linear") used to transform non-Gaussian explained variables to Gaussian variables. Default is 0}
\item{tr_nonneg}{If TRUE, the Box-Cox transfromation used to transform positive non-Gaussian explained variables to Gaussian variables. If tr_num > 0 and tr_nonneg == TRUE, the Box-Cox transformation is applied first. Then, th SAL transformation is applied tr_num times. Default is FALSE}
}
\details{
If tr_num >0, the resf function iterates the SAL transformation tr_num times to transform the explained variables to Gaussian variables. The SAL transformation is defined as SAL(y)=a+b*sinh(c*arcsinh(y)-d) where {a,b,c,d} are parameters. Based on Rois and Tober (2019), an iteration of the SAL transformation approximates a wide variety of non-Gaussian distributions without explicitly assuming data distribution. As a result, our spatial regression approach is applicable to a wide variety of non-Gaussian continuous data too. For non-negative explained variables, a Box-Cox transformation is available prior to the SAL transformations by specifying tr_nonneg >0. tr_num and tr_nonneg can be selected by comparing the BIC (or AIC) value across models. This compositionally-warped spatial regression approach is detailed in Murakami et al. (2020).
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
\item{e}{Vector whose elements are residual standard error (resid_SE), adjusted conditional R2 (adjR2(cond)), restricted log-likelihood (rlogLik), Akaike information criterion (AIC), and Bayesian information criterion (BIC). When method = "ml", restricted log-likelihood (rlogLik) is replaced with log-likelihood (logLik)}
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

Murakami, D., and Griffith, D.A. (2020) Balancing spatial and non-spatial variations in varying coefficient modeling: a remedy for spurious correlation. ArXiv.

Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.

Griffith, D. A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}, \code{\link{besf_vc}}
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

############## SVC modeling1 #################
######## - SVC or constant coefficients on x
######## - Constant coefficients on xconst
res	    <- resf_vc(y=y,x=x,xconst=xconst,meig=meig )
res

plot_s(res,0) # Spatially varying intercept
plot_s(res,1) # 1st SVC (Not shown because the SVC is estimated constant)
plot_s(res,2) # 2nd SVC

############## Compositionally-warped SVC modeling #################
######## - SVC or constant coefficients on x
######## - Constant coefficients on xconst
######## - 2 SAL transformations on y

#res2	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, tr_num = 2 )
#res2                   # tr_num and tr_nonneg can be selected by comparing BIC (or AIC)
#coef_marginal_vc(res2) # marginal effects of x.
                        # The median might be useful as a summary statistic

######## - 2 SAL transformations + Box-Cox transformation on y

#res3	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, tr_num = 2, tr_nonneg=TRUE )
#res3
#coef_marginal_vc(res3)

############## SVC modeling2 #################
######## (SVC on x; Constant coefficients on xconst)
#res4	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_sel = FALSE )

############## SVC modeling3 #################
######## - Group-level SVC or constant coefficients on x
######## - Constant coefficients on xconst
######## - Group effects

#meig_g <- meigen(coords, s_id=xgroup)
#res5	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig_g,xgroup=xgroup)

############## SNVC modeling1 #################
######## - SNVC, SVC, NVC, or constant coefficients on x
######## - Constant coefficients on xconst

#res6	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE)


############## Compositionally-warped SNVC modeling #################
######## - SNVC, SVC, NVC, or constant coefficients on x
######## - 2 SAL transformations + Box-Cox transformation on y

#res7	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE, tr_num = 2, tr_nonneg=TRUE)
#res7       # tr_num and tr_nonneg can be selected by comparing BIC (or AIC)


############## SNVC modeling2 #################
######## - SNVC, SVC, NVC, or constant coefficients on x
######## - NVC or Constant coefficients on xconst

#res8	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE, xconst_nvc=TRUE)
#plot_s(res8,0)            # Spatially varying intercept
#plot_s(res8,1)            # Spatial plot of 1st SNVC (SVC + NVC)
#plot_s(res8,1,btype="svc")# Spatial plot of SVC in the SNVC on x[,1]
#plot_s(res8,1,btype="nvc")# Spatial plot of NVC in the SNVC on x[,1]
#plot_n(res8,1)            # 1D plot of NVC in the SNVC on x[,1]

#plot_s(res8,6,xtype="xconst")# Spatial plot of NVC in the SNVC on xconst[,6]
#plot_n(res8,6,xtype="xconst")# 1D plot of NVC on xconst[,6]

}
