\name{resf_vc}
\alias{resf_vc}
\title{Spatially and non-spatially varying coefficient (SNVC) modeling}
\usage{
resf_vc(y, x, xconst = NULL, xgroup = NULL, x_nvc = FALSE, xconst_nvc = FALSE,
        x_sel = TRUE, x_nvc_sel = TRUE, xconst_nvc_sel = TRUE, nvc_num = 5,
        meig, method = "reml", penalty = "bic", maxiter = 30 )
}
\description{
The model estimates residual spatial dependence, constant coefficients, spatially varying coefficients (SVCs), non-spatially varying coefficients (NVC; coefficients varying with respect to explanatory variable value), SNVC (= SVC + NVC), and group effects. Type of coefficients can be selected through BIC/AIC minimization. By default, it estimates a SVC model.

Note: SNVCs can be mapped just like SVCs. Unlike SVC models, SNVC model is robust against spurious correlation (multicollinearity), so, stable (see Murakami and Griffith, 2020).
}
\arguments{
\item{y}{Vector of explained variables (N x 1)}
\item{x}{Matrix of explanatory variables with spatially varying coefficients (SVC) (N x K)}
\item{xconst}{Matrix of explanatory variables with constant coefficients (N x K_{c}). Default is NULL}
\item{xgroup}{Matrix of group IDs. The IDs may be group numbers or group names (N x K_{g}). Default is NULL}
\item{x_nvc}{If TRUE, SNVCs are assumed on x. Otherwise, SVCs are assumed. Default is FALSE}
\item{xconst_nvc}{If TRUE, NVCs are assumed on xconst. Otherwise, constant coefficients are assumed. Default is FALSE}
\item{x_sel}{If TRUE, type of coefficient (SVC or constant) on x is selected through a BIC (default) or AIC minimization. If FALSE, SVCs are assumed across x. Alternatively, x_sel can be given by column number(s) of x. For example, if x_sel = 2, the coefficient on the second explanatory variable in x is SVC and the other coefficients are constants. The Default is TRUE}
\item{x_nvc_sel}{If TRUE, type of coefficient (NVC or constant) on x is selected through the BIC (default) or AIC minimization. If FALSE, NVCs are assumed across x. Alternatively, x_nvc_sel can be given by column number(s) of x. For example, if x_nvc_sel = 2, the coefficient on the second explanatory variable in x is NVC and the other coefficients are constants. The Default is TRUE}
\item{xconst_nvc_sel}{If TRUE, type of coefficient (NVC or constant) on xconst is selected through the BIC (default) or AIC minimization. If FALSE, NVCs are assumed across xconst. Alternatively, xconst_nvc_sel can be given by column number(s) of xconst. For example, if xconst_nvc_sel = 2, the coefficient on the second explanatory variable in xconst is NVC and the other coefficients are constants.The Default is TRUE}
\item{nvc_num}{Number of basis functions used to model NVC. An intercept and nvc_num natural spline basis functions are used to model each NVC. Default is 5}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{method}{Estimation method. Restricted maximum likelihood method ("reml") and maximum likelihood method ("ml") are available. Default is "reml"}
\item{penalty}{Penalty to select varying coefficients and stablize the estimates. The current options are "bic" for the Baysian information criterion-type penalty (N x log(K)) and "aic" for the Akaike information criterion  (2K) (see Muller et al., 2013). Default is "bic"}
\item{maxiter}{Maximum number of iterations. Default is 30}
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
\item{pred}{Vector of predicted values (N x 1)}
\item{resid}{Vector of residuals (N x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{

Murakami, D., Yoshida, T., Seya, H., Griffith, D.A., and Yamagata, Y. (2017) A Moran coefficient-based mixed effects approach to investigate spatially varying relationships. Spatial Statistics, 19, 68-89.

Murakami, D., and Griffith, D.A. (2019) Spatially varying coefficient modeling for large datasets: Eliminating N from spatial regressions. Spatial Statistics, 30, 39-64.

Murakami, D., and Griffith, D.A. (2020) Balancing spatial and non-spatial variations in varying coefficient modeling: a remedy for spurious correlation. ArXiv.

Griffith, D. A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Muller, S., Scealy, J.L., and Welsh, A.H. (2013) Model selection in linear mixed models. Statistical Science, 28 (2), 136-167.
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
coords  <- boston.c[,c("LAT","LON")]
meig 	  <- meigen(coords=coords)
# meig	<- meigen_f(coords=coords)  ## for large samples

############## SVC modeling1 #################
######## (SVC on x; Constant coefficients on xconst)
res	    <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_sel = FALSE )
res

plot_s(res,0) # Spatially varying intercept
plot_s(res,1) # 1st SVC
plot_s(res,2) # 2nd SVC

############## SVC modeling2 #################
######## (SVC or constant coefficients on x; Constant coefficients on xconst)
#res2	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig )

############## SVC modeling3 #################
######## - Group-level SVC or constant coefficients on x
######## - Constant coefficients on xconst
######## - Group effects

#meig_g <- meigen(coords, s_id=xgroup)
#res3	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig_g,xgroup=xgroup)

############## SNVC modeling1 #################
######## - SNVC, SVC, NVC, or constant coefficients on x
######## - Constant coefficients on xconst

#res4	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE)

############## SNVC modeling2 #################
######## - SNVC, SVC, NVC, or constant coefficients on x
######## - NVC or Constant coefficients on xconst

#res5	  <- resf_vc(y=y,x=x,xconst=xconst,meig=meig, x_nvc =TRUE, xconst_nvc=TRUE)
#plot_s(res5,0)           # Spatially varying intercept
#plot_s(res5,1)           # 1st SNVC
#plot_s(res5,1,snvc=FALSE)# SVC in the 1st SNVC
#plot_n(res5,1,xtype="x") # NVC in the 1st NVC
#plot_n(res5,6,xtype="xconst")

}
