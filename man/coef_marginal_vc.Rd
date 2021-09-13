\name{coef_marginal_vc}
\alias{coef_marginal_vc}
\title{Marginal effects evaluation from models with varying coefficients}
\usage{
coef_marginal_vc( mod )
}
\description{
This function evaluates the marginal effects from x (dy/dx) based on the estimation result of \code{\link{resf_vc}}. This funtion is for non-Gaussian models transforming y using \code{\link{nongauss_y}}.
}
\arguments{
\item{mod}{Output from \code{\link{resf_vc}}}
}
\value{
\item{b_vc}{Matrix of the marginal effects of x (dy/dx) (N x K)}
\item{B_vc_n}{Matrix of the sub-marginal effects of x explained by the spatially varying coefficients (N x K)}
\item{B_vc_s}{Matrix of the sub-marginal effects explained by the non-spatially varying coefficients (N x K)}
\item{c}{Matrix of the marginal effects of xconst (N x K_const)}
\item{other}{List of other outputs, which are internally used}
}
\seealso{
\code{\link{resf_vc}}
}
