\name{coef_marginal_vc}
\alias{coef_marginal_vc}
\title{Marginal effects evaluation from models with varying coefficients}
\usage{
coef_marginal_vc( mod )
}
\description{
This function evaluates the marginal effects of x based on the estimation result of \code{\link{resf_vc}}. This funtion works if y is transformed  (i.e., tr_num > 0 or tr_nonneg =TRUE) in \code{\link{resf}}.
}
\arguments{
\item{mod}{Output from \code{\link{resf}}}
}
\value{
\item{b_vc}{Matrix of the marginal effects of x (N x K)}
\item{B_vc_n}{Matrix of the sub-marginal effects of x explained by the spatially varying coefficients (N x K)}
\item{B_vc_s}{Matrix of the sub-marginal effects explained by the non-spatially varying coefficients (N x K)}
\item{c}{Matrix of the marginal effects of xconst (N x K_const)}
\item{other}{List of other outputs, which are internally used}
}
\seealso{
\code{\link{resf_vc}}
}
