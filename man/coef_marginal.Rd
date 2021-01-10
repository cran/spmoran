\name{coef_marginal}
\alias{coef_marginal}
\title{Marginal effects evaluation}
\usage{
coef_marginal( mod )
}
\description{
This function evaluates the marginal effects of x based on the estimation result of \code{\link{resf}}. This funtion works if y is transformed  (i.e., tr_num > 0 or tr_nonneg =TRUE) in \code{\link{resf}}.
}
\arguments{
\item{mod}{Output from \code{\link{resf}}}
}
\value{
\item{b}{Marginal effects of x}
}
\seealso{
\code{\link{resf}}
}
