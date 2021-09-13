\name{coef_marginal}
\alias{coef_marginal}
\title{Marginal effects evaluation}
\usage{
coef_marginal( mod )
}
\description{
This function evaluates the marginal effects from x (dy/dx) based on the estimation result of \code{\link{resf}}. This funtion is for non-Gaussian models transforming y using \code{\link{nongauss_y}}.
}
\arguments{
\item{mod}{Output from \code{\link{resf}}}
}
\value{
\item{b}{Marginal effects from x (dy/dx)}
}
\seealso{
\code{\link{resf}}
}
