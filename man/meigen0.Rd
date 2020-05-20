\name{meigen0}
\alias{meigen0}
\title{Nystrom extension of Moran eigenvectors}
\usage{
meigen0( meig, coords0, s_id0 = NULL )
}
\description{
This function estimates Moran eigenvectors at unobserved sites using the Nystrom extension.
}
\arguments{
\item{coords0}{Matrix of spatial point coordinates of unobserved sites (N_0 x 2)}
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{s_id0}{Optional. ID specifying groups modeling spatial effects (N_0 x 1). If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
}
\value{
\item{sf}{Matrix of the first L eigenvectors at unobserved sites (N_0 x L)}
\item{ev}{Vector of the first L eigenvalues (L x 1)}
\item{ev_full}{Vector of all eigenvalues (N x 1)}
}
\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}
}
