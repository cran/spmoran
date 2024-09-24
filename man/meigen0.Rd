\name{meigen0}
\alias{meigen0}
\title{Nystrom extension of Moran eigenvectors}
\usage{
meigen0( meig, coords0, coords_z0 = NULL, s_id0 = NULL )
}
\description{
This function estimates Moran eigenvectors at unobserved sites using the Nystrom extension.
}
\arguments{
\item{meig}{Moran eigenvectors and eigenvalues. Output from \code{\link{meigen}} or \code{\link{meigen_f}}}
\item{coords0}{Matrix of spatial point coordinates of prediction sites (N_0 x 2)}
\item{coords_z0}{Optional. One- or two-column matrix whose t-th column represents the t-th temporal coordinate of prediction times (N_0 x 1 or N_0 x 2).}
\item{s_id0}{Optional. ID specifying groups modeling spatial effects (N_0 x 1). If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
}
\value{
\item{sf}{Matrix of the first L eigenvectors at unobserved sites (N_0 x L)}
\item{ev}{Vector of the first L eigenvalues (L x 1)}
\item{sf_z}{List. t-th element is the matrix of the t-th temporal eigenvectors (N x L_t)}
\item{ev_z}{List. t-th element is the vector of the t-th temporal eigenvalues (L_t x 1)}
\item{other}{List of other outputs, which are internally used}
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
