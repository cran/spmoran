\name{meigen_f}
\alias{meigen_f}
\title{Fast approximation of Moran eigenvectors}
\usage{
meigen_f( coords, model = "exp", enum = 200, s_id = NULL, threshold = 0,
          coords_z = NULL, enum_z = 200, interact = TRUE,
          interact_max_dim = 600 )
}
\description{
This function approximates spatial and temporal eigenvectors (i.e., basis functions describing spatial and temporal patterns) computationally efficiently.
}
\arguments{
\item{coords}{Matrix of spatial coordinates (N x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel. Default is "exp"}
\item{enum}{Number of eigenvectors and eigenvalues to be extracted (scalar). Default is 200}
\item{s_id}{Optional. Location/zone ID for modeling inter-group spatial effects. If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
\item{threshold}{Optional. Threshold for the eigenvalues. Suppose that lambda_1 is the maximum eigenvalue, this function extracts eigenvectors whose corresponding eigenvalue is equal or greater than (threshold x lambda_1). threshold must be a value between 0 and 1. Default is zero}
\item{coords_z}{Optional. One- or two-column matrix of temporal coordinates (N x 1 or N x 2).}
\item{enum_z}{Optional. The maximum mumber of temporal eigenvectors to be extracted (scalar)}
\item{interact}{Optional. If TRUE, space-time eigenvectors (space x time) are considered in addition to spatial eigenvectors and temporal eigenvectors}
\item{interact_max_dim}{Optional. The maximum mumber of the space-time eigenvectors to be extracted (scalar)}
}
\details{
This function extracts approximated spatial eigenvectors from MCM. M = I - 11'/N is a centering operator, and C is a spatial connectivity matrix whose (i, j)-th element is given by exp( -d(i,j)/h), where d(i,j) is the Euclidean distance between the sample sites i and j, and h is a range parameter given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006). Following a simulation result in Murakami and Griffith (2019), this function approximates the 200 eigenvectors corresponding to the 200 largest eigenvalues by default (i.e., enum = 200). If enum is given by a smaller value like 100, the computation time will be shorter, but with greater approximation error.

The temporal eigenvectors are extracted in the same way where the spatial distance d(i,j) is replaced with temporal difference. If two temporal coordinates are given, their eigenvectors are evaluated respectively.
}
\value{
\item{sf}{Matrix of the spatial eigenvectors (N x L)}
\item{ev}{Vector of the spatial eigenvalues (L x 1), scaled to have the maximum value of 1}
\item{sf_z}{List. t-th element is the matrix of the t-th temporal eigenvectors (N x L_t)}
\item{ev_z}{List. t-th element is the vector of the t-th temporal eigenvalues (L_t x 1), scaled to have the maximum value of 1}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Dray, S., Legendre, P., and Peres-Neto, P.R. (2006) Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196 (3), 483-493.

Murakami, D. and Griffith, D.A. (2019) Eigenvector spatial filtering for large data sets: fixed and random effects approaches. Geographical Analysis, 51 (1), 23-49.

Murakami, D., Shirota, S., Kajita, S., and Kajita, S. (2024) Fast spatio-temporally varying coefficient modeling with reluctant interaction selection. ArXiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}
}
