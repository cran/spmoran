\name{meigen_f}
\alias{meigen_f}
\title{Fast approximation of Moran eigenvectors}
\usage{
meigen_f( coords, model = "exp", enum = 200, s_id = NULL )
}
\description{
This function performs a fast approximation of Moran eigenvectors and eigenvalues.
}
\arguments{
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel. Default is "exp"}
\item{enum}{Number of eigenvectors and eigenvalues to be extracted (scalar). Default is 200}
\item{s_id}{Optional. ID specifying groups modeling spatial effects (N x 1). If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
}
\details{
This function extracts approximated Moran eigenvectors from MCM. M = I - 11'/N is a centering operator, and C is a spatial connectivity matrix whose (i, j)-th element is given by exp( -d(i,j)/h), where d(i,j) is the Euclidean distance between the sample sites i and j, and h is a range parameter given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006).

Following a simulation result that 200 eigenvectors are sufficient for accurate approximation of ESF models (Murakami and Griffith, 2019), this function approximates the 200 eigenvectors corresponding to the 200 largest eigenvalues by default (i.e., enum = 200). If enum is given by a smaller value like 100, the computation time will be shorter, but with greater approximation error. Eigenvectors corresponding to negative eigenvalues are omitted from the enum eigenvectors.
}
\value{
\item{sf}{Matrix of the first L approximated eigenvectors (N x L)}
\item{ev}{Vector of the first L approximated eigenvalues (L x 1)}
\item{ev_full}{Vector of all approximated eigenvalues (enum x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Dray, S., Legendre, P., and Peres-Neto, P.R. (2006) Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196 (3), 483-493.

Murakami, D. and Griffith, D.A. (2019) Eigenvector spatial filtering for large data sets: fixed and random effects approaches. Geographical Analysis, 51 (1), 23-49.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}
}
