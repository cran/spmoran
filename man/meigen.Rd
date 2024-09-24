\name{meigen}
\alias{meigen}
\title{Extraction of Moran eigenvectors}
\usage{
meigen( coords = NULL,model = "exp", enum = NULL, s_id = NULL,
        threshold = 0, cmat = NULL, coords_z=NULL, enum_z=NULL,
        interact=TRUE, interact_max_dim = 600  )
}
\description{
This function extracts spatial and temporal eigenvectors (i.e., basis functions describing spatial and temporal patterns).
}
\arguments{
\item{coords}{Matrix of spatial coordinates (N x 2). If cmat is specified, it is ignored}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel. Default is "exp"}
\item{enum}{Optional. The maximum mumber of spatial eigenvectors to be extracted (scalar)}
\item{s_id}{Optional. Location/zone ID for modeling inter-group spatial effects. If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
\item{threshold}{Optional. Threshold for the eigenvalues. Suppose that lambda_1 is the maximum eigenvalue, this function extracts eigenvectors whose corresponding eigenvalue is equal or greater than (threshold x lambda_1). threshold must be a value between 0 and 1. Default is zero (see Details)}
\item{cmat}{Optional. A user-specified spatial connectivity matrix (N x N). It must be provided when the user wants to use a spatial connectivity matrix other than the default matrices}
\item{coords_z}{Optional. One- or two-column matrix whose t-th column represents t-th temporal coordinate (N x 1 or N x 2).}
\item{enum_z}{Optional. The maximum mumber of temporal eigenvectors to be extracted (scalar)}
\item{interact}{Optional. If TRUE, space-time eigenvectors (space x time) are considered in addition to spatial eigenvectors and temporal eigenvectors}
\item{interact_max_dim}{Optional. The maximum mumber of the space-time eigenvectors to be extracted (scalar)}
}
\details{
This function extracts spatial eigenvectors from MCM, where M = I - 11'/N is a centering operator. By default, C is a N x N connectivity matrix whose (i, j)-th element equals exp(-d(i,j)/h), where d(i,j) is the spatial Euclidean distance between the sample sites i and j. h is the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006). If cmat is provided, this function performs the same calculation after C is replaced with cmat.

The temporal eigenvectors are extracted in the same way where the spatial distance d(i,j) is replaced with temporal difference. If two temporal coordinates are given, their eigenvectors are evaluated respectively.

If threshold = 0.00 (default), all the eigenvectors corresponding to positive eigenvalues explaining positive spatial/temporal dependence are extracted. threshold = 0.00 or 0.25 are standard assumptions (see Griffith, 2003; Murakami and Griffith, 2015).
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

Griffith, D.A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

Murakami, D., Shirota, S., Kajita, S., and Kajita, S. (2024) Fast spatio-temporally varying coefficient modeling with reluctant interaction selection. ArXiv.
}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen_f}} for fast eigen-decomposition
}
