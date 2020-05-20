\name{meigen}
\alias{meigen}
\title{Extraction of Moran's eigenvectors}
\usage{
meigen( coords, model = "exp", threshold = 0, enum = NULL, cmat = NULL, s_id = NULL )
}
\description{
This function calculates Moran eigenvectors and eigenvalues.
}
\arguments{
\item{coords}{Matrix of spatial point coordinates (N x 2)}
\item{model}{Type of kernel to model spatial dependence. The currently available options are "exp" for the exponential kernel, "gau" for the Gaussian kernel, and "sph" for the spherical kernel. Default is "exp"}
\item{threshold}{Threshold for the eigenvalues (scalar). Suppose that lambda_1 is the maximum eigenvalue, this function extracts eigenvectors whose corresponding eigenvalue is equal or greater than (threshold x lambda_1). threshold must be a value between 0 and 1. Default is zero (see Details)}
\item{enum}{Optional. The muxmum acceptable mumber of eigenvectors to be extracted (scalar)}
\item{cmat}{Optional. A user-specified spatial connectivity matrix (N x N). It must be provided when the user wants to use a spatial connectivity matrix other than the default matrices}
\item{s_id}{Optional. ID specifying groups modeling spatial effects (N x 1). If specified, Moran eigenvectors are extracted by groups. It is useful e.g. for multilevel modeling (s_id is the groups) and panel data modeling (s_id is given by individual location id). Default is NULL}
}
\details{
If cmat is not provided and model = "exp" (default), this function extracts Moran eigenvectors from MCM, where M = I - 11'/N is a centering operator. C is a N x N connectivity matrix whose (i, j)-th element equals exp(-d(i,j)/h), where d(i,j) is the Euclidean distance between the sample sites i and j, and h is given by the maximum length of the minimum spanning tree connecting sample sites (see Dray et al., 2006). If cmat is provided, this function performs the same calculation after C is replaced with cmat.

If threshold is not provided (default), all the eigenvectors corresponding to positive eigenvalue, explaining positive spatial dependence, are extracted to model positive spatial dependence. threshold = 0.00 or 0.25 are standard assumptions (see Griffith, 2003; Murakami and Griffith, 2015).
}
\value{
\item{sf}{Matrix of the first L eigenvectors (N x L)}
\item{ev}{Vector of the first L eigenvalues (L x 1)}
\item{ev_full}{Vector of all eigenvalues (N x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Dray, S., Legendre, P., and Peres-Neto, P.R. (2006) Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). Ecological Modelling, 196 (3), 483-493.

Griffith, D.A. (2003) Spatial autocorrelation and spatial filtering: gaining understanding through theory and scientific visualization. Springer Science & Business Media.

Murakami, D. and Griffith, D.A. (2015) Random effects specifications in eigenvector spatial filtering: a simulation study. Journal of Geographical Systems, 17 (4), 311-331.

}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen_f}} for fast eigen-decomposition
}
