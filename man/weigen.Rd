\name{weigen}
\alias{weigen}
\title{Extract eigenvectors from a spatial weight matrix}
\usage{
weigen( x = NULL, type = "knn", k = 4, threshold = 0.25, enum = NULL )
}
\description{
This function extracts eigenvectors and eigenvalues from a spatial weight matrix.
}
\arguments{
\item{x}{Matrix of spatial point coordinates (N x 2), ShapePolygons object (N spatial units), or an user-specified spatial weight matrix (N x N) (see Details)}
\item{type}{Type of spatial weights. The currently available options are "knn" for the k-nearest neighbor-based weights, and "tri" for the Delaunay triangulation-based weights. If ShapePolygons are provided for x, type is ignored, and the rook-type neighborhood matrix is created}
\item{k}{Number of nearest neighbors. It is used if type ="knn"}
\item{threshold}{Threshold for the eigenvalues (scalar). Suppose that lambda_1 is the maximum eigenvalue. Then, this fucntion extracts eigenvectors whose corresponding eigenvalues are equal or greater than [threshold x lambda_1]. It must be a value between 0 and 1. Default is 0.25 (see Details)}
\item{enum}{Optional. The muximum acceptable mumber of eigenvectors to be used for spatial modeling (scalar)}
}
\details{
If user-specified spatial weight matrix is provided for x, this function returns the eigen-pairs of the matrix. Otherwise, if a SpatialPolygons object is provided to x, the rook-type neighborhood matrix is created using this polygon, and eigen-decomposed. Otherwise, if point coordinats are provided to x, a spatial weight matrix is created according to type, and eigen-decomposed.

By default, the ARPACK routine is implemented for fast eigen-decomposition.

threshold = 0.25 (default) is a standard setting for topology-based ESF (see Tiefelsdorf and Griffith, 2007) while threshold = 0.00 is a usual setting for distance-based ESF.
}
\value{
\item{sf}{Matrix of the first L eigenvectors (N x L)}
\item{ev}{Vector of the first L eigenvalues (L x 1)}
\item{other}{List of other outcomes, which are internally used}
}
\references{
Tiefelsdorf, M. and Griffith, D.A. (2007) Semiparametric filtering of spatial autocorrelation: the eigenvector approach. Environment and Planning A, 39 (5), 1193-1221.

Murakami, D. and Griffith, D.A. (2018) Low rank spatial econometric models. Arxiv, 1810.02956.

}
\author{
Daisuke Murakami
}
\seealso{
\code{\link{meigen}}, \code{\link{meigen_f}}
}
\examples{
require(spdep);library(rgdal)
data(boston)

########## Rook adjacency-based W
poly	    <- readOGR(system.file("shapes/boston_tracts.shp",package="spData")[1])
weig1	    <- weigen( poly )

########## knn-based W
coords    <- boston.c[,c("LAT", "LON")]
weig2	    <- weigen( coords, type = "knn" )

########## Delaunay triangulation-based W
coords    <- boston.c[,c("LAT", "LON")]
weig3	    <- weigen( coords, type = "tri")

########## User-specified W
dmat      <- as.matrix(dist(coords))
cmat 	    <- exp(-dmat)
diag(cmat)<- 0
weig4	    <- weigen( cmat, threshold = 0 )
}
