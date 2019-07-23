\name{predict0}
\alias{predict0}
\title{Spatial prediction using eigenvector spatial filtering (ESF) or random effects ESF}
\usage{
predict0( mod, meig0, x0 = NULL, xgroup0 = NULL )
}
\description{
This function predicts explained variables using eigenvector spatial filtering (ESF) or random effects ESF. The Nystrom extension is used to minimize the expected prediction error
}
\arguments{
\item{mod}{ESF or RE-ESF model estimates. Output from \code{\link{esf}} or \code{\link{resf}}}
\item{meig0}{Moran's eigenvectors at predicted sites. Output from \code{\link{meigen0}}}
\item{x0}{Matrix of explanatory variables at predicted sites (N_0 x K). Default is NULL}
\item{xgroup0}{Matrix of group indeces that may be group IDs (integers) or group names (N_0 x K_{group}). Default is NULL}
}
\value{
\item{pred}{Matrix with the first column for the predicted values. The second and the third columns are the trend component and the spatial component in the predicted values. The remining columns are group effects (if any)}
}
\references{
Drineas, P. and Mahoney, M.W. (2005) On the Nystrom method for approximating a gram matrix for improved kernel-based learning. Journal of Machine Learning Research, 6 (2005), 2153-2175.
}
\seealso{
\code{\link{meigen0}}, \code{\link{predict0_vc}}
}
\examples{
require(spdep)
data(boston)
samp    <- sample( dim( boston.c )[ 1 ], 400)

d       <- boston.c[ samp, ]  ## Data at observed sites
y	      <- d[, "CMEDV"]
x       <- d[,c("ZN","INDUS", "NOX","RM", "AGE", "DIS")]
coords  <- d[,c("LAT","LON")]

d0      <- boston.c[-samp, ]  ## Data at unobserved sites
y0	    <- d0[, "CMEDV"]
x0      <- d0[,c("ZN","INDUS", "NOX","RM", "AGE", "DIS")]
coords0 <- d0[,c("LAT","LON")]

############ Model estimation
meig 	  <- meigen( coords = coords )
mod	    <- resf(y=y, x=x, meig=meig)
## or
# mod   <- esf(y=y,x=x,meig=meig)

############ Spatial prediction

meig0 	<- meigen0( meig = meig, coords0 = coords0 )
pred0   <- predict0( mod = mod, x0 = x0, meig0 = meig0 )
pred0[1:10,]
}
