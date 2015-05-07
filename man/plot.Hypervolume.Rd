\name{plot.HypervolumeList}
\alias{plot.Hypervolume}
\alias{plot.HypervolumeList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a hypervolume or list of hypervolumes
}
\description{
Plots a single hypervolume or multiple hypervolumes as either a pairs plot (all axes) or a 3D plot (a subset of axes). The hypervolume is drawn as a uniformly random set of points guaranteed to be in the hypervolume.
}
\usage{
\method{plot}{HypervolumeList}(x, 
npmax_data = 1000, npmax_random = 10000, 
colors=rainbow(length(x@HVList),alpha=0.8), names=NULL, 
reshuffle=TRUE, showrandom=TRUE, showdensity=TRUE,showdata=TRUE,darkfactor=0.5,
cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
legend=TRUE, varlims=NULL, 
showcontour=TRUE, contour.lwd=1, contour.filled=FALSE,contour.filled.alpha=0.5,
 showcentroid=TRUE, cex.centroid=3,
pairplot=TRUE,whichaxes=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A Hypervolume or HypervolumeList object. The objects to be plotted.
}
  \item{npmax_data}{
An integer indicating the maximum number of data points to be sampled from each hypervolume. Lower values result in faster plotting and smaller file sizes but less accuracy.
}
  \item{npmax_random}{
An integer indicating the maximum number of data points to be sampled from each hypervolume. Lower values result in faster plotting and smaller file sizes but less accuracy.
}
  \item{colors}{
A vector of colors to be used to plot each hypervolume, in the same order as the input hypervolumes.
}
  \item{names}{
A vector of strings in the same order as the input hypervolumes. Used to draw the axes labels.
}
  \item{reshuffle}{
A logical value relevant when pair=TRUE. If false, each hypervolume is drawn on top of the previous hypervolume; if true, all points of all hypervolumes are randomly shuffled so no hypervolume is given visual preference during plotting.
}
  \item{showrandom}{
A logical value indicating whether the random points comprising the hypervolume should be drawn.
}
  \item{showdensity}{
A logical value indicating whether the probability density of each hypervolume should be drawn by modulating alpha values (more transparent for lower probability density). Note that this has no effect when probability density is not relevant, i.e. for hypervolumes that are the output of set operations. Used only if \code{showrandom=TRUE}.
}
  \item{showdata}{
A logical value indicating whether the original data should be drawn on top of the uniformly random points. Note that this has no effect if the hypervolume is not associated with data points, e.g. for those that are the output of set operations.
}
  \item{darkfactor}{
A value in [0,1] that modulates the color of data points, if shown. Values closer to 0 make data points more black, while values closer to 1 make data points closer to the input color.
}
  \item{cex.random}{
cex value for uniformly random points.
}
  \item{cex.data}{
cex value for data points.
}
  \item{cex.axis}{
cex value for axes, if pair=T.
}
  \item{cex.names}{
cex value for variable names printed on the diagonal, if pair=T.
}
  \item{cex.legend}{
cex value for the legend text
}
  \item{legend}{
Logical value indicating if a legend should be plotted, if pair=T.
}
  \item{varlims}{
A list of two-element vectors corresponding to the axes limits for each dimension. If a single two-element vector is provided it is re-used for all axes.
}
\item{showcontour}{
A logical value indicating whether a boundary line should be drawn around each two-dimensional projection. Ignored if \code{pair=FALSE}.}
\item{contour.lwd}{
Line width used for contour lines. Ignored if \code{showcontour=FALSE}.
}
\item{contour.filled}{
A logical value indicating whether contours should be shaded. Ignored if \code{showcontour=FALSE}.
}
\item{contour.filled.alpha}{
An alpha value used to determine the transparency of shaded contours. Ignored if \code{contour.filled=FALSE}.
}
\item{showcentroid}{
<<<<<<< HEAD
A logical value indicating whether a solid point indicating the centroid for each hypervolume should be drawn. Centroids are randomly jittered by 1 percent to improve visual clarity.
=======
A logical value indicating whether a solid point indicating the centroid for each hypervolume should be drawn.
>>>>>>> e9d82c5f24f20034998ea2ffd2b04e1853877bc0
}
\item{cex.centroid}{
Expansion factor for the centroid symbol.
}
\item{pairplot}{
If true, a pair plot is produced. If false, a 3D plot is produced.
}
\item{whichaxes}{
A length-three vector of integer IDs corresponding to the axes to be plotted when pair=F.
}
\item{...}{
 
}
}
\value{
None; used for the side-effect of producing a plot.
}
\note{
Note that contour lines are drawn only for convenience and visual presentation. They are calculated using \code{MASS::kde2d} using projected data with a five percent quantile threshold. The true hypervolume boundary (on which volumes and set operations are calculated) is reflected by the less smooth uniformly random set of points drawn by \code{showrandom=TRUE}. 
}

\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2)

# choose fixed axes
plot(hv1, pair=TRUE, varlims=list(c(3,6),c(2,5),c(0,3),c(-1,1)))
}