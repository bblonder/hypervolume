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
   show.3d=FALSE,plot.3d.axes.id=NULL,
   show.axes=TRUE, show.frame=TRUE,
   show.random=TRUE, show.density=TRUE,show.data=TRUE,
   names=NULL, show.legend=TRUE, limits=NULL, 
   show.contour=TRUE, contour.lwd=1, 
    contour.filled=FALSE,contour.filled.alpha=0.5,contour.factor=0.5,
   show.centroid=TRUE, cex.centroid=3,
   colors=rainbow(floor(length(x@HVList)*1.5),alpha=0.8), 
   point.alpha.min=0.2, point.dark.factor=0.5,
   cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
   num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE, 
   verbose=FALSE,
   ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A Hypervolume or HypervolumeList object. The objects to be plotted.
}
\item{show.3d}{
If \code{TRUE}, makes a three-dimensional plot of a subset of axes determined by \code{plot.3d.axes.id}; otherwise, a pairs plot of all axes.
}
\item{plot.3d.axes.id}{
Numeric identities of axes to plot in three dimensions. Defaults to 1:3 if set to \code{NULL}.
}
\item{show.axes}{
If \code{TRUE}, draws axes on the plot.
}
\item{show.frame}{
If \code{TRUE}, frames the plot with a box.
}
\item{show.random}{
If \code{TRUE}, shows random points from the hypervolume.
}
\item{show.density}{
If \code{TRUE}, draws random points with alpha level proportional to their unit-scaled probability density. Note that this has no effect when probability density is not relevant, i.e. for hypervolumes that are the output of set operations.
}
\item{show.data}{
If \code{TRUE}, draws data points from the hypervolume. Note that this has no effect if the hypervolume is not associated with data points, e.g. for those that are the output of set operations.
}
  \item{names}{
A vector of strings in the same order as the input hypervolumes. Used to draw the axes labels.
}
  \item{show.legend}{
If \code{TRUE}, draws a color legend.
}
  \item{limits}{
A list of two-element vectors corresponding to the axes limits for each dimension. If a single two-element vector is provided it is re-used for all axes.
}
\item{show.contour}{
If \code{TRUE}, draws a boundary line saround each two-dimensional projection. Ignored if \code{show.3d=TRUE}.}

\item{contour.lwd}{
Line width used for contour lines. Ignored if \code{show.contour=FALSE}.
}
\item{contour.filled}{
If \code{TRUE}, draws contours with interior shading. Ignored if \code{show.contour=FALSE}.
}
\item{contour.filled.alpha}{
An alpha value used to determine the transparency of shaded contours. Ignored if \code{contour.filled=FALSE}.
}
\item{contour.factor}{
A numeric value that controls the threshold level at which the contour is drawn. Scaled so that 0 includes all data and 1 includes no data. Ignored if \code{showcontour=FALSE}.
}
\item{show.centroid}{
If \code{TRUE}, draws a colored point indicating the centroid for each hypervolume.
}
\item{cex.centroid}{
Expansion factor for the centroid symbol.
}
  \item{colors}{
A vector of colors to be used to plot each hypervolume, in the same order as the input hypervolumes.
}
  \item{point.alpha.min}{
Fractional value corresponding to the most transparent value for plotting random points. 0 corresponds to full transparency.
}
  \item{point.dark.factor}{
Fractional value corresponding to the darkening factor for plotting data points. 0 corresponds to fully black.
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

  \item{num.points.max.data}{
An integer indicating the maximum number of data points to be sampled from each hypervolume. Lower values result in faster plotting and smaller file sizes but less accuracy.
}
  \item{num.points.max.random}{
An integer indicating the maximum number of random points to be sampled from each hypervolume. Lower values result in faster plotting and smaller file sizes but less accuracy.
}
  \item{reshuffle}{
A logical value relevant when pair=TRUE. If false, each hypervolume is drawn on top of the previous hypervolume; if true, all points of all hypervolumes are randomly shuffled so no hypervolume is given visual preference during plotting.
}
\item{verbose}{
If \code{TRUE}, prints diagnostic information about the number of points being plotted
}




\item{...}{
Additional arguments to be passed to \code{rgl::plot3d}. 
}
}
\value{
None; used for the side-effect of producing a plot.
}
\note{
Note that contour lines are drawn only for convenience and visual presentation. They are calculated using \code{MASS::kde2d} using a hexagonal binning of projected data with a five percent quantile threshold. The true hypervolume boundary (on which volumes and set operations are calculated) is reflected by the less smooth uniformly random set of points drawn by \code{showrandom=TRUE}. 
}
\seealso{
\code{\link{hypervolume_save_animated_gif}}
}
\examples{
data(iris)
hv = hypervolume_gaussian(subset(iris, Species=="setosa")[,1:3])

plot(hv, show.3d=TRUE)
plot(hv, show.3d=FALSE)
}