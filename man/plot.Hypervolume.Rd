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
   show.contour=TRUE, contour.lwd=1.5, 
    contour.type='kde', 
    contour.alphahull.alpha=0.25,
    contour.ball.radius.factor=1, 
    contour.kde.level=0.01,
    contour.raster.resolution=100,
   show.centroid=TRUE, cex.centroid=2,
   colors=rainbow(floor(length(x@HVList)*1.5),alpha=0.8), 
   point.alpha.min=0.2, point.dark.factor=0.5,
   cex.random=0.5,cex.data=0.75,cex.axis=0.75,cex.names=1.0,cex.legend=0.75,
   num.points.max.data = 1000, num.points.max.random = 2000, reshuffle=TRUE,
   plot.function.additional=NULL,
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
\item{contour.type}{
Type of contour boundary: any of \code{"alphahull"} (alpha hull), \code{"ball"} (experimental ball covering), \code{"kde"} (2D KDE smoothing), or \code{"raster"} (grid-based rasterization).
}

\item{contour.alphahull.alpha}{
Value of the alpha parameter for a  \code{"alphahull"} contour. Can be increased to provide smoother contours.
}
\item{contour.ball.radius.factor}{
Factor used to multiply radius of ball surrounding each random point for a \code{"ball"} contour.
}
\item{contour.kde.level}{
Probability level used to delineate edges for a \code{"kde"} contour.
}
\item{contour.raster.resolution}{
Grid resolution for a \code{"raster"} contour.
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
\item{plot.function.additional}{
Any \code{function(i,j)} that will add additional plotting commands for column \code{i} and row \code{j} of the pairs plot. Should not create new plots or change \code{par()} settings. Has no effect if \code{show.3d=TRUE}.

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
Contour line plotting with \code{alphahull} requires the non-FOSS \code{alphahull} package to be installed. Please do so in order to use this functionality!
}
\seealso{
\code{\link{hypervolume_save_animated_gif}}
}
\examples{
\dontrun{
data(iris)
hv = hypervolume_gaussian(iris[,1:3])

plot(hv, show.3d=TRUE)
plot(hv, show.3d=FALSE)

plot(hv,plot.function.additional=function(i,j){
  points( x=iris[iris$Species=="setosa",i],
          y=iris[iris$Species=="setosa",j],col='purple') 
  })
}
}