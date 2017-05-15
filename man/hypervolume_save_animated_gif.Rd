\name{hypervolume_save_animated_gif}
\alias{hypervolume_save_animated_gif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saves animated GIF of three-dimensional hypervolume plot.
}
\description{
Assumes there is an open RGL plot (e.g. from calling \code{plot(hv, show.3d=TRUE)}). Rotates the plot around an axis at a given speed and saves results as a series of GIFs. If the \code{magick} package is available, combines these GIFs into a single animation.
}
\usage{
hypervolume_save_animated_gif(image.size = 400, 
                              axis = c(0, 0, 1), rpm = 4, duration = 15, fps = 10, 
                              file.name = "movie", directory.output = ".", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image.size}{
Number of pixels on each side of the animated image.
}
  \item{axis}{
A three-element vector describing the rotation axis.
}
  \item{rpm}{
Animation speed in rotations per minute.
}
  \item{duration}{
Animation duration in seconds.
}
  \item{fps}{
Animation speed in frames per second.
}
  \item{file.name}{
A base name (no extension) for the GIFs.
}
  \item{directory.output}{
The folder in which output should be located.
}
  \item{\dots}{
Other arguments to be passed to \code{rgl::movie3d}.
}
}

\value{
None; used for the side-effect of producing files.
}
\examples{
# not run for speed - uncomment to try!
#data(iris)
# hv = hypervolume_gaussian(iris[,1:3])
# plot(hv, show.3d=TRUE)
# hypervolume_save_animated_gif()

}