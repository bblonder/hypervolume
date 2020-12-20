\name{hypervolume_thin}
\alias{hypervolume_thin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reduces the number of random points in a hypervolume
}
\description{
Many hypervolume algorithms have computational complexities that scale with the number of random points used to characterize a hypervolume (\code{@RandomPoints}). This value can be reduced to improve runtimes at the cost of lower resolution.
}
\usage{
hypervolume_thin(hv, factor = NULL, num.points = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An object of class \code{Hypervolume}
}
  \item{factor}{
A number in (0,1) describing the fraction of random points to keep.
}
  \item{num.points}{
A number describing the number random points to keep.
}
}
\details{
Either \code{factor} or \code{npoints} (but not both) must be specified.
}
\value{
A \code{Hypervolume} object
}
\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
  
hv = hypervolume_box(penguins_adelie,name='Adelie')

# downsample to 1000 random points
hv_thinned = hypervolume_thin(hv, num.points=1000)
hv_thinned
}