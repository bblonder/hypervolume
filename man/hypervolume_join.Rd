\name{hypervolume_join}
\alias{hypervolume_join}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concatenate hypervolumes
}
\description{
Combines multiple hypervolumes or hypervolume lists into a single HypervolumeList suitable for analysis or plotting.
}
\usage{
hypervolume_join(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
One or more objects of class \code{Hypervolume} or \code{HypervolumeList}
}
}

\value{
A \code{HypervolumeList} containing all hypervolumes in all arguments.
}

\author{
Benjamin Blonder
}

\examples{
# UNCOMMENT TO RUN - disabled to meet CRAN runtime guidelines

#data(iris)
#hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],bandwidth=0.2,name='setosa')
#hv2 = hypervolume(subset(iris, Species=="versicolor")[,1:4],bandwidth=0.2,name='versicolor')

#hv_all = hypervolume_join(hv1, hv2)

#plot(hv_all)
}
