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
One or more objects of class \code{Hypervolume} or \code{HypervolumeList}, or a \code{list()} of \code{Hypervolume} objects.
}
}

\value{
A \code{HypervolumeList} containing all hypervolumes in all arguments.
}

\examples{
# data(iris)
# data_split = split(iris[,1:3],iris$Species)
# hvs_split = lapply(data_split, hypervolume); 
# hvs_joined = hypervolume_join(hvs_split)
}
