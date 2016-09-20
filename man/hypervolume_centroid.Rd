\name{hypervolume_centroid}
\alias{hypervolume_centroid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume centroid
}
\description{
Computes the position of the center of mass (centroid) of a hypervolume
}
\usage{
hypervolume_centroid(hv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
An input hypervolume
}
}
\value{
A vector with length equal to the dimensionality of the input hypervolume.
}
\examples{
data(iris)
hv = hypervolume_box(data=subset(iris, Species=="setosa")[,1:3],name='setosa')
hypervolume_centroid(hv)
}