\name{hypervolume_sorensen_overlap}
\alias{hypervolume_sorensen_overlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sorensen similarity index for hypervolume set operations
}
\description{
Calculates the Sorensen index, i.e. for hypervolumes A and B, 2*|A int B| / (|A| + |B|). Note that this definition has been changed since version 1.2 of the package.
}
\usage{
hypervolume_sorensen_overlap(hvlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hvlist}{
A \code{HypervolumeList} object corresponding to the output of a \code{hypervolume_set} command.
}
}

\value{
A number between 0 and 1 indicating how similar two input hypervolumes are.
}

\examples{
data(iris)
hv2 = hypervolume(subset(iris, Species=="virginica")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)
hv3 = hypervolume(subset(iris, Species=="versicolor")[,1:4],reps=1000,bandwidth=0.2,warn=FALSE)

hv_set23 = hypervolume_set(hv2, hv3, check_memory=FALSE)

# no overlap found between setosa and virginica
hypervolume_sorensen_overlap(hv_set23)
}