\name{hypervolume_variable_importance}
\alias{hypervolume_variable_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume variable importance
}
\description{
Assesses the contribution of each variable to the total hypervolume as a rough metric of variable importance.
}
\usage{
hypervolume_variable_importance(hv, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
A hypervolume for which the importance of each variable should be calculated.

}
  \item{verbose}{
If \code{TRUE}, prints diagnostic progress messages.
}
}
\details{
The algorithm proceeds by comparing the n-dimensional input hypervolume's volume to all possible n-1 dimensional hypervolumes where each variable of interest has been deleted. The importance score reported is the ratio of the n-dimensional hypervolume relative to each of the n-1 dimensional hypervolumes. Larger values indicate that a variable makes a proportionally higher contribution to the overall volume.

The algorithm can only be used on Hypervolumes that have a \code{Data} and \code{Method} value, because the variable deletion process is not well defined for objects that are not associated with a particular set of observations and construction method.

}
\value{
A named vector with importance scores for each axis. Note that these scores are not dimensionless but rather have units corresponding to the original units of each variable.
}
\examples{
# low parameter values for speed
data(iris)
hv = hypervolume_gaussian(subset(iris, Species=="versicolor")[,1:2],samples.per.point=10)
varimp = hypervolume_variable_importance(hv,verbose=FALSE)
barplot(varimp)
}