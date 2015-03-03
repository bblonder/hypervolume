\name{hypervolume_importance}
\alias{hypervolume_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume variable importance
}
\description{
Assesses the contribution of each variable to the total hypervolume as a rough metric of variable importance. 
}
\usage{
hypervolume_importance(hv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
A hypervolume for which the importance of each variable should be calculated.
}
}
\details{
The algorithm proceeds by comparing the n-dimensional input hypervolume's volume to all possible n-1 dimensional hypervolumes where each variable of interest has been deleted. The importance score reported is the ratio of the n-dimensional hypervolume relative to each of the n-1 dimensional hypervolumes. Larger values indicate that a variable makes a proportionally higher contribution to the overall volume.

This function is currently experimental. Other definitions of variable importance may also be possible. The algorithm currently cannot be used on set operations hypervolumes because the variable deletion process is not well defined for objects that are not associated with a particular set of observations.
}
\value{
A named vector with importance scores for each axis. Note that these scores are not dimensionless but rather have linear units.
}
\examples{
data(iris)
hv1 = hypervolume(subset(iris, Species=="setosa")[,1:4],reps=1000,bandwidth=0.2)
imp = hypervolume_importance(hv1)

print(imp)
}
