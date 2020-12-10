\name{hypervolume_overlap_confidence}
\alias{hypervolume_overlap_confidence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Confidence intervals for overlap statistics
}
\description{
Generates confidence intervals of four different overlap statistics. In order to find the confidence interval for the overlap statistics of two hypervolumes, use \code{hypervolume_resample} twice to generate bootstraps. The function takes in paths to two sets of bootstrapped hypervolumes and gets overlap statistics for each possible pair. Confidence interval is calculated by taking a quantile of generated overlap statistics.
}
\usage{
hypervolume_overlap_confidence(path1, path2, CI = .95, cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path1}{
A path to a directory of bootstrapped hypervolumes
}
\item{path2}{
A path to a directory of bootstrapped hypervolumes
}
\item{CI}{
Desired confidence interval proportion
}
  \item{cores}{
Number of logical cores to use while generating overlap statistics. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in cores.
}
}
\details{
The four overlap statistics are Sorensen, Jaccard, frac_unique_1, frac_unique_2. See \code{\link{hypervolume_overlap_statistics}}

Each hypervolume from path1 is overlapped with each hypervolume from path2 using \code{hypervolume_set}. The four overlap statistics are calculated for each overlap.
}
\value{
  \item{jaccard}{
Confidence interval for jaccard similarity score
}
  \item{sorensen}{
Confidence interval for sorensen similarity score
}
  \item{frac_unique_1}{
Confidence interval for fraction of first hypervolume that is unique
}
  \item{frac_unique_2}{
Confidence interval for fraction of second hypervolume that is unique
}
  \item{distribution}{
a matrix of overlap statistics used to generate the confidence intervals
}
}
\seealso{
\code{\link{hypervolume_resample}}
}
\examples{
\dontrun{
# Let us overlap two hypervolumes generated from multivariate nomral 
# distributions with different means and same covariance matrices.
sample1 = rmvnorm(150, mean = c(0, 0))
sample2 = rmvnorm(150, mean = c(0.5, 0.5))

hv1 = hypervolume(sample1)
hv2 = hypervolume(sample2)

# generates confidence intervals from quantiles of 20*20 overlaps
path1 = hypervolume_resample("mean_0_0", hv1, n = 20)
path2 = hypervolume_resample("mean_0.5_0.5", hv2, n = 20)

result = hypervolume_overlap_confidence(path1, path2)
# confidence index of Sorensen coefficient
print(result["sorensen"])

}
}