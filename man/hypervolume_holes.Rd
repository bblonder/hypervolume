\name{hypervolume_holes}
\alias{hypervolume_holes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hole detection
}
\description{
Detects the holes in an observed hypervolume relative to an expectation
}
\usage{
hypervolume_holes(hv.obs, hv.exp, set.num.points.max = NULL, set.check.memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv.obs}{
The observed hypervolume whose holes are to be detected
}
  \item{hv.exp}{
The expected hypervolume that provides a baseline expectation geometry
}
  \item{set.num.points.max}{
Maximum number of points to be used for set operations comparing \code{hv_obs} to \code{hv_exp}. Defaults to 10^(3+sqrt(n)), where n is the dimensionality of the input hypervolumes.
}
  \item{set.check.memory}{
If \code{TRUE}, estimates the memory usage required to perform set operations, then exits. If \code{FALSE}, prints resource usage and continues algorithm. It is useful for preventing crashes to check the estimated memory usage on large or high dimensional datasets before running the full algorithm.
}
}
\details{
This algorithm has a good Type I error rate (rarely detects holes that do not actually exist). However it can have a high Type II error rate (failure to find holes when they do exist). To reduce this error rate, make sure to re-run the algorithm with input hypervolumes with higher values of \code{@PointDensity}, or increase \code{set.num.points.max}.

The algorithm performs the set difference between the observed and expected hypervolumes, then removes stray points in this hypervolume by deleting any random point whose distance from any other random point is greater than expected.

A 'rule of thumb' is that algorithm has acceptable statistical performance when log_e(m) > n, where m is the number of data points and n is the dimensionality.
}
\value{
A \code{Hypervolume} object containing a uniformly random set of points describing the holes in \code{hv_obs}. Note that the point density of this object is likely to be much lower than that of the input hypervolumes due to the stochastic geometry algorithms used.
}

\examples{
\dontrun{
# generate annulus data
data_annulus <- data.frame(matrix(data=runif(4000),ncol=2))
names(data_annulus) <- c("x","y")
data_annulus  <- subset(data_annulus, 
sqrt((x-0.5)^2+(y-0.5)^2) > 0.4 & sqrt((x-0.5)^2+(y-0.5)^2) < 0.5)

# MAKE HYPERVOLUME (low reps for fast execution)
hv_annulus <- hypervolume_gaussian(data_annulus,
                kde.bandwidth=0.05,name='annulus',samples.per.point=1)

# GET CONVEX EXPECTATION
hv_convex <- expectation_convex(hypervolume_thin(hv_annulus,num.samples=500),
                check.memory=FALSE,use.random=TRUE)

# DETECT HOLES (low npoints for fast execution)
features_annulus <- hypervolume_holes(
                      hv.obs=hv_annulus, 
                      hv.exp=hv_convex,
                      set.check.memory=FALSE)

# CLEAN UP RESULTS
features_segmented <- hypervolume_segment(features_annulus, 
                        check.memory=FALSE,distance.factor=2)
features_segmented_pruned <- hypervolume_prune(features_segmented, 
                                volume.min=0.02)

# PLOT RETAINED HOLE(S)
plot(hypervolume_join(hv_annulus, features_segmented_pruned))
}
}
