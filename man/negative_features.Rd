\name{negative_features}
\alias{negative_features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Negative feature detection
}
\description{
Detects the negative features of an observed hypervolume relative to an expectation
}
\usage{
negative_features(hv_obs, hv_exp, set_npoints_max = NULL, set_check_memory = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_obs}{
The observed hypervolume whose negative features are to be detected
}
  \item{hv_exp}{
The expected hypervolume that provides a baseline expectation geometry
}
  \item{set_npoints_max}{
Maximum number of points to be used for set operations comparing \code{hv_obs} to \code{hv_exp}. Defaults to 100*10^n, where n is the dimensionality of the input hypervolumes.
}
  \item{set_check_memory}{
If \code{TRUE}, estimates the memory usage required to perform set operations, then exits. If \code{FALSE}, prints resource usage and continues algorithm. It is useful for preventing crashes to check the estimated memory usage on large or high dimensional datasets before running the full algorithm.
}
}
\details{
This algorithm has a good Type I error rate (rarely detects negative features that do not actually exist). However it can have a high Type II error rate (failure to find negative features when they do exist). To reduce this error rate, make sure to re-run the algorithm with input hypervolumes with higher values of \code{@PointDensity}, or increase \code{set_npoints_max}.

The algorithm performs the set difference between the observed and expected hypervolumes, then removes stray points in this hypervolume by deleting any random point whose distance from any other random point is greater than expected.
}
\value{
A \code{Hypervolume} object containing a uniformly random set of points describing the negative features in \code{hv_obs}. Note that the point density of this object is likely to be much lower than that of the input hypervolumes due to the stochastic geometry algorithms used.
}
\author{
Benjamin Blonder
}


\examples{
# generate annulus data
data_annulus <- data.frame(matrix(data=runif(4000),ncol=2))
names(data_annulus) <- c("x","y")
data_annulus  <- subset(data_annulus, 
sqrt((x-0.5)^2+(y-0.5)^2) > 0.4 & sqrt((x-0.5)^2+(y-0.5)^2) < 0.5)

# MAKE HYPERVOLUME (low reps for fast execution)
hv_annulus <- hypervolume(data_annulus,bandwidth=0.1,name='annulus',reps=500)

# GET CONVEX EXPECTATION
hv_convex <- expectation_convex(hv_annulus, check_memory=FALSE)

# DETECT FEATURES (low npoints for fast execution)
features_annulus <- negative_features(
                      hv_obs=hv_annulus, 
                      hv_exp=hv_convex,
                      set_check_memory=FALSE)


# PLOT RESULTS
plot(hypervolume_join(hv_annulus, features_annulus))
}
