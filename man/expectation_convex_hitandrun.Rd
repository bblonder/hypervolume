\name{expectation_convex_hitandrun}
\alias{expectation_convex_hitandrun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random sampling of convex polytope using hit-and-run sampling
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
expectation_convex_hitandrun(data, point_density, hitandrun.chunksize = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{point_density}{
%%     ~~Describe \code{point_density} here~~
}
  \item{hitandrun.chunksize}{
%%     ~~Describe \code{hitandrun.chunksize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, point_density, hitandrun.chunksize = 1e+05) 
{
    data <- as.matrix(data)
    cat("Calculating convex hull...")
    hull <- convhulln(data, options = "FA")
    cat("done.\n")
    cat("Calculating linear constraints...\n")
    constraints <- hullToConstr(data, hull$hull)
    cat("done.\n")
    volume_convexhull <- hull$vol
    np <- ceiling(point_density * volume_convexhull)
    cat(sprintf("Sampling \%d random points via hit-and-run, \%d per chunk...\n", 
        np, hitandrun.chunksize))
    num.samples.completed <- 0
    num.chunks <- ceiling(np/hitandrun.chunksize)
    samples <- vector(mode = "list", length = num.chunks)
    for (i in 1:num.chunks) {
        cat(".")
        num.samples.to.take <- min(hitandrun.chunksize, np - 
            num.samples.completed)
        samples_this <- hitandrun::hitandrun(constraints, n.samples = num.samples.to.take)
        samples[[i]] <- samples_this
        num.samples.completed <- num.samples.completed + num.samples.to.take
    }
    samples <- do.call("rbind", samples)
    dimnames(samples) <- list(NULL, dimnames(data)[[2]])
    cat("done.\n")
    hv_hull <- new("Hypervolume", Data = data, RandomUniformPointsThresholded = samples, 
        PointDensity = point_density, Volume = volume_convexhull, 
        Dimensionality = ncol(samples), ProbabilityDensityAtRandomUniformPoints = rep(1/nrow(samples), 
            nrow(samples)), Name = "Convex hull")
    return(hv_hull)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
