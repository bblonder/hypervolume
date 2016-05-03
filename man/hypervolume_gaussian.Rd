\name{hypervolume_gaussian}
\alias{hypervolume_gaussian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct hypervolume using Gaussian kernel density estimate
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hypervolume_gaussian(data, output.density = 10^(ncol(data)), expectation.num.shifts = 1, expectation.bin.widths = 2 * estimate_bandwidth(data), kde.bandwidth = estimate_bandwidth(data)/2, kde.chunksize = 10000, output.threshold = estimate_threshold_gaussian(sd.count = 1, dim = ncol(data)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{output.density}{
%%     ~~Describe \code{output.density} here~~
}
  \item{expectation.num.shifts}{
%%     ~~Describe \code{expectation.num.shifts} here~~
}
  \item{expectation.bin.widths}{
%%     ~~Describe \code{expectation.bin.widths} here~~
}
  \item{kde.bandwidth}{
%%     ~~Describe \code{kde.bandwidth} here~~
}
  \item{kde.chunksize}{
%%     ~~Describe \code{kde.chunksize} here~~
}
  \item{output.threshold}{
%%     ~~Describe \code{output.threshold} here~~
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
function (data, output.density = 10^(ncol(data)), expectation.num.shifts = 1, 
    expectation.bin.widths = 2 * estimate_bandwidth(data), kde.bandwidth = estimate_bandwidth(data)/2, 
    kde.chunksize = 10000, output.threshold = estimate_threshold_gaussian(sd.count = 1, 
        dim = ncol(data))) 
{
    expectation <- expectation_adaptive_box(data, density = output.density, 
        num.shifts = expectation.num.shifts, bin.widths = expectation.bin.widths)
    cat("running KDE...\n")
    if (ncol(data) == 1) {
        Hmatrix <- diag(1, nrow = ncol(data), ncol = ncol(data)) * 
            (kde.bandwidth)
    }
    else {
        Hmatrix <- diag(1, nrow = ncol(data), ncol = ncol(data)) * 
            (kde.bandwidth)^2
    }
    np <- nrow(expectation@RandomUniformPointsThresholded)
    cat(sprintf("Sampling \%d random points from kernel density estimate, \%d per chunk...\n", 
        np, kde.chunksize))
    num.samples.completed <- 0
    num.chunks <- ceiling(np/kde.chunksize)
    kde.probs <- vector(mode = "list", length = num.chunks)
    for (i in 1:num.chunks) {
        cat(".")
        num.samples.to.take <- min(kde.chunksize, np - num.samples.completed)
        kde.probs.this <- ks::kde(x = data, H = Hmatrix, eval.points = expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed + 
            num.samples.to.take - 1), ], verbose = T)$estimate
        kde.probs[[i]] <- kde.probs.this
        num.samples.completed <- num.samples.completed + num.samples.to.take
    }
    kde.probs <- do.call("c", kde.probs)
    cat("done.\n")
    points_thresholded <- expectation@RandomUniformPointsThresholded[kde.probs > 
        output.threshold, , drop = F]
    probs_thresholded <- kde.probs[kde.probs > output.threshold]
    numpoints_kde <- nrow(points_thresholded)
    vol_kde <- numpoints_kde/expectation@PointDensity
    hv_kde <- new("Hypervolume", Data = as.matrix(data), RandomUniformPointsThresholded = points_thresholded, 
        PointDensity = expectation@PointDensity, Volume = vol_kde, 
        Dimensionality = ncol(data), ProbabilityDensityAtRandomUniformPoints = probs_thresholded, 
        Name = "Gaussian kernel density estimate")
    return(hv_kde)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
