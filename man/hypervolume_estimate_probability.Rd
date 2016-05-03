\name{hypervolume_estimate_probability}
\alias{hypervolume_estimate_probability}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate probability within a hypervolume at an arbitrary point
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hypervolume_estimate_probability(hv, points, reduction_factor = 1, verbose = T, distance_factor = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
%%     ~~Describe \code{hv} here~~
}
  \item{points}{
%%     ~~Describe \code{points} here~~
}
  \item{reduction_factor}{
%%     ~~Describe \code{reduction_factor} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{distance_factor}{
%%     ~~Describe \code{distance_factor} here~~
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
function (hv, points, reduction_factor = 1, verbose = T, distance_factor = 1) 
{
    np = nrow(hv@RandomUniformPointsThresholded)
    dimhv = ncol(hv@RandomUniformPointsThresholded)
    dimp = ncol(points)
    if (dimp != dimhv) {
        stop("Dimensionality of hypervolume and points is not the same.")
    }
    if (reduction_factor <= 0 | reduction_factor > 1) {
        stop("Reduction factor is not in (0,1].")
    }
    numpointstokeep_hv = floor(np * reduction_factor)
    if (reduction_factor < 1) {
        hv_points_ss = hv@RandomUniformPointsThresholded[sample(1:np, 
            size = numpointstokeep_hv), ]
    }
    else {
        hv_points_ss = hv@RandomUniformPointsThresholded
    }
    if (verbose == TRUE) {
        cat(sprintf("Retaining \%d points for \%d inclusion tests.\n", 
            numpointstokeep_hv, nrow(points)))
    }
    point_density = nrow(hv_points_ss)/hv@Volume
    cutoff_dist = point_density^(-1/dimhv) * distance_factor
    points_in_hv_all_list = evalfspherical(hv_points_ss, cutoff_dist, 
        points, verbose = verbose, getid.nearestneighbor = TRUE)
    points_in_hv_all_list_probs <- points_in_hv_all_list[[2]]
    if (verbose == TRUE) {
        cat(sprintf("Estimating probability at \%d points...", 
            length(points_in_hv_all_list_probs)))
    }
    df_probs <- as.data.frame(do.call("rbind", lapply(1:length(points_in_hv_all_list_probs), 
        function(i) {
            cbind(index.point = i, index.prob = points_in_hv_all_list_probs[[i]])
        })))
    df_probs <- df_probs[df_probs$index.prob > -1, ]
    df_probs$prob <- hv@ProbabilityDensityAtRandomUniformPoints[df_probs$index.prob + 
        1]
    probabilities <- tapply(df_probs$prob, df_probs$index.point, 
        mean)
    probabilities_out <- rep(0, length(points_in_hv_all_list_probs))
    print(str(probabilities_out))
    probabilities_out[as.numeric(names(probabilities))] <- probabilities
    print(str(probabilities_out))
    return(probabilities_out)
    probabilities <- sapply(points_in_hv_all_list_probs, function(indices) {
        cat(".")
        indices_startatone <- indices + 1
        mean_probability <- mean(hv@ProbabilityDensityAtRandomUniformPoints[indices_startatone])
        return(mean_probability)
    })
    if (verbose == TRUE) {
        cat(sprintf("done.\n"))
    }
    probabilities[points_in_hv_all_list[[1]] == 0] <- 0
    return(probabilities)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
