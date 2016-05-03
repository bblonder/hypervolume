\name{hypervolume_threshold}
\alias{hypervolume_threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Threshold hypervolume at single or multiple levels
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hypervolume_threshold(hv, thresholds = NULL, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
%%     ~~Describe \code{hv} here~~
}
  \item{thresholds}{
%%     ~~Describe \code{thresholds} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (hv, thresholds = NULL, plot = TRUE) 
{
    if (is.null(thresholds)) {
        thresholds <- seq(min(hv@ProbabilityDensityAtRandomUniformPoints), 
            max(hv@ProbabilityDensityAtRandomUniformPoints), 
            length.out = 100)
    }
    result <- vector(mode = "list", length = length(thresholds))
    for (i in 1:length(thresholds)) {
        hv_new <- hv
        hv_new@RandomUniformPointsThresholded <- hv@RandomUniformPointsThresholded[hv@ProbabilityDensityAtRandomUniformPoints > 
            thresholds[i], , drop = F]
        hv_new@Volume <- length(which(hv@ProbabilityDensityAtRandomUniformPoints > 
            thresholds[i]))/hv@PointDensity
        result[[i]] <- hv_new
    }
    names(result) <- thresholds
    volumes <- sapply(result, function(x) {
        x@Volume
    })
    if (plot == TRUE) {
        plot(volumes ~ thresholds, type = "l", xlab = "Threshold value", 
            ylab = "Hypervolume")
        axis(side = 3, at = seq(thresholds[1], thresholds[length(thresholds)], 
            length.out = 10), labels = seq(1, length(thresholds), 
            length.out = 10))
        mtext(side = 3, line = 3, "Threshold index")
        abline(v = estimate_threshold_gaussian(sd.count = 1, 
            dim = hv@Dimensionality), col = "red")
        abline(v = estimate_threshold_gaussian(sd.count = 2, 
            dim = hv@Dimensionality), col = "blue")
        legend("topright", c("Probability at 1 s.d. distance", 
            "Probability at 2 s.d. distance"), col = c("red", 
            "blue"), cex = 0.5, lty = 1)
    }
    return(list(HypervolumesThresholded = result, Volumes = data.frame(threshold = thresholds, 
        volume = volumes)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
