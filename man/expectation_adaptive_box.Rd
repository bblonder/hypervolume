\name{expectation_adaptive_box}
\alias{expectation_adaptive_box}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates adaptive hyperbox expectation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
expectation_adaptive_box(data, bin.widths, num.shifts = 1, density = 10^ncol(data))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{bin.widths}{
%%     ~~Describe \code{bin.widths} here~~
}
  \item{num.shifts}{
%%     ~~Describe \code{num.shifts} here~~
}
  \item{density}{
%%     ~~Describe \code{density} here~~
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
function (data, bin.widths, num.shifts = 1, density = 10^ncol(data)) 
{
    binneddata <- bindata_padded(data, bin.widths = bin.widths, 
        num.shifts = num.shifts)
    data_cut <- unique(binneddata$data_cut)
    cuts = binneddata$cuts
    nbins.filled <- nrow(data_cut)
    nbins.possible <- nrow(cuts)^ncol(cuts)
    bin_fraction = nbins.filled/nbins.possible
    cat(sprintf("Exploring \%d / \%.0f bins\n", nbins.filled, nbins.possible))
    rp_all <- vector(mode = "list", length = nrow(data_cut))
    volume_all <- rep(NA, nrow(data_cut))
    cat(sprintf("Random sampling within \%d hyperboxes", nrow(data_cut)))
    for (i in 1:nrow(data_cut)) {
        cat(".")
        lower <- sapply(1:length(cuts), function(dim) {
            cuts[[dim]][data_cut[i, dim]]
        })
        higher <- sapply(1:length(cuts), function(dim) {
            cuts[[dim]][data_cut[i, dim] + 1]
        })
        volume <- prod(higher - lower)
        np <- ceiling(density * volume)
        randpoints <- replicate(length(cuts), runif(n = np))
        for (column_index in 1:ncol(randpoints)) {
            randpoints[, column_index] <- randpoints[, column_index] * 
                (higher[column_index] - lower[column_index]) + 
                lower[column_index]
        }
        rp_all[[i]] <- randpoints
        volume_all[i] <- volume
    }
    cat("done\n")
    cat(sprintf("Binding rows..."))
    rp_all <- do.call("rbind", rp_all)
    cat(sprintf("done\n"))
    volume_all <- sum(volume_all)
    dimnames(rp_all) <- list(NULL, dimnames(data_cut)[[2]])
    result <- new("Hypervolume", Data = as.matrix(data), RandomUniformPointsThresholded = rp_all, 
        PointDensity = density, Volume = volume_all, Dimensionality = ncol(data), 
        ProbabilityDensityAtRandomUniformPoints = rep(1/nrow(rp_all), 
            nrow(rp_all)), Name = sprintf("Expectation - adaptive box (\%d/\%d bins)", 
            nbins.filled, nbins.possible))
    return(result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
