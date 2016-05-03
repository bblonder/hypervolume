\name{hypervolume_svm}
\alias{hypervolume_svm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via one-class support vector machine with radial kernel
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hypervolume_svm(data, output.density = 10^ncol(data), expectation.num.shifts = 1, expectation.bin.widths = 2 * estimate_bandwidth(data), svm.nu = 0.01, svm.gamma = 0.5, svm.chunksize = 10000)
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
  \item{svm.nu}{
%%     ~~Describe \code{svm.nu} here~~
}
  \item{svm.gamma}{
%%     ~~Describe \code{svm.gamma} here~~
}
  \item{svm.chunksize}{
%%     ~~Describe \code{svm.chunksize} here~~
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
function (data, output.density = 10^ncol(data), expectation.num.shifts = 1, 
    expectation.bin.widths = 2 * estimate_bandwidth(data), svm.nu = 0.01, 
    svm.gamma = 0.5, svm.chunksize = 10000) 
{
    data <- as.matrix(data)
    expectation <- expectation_adaptive_box(data, density = output.density, 
        num.shifts = expectation.num.shifts, bin.widths = expectation.bin.widths)
    cat("running SVM")
    svm.model <- e1071::svm(data, y = NULL, type = "one-classification", 
        nu = svm.nu, gamma = svm.gamma, scale = TRUE, kernel = "radial")
    cat("...done\n")
    np <- nrow(expectation@RandomUniformPointsThresholded)
    cat(sprintf("Sampling \%d random points from support vector machine, \%d per chunk...\n", 
        np, svm.chunksize))
    num.samples.completed <- 0
    num.chunks <- ceiling(np/svm.chunksize)
    svm.probs <- vector(mode = "list", length = num.chunks)
    for (i in 1:num.chunks) {
        cat(".")
        num.samples.to.take <- min(svm.chunksize, np - num.samples.completed)
        inputdata.this <- expectation@RandomUniformPointsThresholded[num.samples.completed:(num.samples.completed + 
            num.samples.to.take - 1), ]
        svm.pred.this <- predict(svm.model, inputdata.this)
        svm.pred.this.pos <- inputdata.this[svm.pred.this == 
            TRUE, ]
        svm.probs[[i]] <- svm.pred.this.pos
        num.samples.completed <- num.samples.completed + num.samples.to.take
    }
    svm.probs <- do.call("rbind", svm.probs)
    dimnames(svm.probs) <- list(NULL, dimnames(data)[[2]])
    cat("done.\n")
    numpoints_svm <- nrow(svm.probs)
    vol_svm <- numpoints_svm/expectation@PointDensity
    hv_svm <- new("Hypervolume", Data = data, RandomUniformPointsThresholded = svm.probs, 
        PointDensity = expectation@PointDensity, Volume = vol_svm, 
        Dimensionality = ncol(svm.probs), ProbabilityDensityAtRandomUniformPoints = rep(1/nrow(svm.probs), 
            nrow(svm.probs)), Name = "Support vector machine")
    return(hv_svm)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
