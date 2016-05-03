\name{hypervolume_project}
\alias{hypervolume_project}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Species distribution modeling - geographic projection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hypervolume_project(hv, rasters, do.probability = TRUE, normalize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
%%     ~~Describe \code{hv} here~~
}
  \item{rasters}{
%%     ~~Describe \code{rasters} here~~
}
  \item{do.probability}{
%%     ~~Describe \code{do.probability} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
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
function (hv, rasters, do.probability = TRUE, normalize = TRUE) 
{
    raster.values = data.frame(raster::getValues(rasters))
    if (do.probability == TRUE) {
        projected.values <- hypervolume_estimate_probability(hv, 
            raster.values)
        print(str(projected.values))
        if (normalize == TRUE) {
            projected.values <- projected.values/sum(projected.values)
        }
    }
    else {
        projected.values = hypervolume_inclusion_test(hv, raster.values)
    }
    raster.out = raster(rasters[[1]])
    values(raster.out) <- projected.values
    return(raster.out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
