\name{sample_model_rejection}
\alias{sample_model_rejection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rejection sampling of model
}
\description{
Uses rejection sampling to obtain predicted values of a model object at uniformly random points within a range box
}
\usage{
sample_model_rejection(model, range, N.samples, 
    chunksize = 1000, verbose = TRUE, min.value = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
Any model object which can be used within a \code{predict(model, newdata, ...) call.}
}
  \item{range}{
A 2 x n matrix, where n is the number of dimensions of the data, and the first row corresponds to a lower limit and the second row to an upper limit. Each column is thus the low and high limits of the range box along each axis. Can be generated via \code{\link{padded_range}}.
}
  \item{N.samples}{
Number of samples to draw from the range box.
}
  \item{chunksize}{
Number of samples to process in each \code{predict} call. Changing this value may affect the speed of function return but not the returned values.
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic output.
}
  \item{min.value}{
If \code{TRUE}, discards sampled values below this threshold. Effectively used to set hypervolume boundaries.
}
  \item{\dots}{
Other arguments to be passed to \code{predict}, e.g. \code{type='response'}.
}
}

\value{
To be written
}
