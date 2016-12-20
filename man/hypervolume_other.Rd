\name{hypervolume_other}
\alias{hypervolume_other}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimates hypervolume for any model with binary or continuous response.
}
\description{
Uses random points within a n-dimensional box to evaluate model predictions, then thresholds predictions to yield a semi-probabilistic geometric shape. It is not recommended to use this function in high dimensionalities because of the high computational costs.
}
\usage{
hypervolume_other(model, range, data = NULL, 
  data_pad_fraction = 1, name = NULL, verbose = TRUE, 
  output.density = 10^ncol(data), threshold_initial = 0, 
  normalize = TRUE, chunksize = 10000, 
  check_memory = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A fitted model on which \code{predict(model, type="response", ...)} can be called.
}
  \item{range}{
A 2 x n element matrix whose columns correspond to predictor variables and whose rows correspond to minimum and maximum range box values over which the function should be evaluated. Could be created as e.g. \code{diff(apply(data, 2, range))}
}
  \item{data}{
If not \code{NULL}, used to determine a minimum range box enclosing the data. These data are used no further in the analysis.
}
  \item{data_pad_fraction}{
A \code{data} is not \code{NULL}, a fraction by which the range box should be extended along each axis.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
  \item{output.density}{
Point density of output hypervolume.
}
  \item{threshold_initial}{
A minimum value at which random points are retained in the hypervolume. Further thresholding beyond this initial value is possible with \code{\link{hypervolume_quantile_threshold}}.
}
  \item{normalize}{
If \code{TRUE}, predicted response values are normalized to probability density values, so that the integral of the response over the whole space is equal to one (consistent with other methods). If \code{FALSE}, response values are unmodified.
}
  \item{chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
  \item{check_memory}{
If \code{TRUE}, estimates the memory usage required to create the n-dimensional random points, then exits. If \code{FALSE}, continues algorithm. It is useful for preventing crashes to check the estimated memory usage on large or high dimensional datasets before running the full algorithm.
}
  \item{\dots}{
Other arguments to be passed to \code{predict}.
}
}
\value{
A \code{Hypervolume} class object corresponding to model predictions above the specified threshold.
}

\seealso{
\code{\link{hypervolume_quantile_threshold}}, \code{\link{expectation_box}}, 
}
\examples{
require(randomForest)
data(iris)
iris$Species <- iris$Species=="versicolor"
m_rf <- randomForest::randomForest(Species~.,data=iris)
hv_rf <- hypervolume_other(model=m_rf,
  data=iris[,1:4],data_pad_fraction=1,output.density=10,
  check_memory=FALSE,threshold_initial = 0.5)
hv_rf_threshold_series <- hypervolume_quantile_threshold(hv_rf,num.thresholds=10)
plot(hv_rf_threshold_series$HypervolumesThresholded,
  showcentroid=FALSE,npmax_random=1000)

}
