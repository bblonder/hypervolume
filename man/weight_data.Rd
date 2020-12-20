\name{weight_data}
\alias{weight_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Abundance weighting and prior  of data for hypervolume input
}
\description{
Resamples input data for hypervolume construction, so that some data points can be weighted more strongly than others in kernel density estimation. Also allows a multidimensional normal prior distribution to be placed on each data point to enable simulation of uncertainty or variation within each observed data point. 

Note that this algorithm will change the number of data points and may thus lead to changes in the inferred hypervolume if the selected algorithm (e.g. for bandwidth selection) depends on sample size.

A direct weighting approach (which does not artificially change the sample size, and thus the kernel bandwidth estimate) is available for Gaussian hypervolumes within \code{\link{hypervolume_gaussian}}.
}
\usage{
weight_data(data, weights, jitter.sd = matrix(0, nrow = nrow(data), ncol = ncol(data)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame or matrix of unweighted data. Must only contain numeric values.
}
  \item{weights}{
A vector of weights with the same length as the number of rows in \code{data}. All values must take positive integer values.
}
  \item{jitter.sd}{
A matrix of the same size as \code{data} corresponding to the standard deviation of a normal distribution with mean equal to that of the observed data. If a vector of length equal to 1 or the number of columns of \code{data}, is repeated for all observations.
}
}
\details{
Each data point is jittered a single time. To sample many points from a distribution around each observed data point, multiply all weights by a large number.
}
\value{
A data frame with the rows of \code{data} repeated by \code{weights}, potentially with noise added. The output has the same columns as the input but \code{sum(weights)} total rows.
}
\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

weighted_data <- weight_data(penguins_adelie,
  weights=1+rpois(n=nrow(penguins_adelie),lambda=3))
# color points by alpha to show overlaps
pairs(weighted_data,col=rgb(1,0,0,alpha=0.15)) 

weighted_noisy_data <- weight_data(penguins_adelie,
  weights=1+rpois(n=nrow(penguins_adelie),lambda=3),jitter.sd=0.5)
# color points by alpha to show overlaps
pairs(weighted_noisy_data,col=rgb(1,0,0,alpha=0.15)) 
}
\seealso{\code{\link{hypervolume_gaussian}}}
