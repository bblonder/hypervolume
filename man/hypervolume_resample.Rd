\name{hypervolume_resample}
\alias{hypervolume_resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume resampling methods
}
\description{
\code{hypervolume_resample} generates new hyperolumes based on the method input. Outputs written to file.

  - \code{"bootstrap"}:        Generates n hypervolumes using data bootstrapped from original data
  
  - \code{"bootstrap seq"}:    Generates n hypervolumes for each sample size in sequence specified by user
  
  - \code{"biased bootstrap"}: Bootstraps input hypervolume with biases applied through multivariate normal weights or user specified weights
}
\usage{
hypervolume_resample(name, 
                      hv, 
                      method, 
                      n = 10, 
                      points_per_resample = "sample_size", 
                      seq = 3:nrow(hv@Data), 
                      k = 5,
                      cores = 1,
                      verbose = TRUE, 
                      mu = NULL, 
                      sigma = NULL, 
                      cols_to_bias = 1:ncol(hv@Data), 
                      weight_func = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
File name; The function writes hypervolumes to file in ./Objects/<name>
}
  \item{hv}{
A hypervolume object
}
  \item{method}{
String input; options are \code{"bootstrap"}, \code{"bootstrap seq"}, and \code{"biased bootstrap"}.
}
  \item{n}{
Number of resamples to take.
Used for every method.
}
  \item{points_per_resample}{
Number of points in each resample. If the input is "sample_size", then the same number of points as the original sample is used.
Used for \code{method = "bootstrap"} and \code{method = "biased bootstrap"}.
}
  \item{seq}{
Sequence of sample sizes. If \code{method = "bootstrap seq"} then the function generates n bootstrapped hypervolumes for each sample size in seq.
Used for \code{method = "bootstrap seq"}.
}
  \item{k}{
Number of splits. Used only for \code{method = "k_split"}.
}
  \item{cores}{
Number of logical cores to use while generating bootstraped hypervolumes. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in cores.
}
  \item{verbose}{
Logical value; If function is being run sequentially, outputs progress bar in console.
}
  \item{mu}{
Array of values specifying the mean of multivariate normal weights.
Used for \code{method = "biased bootstrap"}.
}
  \item{sigma}{
Array of values specifying the variance in each dimension. (Lower variance corresponds to stronger bias)
Used for \code{method = "biased bootstrap"}.
}
  \item{cols_to_bias}{
Array of column indices; must be same length as mu and sigma.  
Used for \code{method = "biased bootstrap"}.
}
  \item{weight_func}{
Custom weight function that takes in a matrix of values and returns desired weights for each row
Used for \code{method = "biased bootstrap"}.
}
}
\details{
\code{hypervolume_resample} creates a directory called Objects in the current working directory if a directory of that name doesn't already exist. Returns an absolute path to directory with resampled hypervolumes. rds files are stored in different file structures depending on which method is called.

Use \code{to_hv_list} to extract every hypervolume object in a directory into a HypervolumeList object.
It is also possible to access the hypervolumes by using readRDS to read the hypervolume objects in one by one.

The resampled hypervolumes are generated using the same parameters used to generate the input hypervolume. The only exception is that the bandwidth is re-estimated if \code{method = "gaussian"} or \code{method = "box"}. See \code{\link{copy_param_hypervolume}} for more details.
}
\value{
returns a string containing an absolute path equivalent to ./Objects/<name>
}

\seealso{
\code{\link{to_hv_list}}, \code{\link{hypervolume_overlap_test}}, \code{\link{hypervolume_funnel}}, \code{\link{hypervolume_overlap_confidence}}
}
\examples{
\dontrun{
library(palmerpenguins)
data(penguins)
bill_data = na.omit(penguins[,3:4])
hv = hypervolume(bill_data)

# Example 1: Get 50 resampled hypervolumes
# Use detectCores to see how many cores are availible in current environment
# Set cores = 1 to run sequentially (default)
path = hypervolume_resample("example_bootstrap", 
                              hv, 
                              method = "bootstrap", 
                              n = 50, 
                              cores = 12)
hvs = to_hv_list(path)

# Example 2: Get resample with applied bias
# Get maximum bill length
max_bill = max(bill_data$bill_length_mm)
# Make data with larger bill length slightly more likley to be resampled
biased_path = hypervolume_resample("biased test", 
                                    hv, 
                                    method = "biased bootstrap", 
                                    n = 50, 
                                    cores = 12, 
                                    mu = max_bill, 
                                    sigma = 90, 
                                    cols_to_bias = "bill_length_mm")
hvs_biased = to_hv_list(biased_path)
}
}

