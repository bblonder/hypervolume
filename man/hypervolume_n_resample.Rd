\name{hypervolume_n_resample}
\alias{hypervolume_n_resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bootstrap n hypervolumes
}
\description{
The function \code{hypervolume_n_resample()} generates n hypervolumes using data bootstrapped from original data of the input hypervolumes.
}
\usage{
hypervolume_n_resample(name,
                       hv_list,
                       n = 10,
                       points_per_resample = 'sample_size',
                       cores = 1,
                       verbose = TRUE,
                       seed = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
File name; The function writes hypervolumes to file in "./Objects/<name>""
}
  \item{hv_list}{
A \code{Hypervolume} or \code{HypervolumeList} object.
}
  \item{n}{
Number of resamples to take.
Used for every method.
}
  \item{points_per_resample}{
Number of points in each resample. If the input is \code{sample_size}, then the same number of points as the original sample is used.
}
  \item{cores}{
Number of logical cores to use while generating bootstraped hypervolumes. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in \code{cores}.
}
  \item{verbose}{
Logical value; If function is being run sequentially, outputs progress bar in console.
}
  \item{seed}{
Set seed for random number generation.
}
}
\value{
Returns a string containing an absolute path equivalent to ./Objects/<name>
}

\details{
\code{hypervolume_n_resample()} creates a directory called Objects in the current working directory if a directory of that name doesn't already exist. A directory is then created for each hypervolume in \code{hv_list}. Returns an absolute path to directory with resampled hypervolumes. \cr
It is possible to access the hypervolumes by using readRDS to read the hypervolume objects one by one. \cr
The resampled hypervolumes are generated using the same parameters used to generate the input hypervolume. The only exception is that the bandwidth is re-estimated if \code{method = "gaussian"} or \code{method = "box"}. See \code{\link{copy_param_hypervolume}} for more details.
}

\seealso{
\code{\link{hypervolume_n_occupancy_bootstrap}}
}

\examples{
\dontrun{

library(palmerpenguins)
data(penguins)
bill_data = na.omit(penguins[,3:4])
hv = hypervolume(bill_data)

# Example 1: get 50 resampled hypervolumes for each input hypervolume
# Use detectCores to see how many cores are availible in current environment
# Set cores = 1 to run sequentially (default)
# bootstrap the hypervolumes
hv_list_boot = hypervolume_n_resample(name = "example", hv_list, n = 50)

}
}