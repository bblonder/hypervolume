\name{copy_param_hypervolume}
\alias{copy_param_hypervolume}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate hypervolumes using pre-existing parameters
}
\description{
\code{copy_param_hypervolume} takes in a hypervolume and data. After detecting the method used to generate the input hypervolume, the function returns a new hypervolume generated from the data using the same method and parameters as the input hypervolume. 
}
\usage{
copy_param_hypervolume(hv, data, name = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
hypervolume object
}
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{name}{
String name of hypervolume
}
}
\details{
\code{copy_param_hypervolume} only works if the input hypervolume was generated using \code{method = "box"}, \code{method = "gaussian"}, or \code{method = "svm"}. Calling this function on hypervolumes generated from \code{hypervolume_set} will result in an error. 
Note that \code{kde.bandwidth} is affected by size of the data and will be re-estimated using whichever method was used to generate the original bandwidth if \code{method = "gaussian"} or \code{method = "box"}.
Use \code{hv@Parameters} to see what parameters are copied from the input hypervolume.
}
\value{
hypervolume object
}
\examples{
\dontrun{
library(palmerpenguins)
data("penguins")
bill_data = na.omit(penguins[,3:4])
hv = hypervolume(data = bill_data, 
                  method = "gaussian", 
                  quantile.requested = .9, 
                  quantile.requested.type = "volume")

# Generates a new hypervolume using the same hypervolume and data
hv_copy = copy_param_hypervolume(hv, hv@Data)
# Check to see that the information of the two hypervolumes is the same
print(hv)
print(hv_copy)
}
}

