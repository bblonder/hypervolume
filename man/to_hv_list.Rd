\name{to_hv_list}
\alias{to_hv_list}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read hypervolumes from directory
}
\description{
Takes a path to a directory containing only rds files and reads them into a HypervolumeList object.
}
\usage{
to_hv_list(path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
absolute or relative path to directory containing rds files
}
}
\details{
Use to_hv_list on the output from \code{hypervolume_resample} when method = "bootstrap" to read bootstrapped hypervolumes into memory.
}
\value{
\code{HypervolumeList} object
}
\examples{
\dontrun{
library(palmerpenguins)
data(penguins)
bill_data = na.omit(penguins[,3:4])
hv = hypervolume(bill_data)

# Use detectCores to see how many cores are availible in current environment
path = hypervolume_resample("example_bootstrap", hv, method = "bootstrap", n = 50, cores = 12)
hvs = to_hv_list(path)
}
}
