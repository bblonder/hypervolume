\name{hypervolume_funnel}
\alias{hypervolume_funnel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolumes at different sample sizes
}
\description{
This function takes in hypervolumes bootstrapped at different sample sizes applies a function to each hypervolume. The output of the function can either be a plot of nonparametric confidence intervals or a table of the mean and quantiles.
}
\usage{
hypervolume_funnel(input_path, 
                    title = NULL, 
                    func = get_volume, 
                    CI = .95, 
                    as_table = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_path}{
output of \code{resample} with method = "bootstrap seq"; path to a sequence of different sample sized bootstraps
}
  \item{title}{
title of output plot, ignore if outputting as table
}
  \item{func}{
a function that takes a single parameter which is a hypervolume and returns a numerical value.
}
  \item{CI}{
Confidence interval is taken by using the the (1-CI)/2 and (1+CI)/2 quantile
}
  \item{as_table}{
If TRUE, returns a table with columns upper quantile, mean, lower quantile
}
}
\details{
This function is used to evaluate the behavior of hypervolumes at different sample sizes and determine bias. Statistics such as volume are affected by sample size especially when the hypervolumes are constructed with method = "gaussian" since the bandwidth estimate is dependent on sample size.
}
\value{
ggplot object, or dataframe object
}
\examples{
\dontrun{
# 3000 data point hypervolume
data(quercus)
hv_quercus = hypervolume(quercus[,c(2,3)])

# the seq argument is equivalent to a length 30 vector {10, 139, ... , 3649, 3779}
# 6hr sequential runtime
quercus_bootstrap_seq <- resample('quercus_bootstrap_seq', 
                                    hv_quercus, 
                                    method = 'bootstrap seq', 
                                    points_per_resample = "sample_size", 
                                    seq = floor(seq(10, 3779, length.out = 30)),
                                    cores = 20)

# Compatible with ggplot syntax when used with as_table = FALSE
hypervolume_funnel(quercus_bootstrap_seq, 
                    title = 'Resampled volumes of Quercus', 
                    func = get_volume) + 
  geom_line(aes(y = get_volume(hv_quercus))) +
  ylab("Volume")
}
}

