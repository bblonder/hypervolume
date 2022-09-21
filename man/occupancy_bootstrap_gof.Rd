\name{occupancy_bootstrap_gof}
\alias{occupancy_bootstrap_gof}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Goodness of fit metrics for bootstrapped occupancy objects
}
\description{
The \code{occupancy_bootstrap_gof()} function calculates goodness of fit metrics for objects generated with \code{hypervolume_n_occupancy_bootstrap()}.
}
\usage{
occupancy_bootstrap_gof(path, FUN)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
  A path to a directory of bootstrapped hypervolumes generated with \code{hypervolume_n_occupancy_bootstrap()}.
}
  \item{FUN}{
Function to calculate the goodness of fit. It can be \code{mae} for the mean absolute error, \code{rmse} for the root mean square error or a function provided by the user.  
 }
}
\details{
Goodness of fit metrics are calculated on the difference between input and recomputed volumes for each bootstrapped element (set with \code{n} in \code{hypervolume_n_resample()}). See \code{\link[=hypervolume_n_occupancy]{hypervolume_n_occupancy()}} for details on the meaning of input and recomputed hypervolumes. 
}

\value{
A one row \code{data.frame} reporting mean, standard deviation, minimum, maximum, median, 2.5\%, 25\%, 75\% ans 97.5\% quantiles. 
}

\examples{
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))

# split the dataset on species and sex
penguins_no_na_split = split(penguins_no_na, 
paste(penguins_no_na$species, penguins_no_na$sex, sep = "_"))

# calculate the hypervolume for each element of the splitted dataset
hv_list = mapply(function(x, y) 
  hypervolume_gaussian(x[, c("bill_length_mm","bill_depth_mm","flipper_length_mm")],
                       samples.per.point=100, name = y), 
                       x = penguins_no_na_split, 
                       y = names(penguins_no_na_split))


# transform the list into an HypervolumeList
hv_list = hypervolume_join(hv_list)

# bootstrap hypervolumes based on sex
hv_resample = hypervolume_n_resample(hv_list, name = "boot_example")

# calculate occupancy for each bootstrap
hv_occupancy_bootstrap = hypervolume_n_occupancy_bootstrap(hv_resample, 
                                    name = "occupancy_example", 
                                    classification = rep(c("female", "male"), 3))
                
                
occupancy_bootstrap_gof(hv_occupancy_bootstrap, FUN = "rmse")

}
}