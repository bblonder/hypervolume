\name{get_occupancy_stats}
\alias{get_occupancy_stats}
\alias{get_occupancy_stats_bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stats from occupancy objects
}
\description{
Functions \code{get_occupancy_stats()} and \code{get_occupancy_stats_bootstrap()} return the results of a function applied to hypervolumes generated with \code{hypervolume_n_occupancy()}, \cr
\code{hypervolume_n_occupancy_bootstrap()} or \code{hypervolume_n_occupancy_test()}.
}
\usage{
get_occupancy_stats(hv, FUN, remove_zeroes = TRUE)

get_occupancy_stats_bootstrap(path,
                              FUN,
                              remove_zeroes = TRUE,
                              method = "pairwise",
                              res_type = "summary",
                              verbose = TRUE,
                              cores = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
  A \code{Hypervolume} or \code{HypervolumeList} object generated with \cr
  \code{hypervolume_n_occupancy()} or \code{hypervolume_n_occupancy_test()}.
}
  \item{FUN}{
  The function to be applied.
}
  \item{remove_zeroes}{
  Remove zeroes before the calculation. See Details for further information.
}
  \item{path}{
  A path to a directory of bootstrapped hypervolumes obtained with \cr
  \code{hypervolume_n_occupancy_bootstrap()}.
}
  \item{method}{
  If \code{all} compares the target hypervolume with all the other hypervolume. If \code{pairwise} compare all the pairwise combinations of hypervolumes.
}
  \item{res_type}{
  If \code{raw} returns the raw values obtained with \code{FUN}. If \code{summary} returns 
}
  \item{verbose}{
  Logical value; print diagnostic output if \code{TRUE}.
}
  \item{cores}{
  Number of logical cores to use while generating permuted hypervolumes. If parallel backend already registered to doParallel, function will use that backend and ignore the argument in cores.
}


}
\value{
Either a \code{vector}, a \code{matrix} or a \code{data.frame} with the results of the applied function. 
}

\details{
These functions remove occupancy values equal to 0 by default. These values are generated
during the occupancy routine when a random point is included 
in some group of hypervolumes but not in others. A tipical usage of
\code{get_occupancy_stats()} or \code{get_occupancy_stats_bootstrap()} should remove zeroes before applying a function (the default).
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}, \cr \code{\link{hypervolume_n_occupancy_test}}
}

\examples{
\dontrun{

##### single occupancy object #####
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

# calculate occupancy based on sex
hv_occupancy_list_sex = hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), 3))


# calculate the mean occupancy value
get_occupancy_stats(hv_occupancy_list_sex, mean)

##### bootstrapped occupancy objects #####

# bootstrap input hypervolumes
hv_boot = hypervolume_n_resample(name = "example", hv_list = hv_list,  n = 9)

# calculate occupancy on bootstrapped hypervolumes
hv_boot_occ = hypervolume_n_occupancy_bootstrap(hv_boot, name = "example_occ",
                                  classification = rep(c("female", "male"), 3))

# calculate summary statistics for pairwise comparisons
get_occupancy_stats_bootstrap(hv_boot_occ, FUN = mean)

}
}