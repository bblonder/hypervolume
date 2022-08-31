\name{get_occupancy_volume_bootstrap}
\alias{get_occupancy_volume_bootstrap}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the volume from occupancy bootstrap objects
}
\description{
The function \code{get_occupancy_volume_bootstrap()} extract the volume from objects generated with \code{hypervolume_n_resample()}. It can provide raw values or summary statistics for both single groups or pairwise comparisons.
}
\usage{
get_occupancy_volume_bootstrap(path,
                               method = "all",
                               res_type = "raw",
                               relative = FALSE,
                               tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
  A path to a directory containing bootstrapped hypervolumes generated with \code{get_occupancy_bootstrap_volume()}.
}
  \item{method}{
If \code{all} the function computes returns the volume of single bootstrapped hypervolumes. If \code{pairwise} returns the comparison among pairwise combination of hypervolumes.
 }
  \item{res_type}{
  It can be \code{raw} or \code{pairwise}. See details.
}
  \item{relative}{
  If \code{TRUE} it computes relative instead of absolute volumes.
}
  \item{tol}{
  Set the tolerance for reconstructing whole volume.
 }
}
\details{
  The function \code{get_occupancy_volume_bootstrap()} returns the volume of bootstrapped hypervolumes if \code{res_type} is equal to \code{raw} and \code{method} to \code{all}. When \code{res_type} is equal to \code{summary} and method to \code{all} this function returns the mean volume as well as standard deviation, median, minimum, maximum and two quantiles of each group. These summary statistics are calculated for among groups pairwise comparisons when \code{res_type} is equal to \code{summary} and \code{method} to \code{pairwise}.
}
\value{
A \code{data.frame} with bootrapped volumes or summary statistics for single groups or pairwise comparisons.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}, \cr \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_intersection}},  \cr
\code{\link{occupancy_to_union}} \code{\link{occupancy_filter}}
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

# bootstrap the hypervolumes
hv_list_boot = hypervolume_n_resample(name = "example", hv_list)

# calculate occupancy on bootstrapped hypervolumes
hv_occupancy_boot_sex = hypervolume_n_occupancy_bootstrap(path = hv_list_boot,
                                    name = "example_occ",
                                    classification = rep(c("female", "male"), 3))

# get the volume of the bootstrapped hypervolumes
get_occupancy_volume_bootstrap(hv_occupancy_boot_sex)
}
}