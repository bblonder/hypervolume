\name{get_occupancy_volume_bootstrap}
\alias{get_occupancy_volume_bootstrap}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the volume from occupancy bootstrap objects
}
\description{
The function \code{get_occupancy_volume_bootstrap()} extract the volume from objects generated with \code{hypervolume_n_occupancy_bootstrap()}. It provides raw values or summary statistics for both single hypervolumes or their pairwise comparisons.
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
  A path to a directory containing bootstrapped occupancy objects generated with \code{hypervolume_n_occupancy_bootstrap()}.
}
  \item{method}{
If \code{all} the function returns the volume of each bootstrapped hypervolume for each bootstrapped occupancy object in \code{path}. If \code{pairwise} returns the volume difference for each pairwise combination of hypervolumes within the bootstrapped occupancy objects found in \code{path}.
 }
  \item{res_type}{
  It can be \code{raw} or \code{pairwise}. See details.
}
  \item{relative}{
  If \code{TRUE} it computes relative instead of absolute volumes.
}
  \item{tol}{
  Set the tolerance for reconstructing whole volume. See details.
 }
}
\details{
The function \code{get_occupancy_volume_bootstrap()} returns the volume for each bootstrapped hypervolume if \code{res_type = "raw"} and \code{method = "all"}. When \code{res_type = "summary"} and \code{method = "all"} this function returns the mean volume as well as the standard deviation, median, minimum, maximum, 2.5\% and 97.5\% quantiles, skewness and kurtosis for each of hypervolume. The same summary statistics are calculated for the difference of volume for each pairwise combination of hypervolumes when \code{res_type = "summary"} and \code{method = "pairwise"}. \cr
When \code{relative = TRUE} relative instead of absolute volumes are returned. The relative volume is calculated as the ratio between the volume of an hypervolume and the volume of the union of all the hypervolumes. The \code{get_occupancy_volume_bootstrap()} function attempts to reconstruct the volume of the union of all the hypervolumes from each bootstrapped hypervolume. At first, the volume of the union of hypervolumes is calculated for each hypervolume of the jth bootstrapped occupancy_object as the the ratio between the total number of random points and the number of random points of the ith hypervolume of the jth bootstrapped occupancy_object, multiplied by the volume of the ith hypervolume of the jth bootstrapped occupancy_object. This step results in a number of reconstructed volumes equal to the number of hypervolumes in the jth bootstrapped occupancy_object. Reconstructed volumes are then compared among each other to ensure the consistency of the reconstruction. To do this, the distance among reconstructed volumes is calculated with the \code{dist()} function of the \code{stats} package. If at least one of the distances is greater than \code{tol} the computation is stopped and some suggestions are returned.
}

\value{
A \code{data.frame} with bootstrapped volumes or summary statistics for single hypervolumes or their pairwise comparisons.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}
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