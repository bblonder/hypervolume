\name{get_occupancy_unshared_bootstrap}
\alias{get_occupancy_unshared_bootstrap}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Volume of the unshared fraction of a bootstrapped occupancy object
}
\description{
The \code{get_occupancy_unshared_bootstrap()} function is used to get the volume of the unshared fraction of an object generated with \code{hypervolume_n_occupancy_bootstrap()}. It provides raw values or summary statistics for both individual hypervolumes or their pairwise comparisons.
}
\usage{
get_occupancy_unshared_bootstrap(path,
                                 method = "pairwise",
                                 res_type = "summary",
                                 relative = FALSE,
                                 tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
  A path to a directory of bootstrapped occupancy objects obtained with \cr
  \code{hypervolume_n_occupancy_bootstrap()}.
}
  \item{method}{
If \code{all} compute the volume of the unique fraction of each hypervolume compared to all the hypervolumes for each occupancy object in \code{path}. If \code{pairwise} compute the difference of the volume of the unshared fraction for each pairwise combination of hypervolumes within the bootstrapped occupancy objects found in \code{path}.
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
The function \code{get_occupancy_unshared_bootstrap()} returns the volume of the unshared fraction for each hypervolume in the bootstrapped occupancy object if \code{res_type = "raw"} and \code{method = "all"}. When \code{res_type = "summary"} and \code{method = "all"} this function returns the mean volume as well as the standard deviation, median, minimum, maximum, 2.5\% and 97.5\% quantiles, skewness and kurtosis of the unshared fraction for each hypervolume. The same summary statistics are calculated for the difference of volume of the unshared fraction for each pairwise combination of hypervolumes when \code{res_type = "summary"} and \code{method = "pairwise"}. The unshared fraction is calculated by finding the set of random points that are not shared with other hypervolumes or pairwise combinations of hypervolumes in each bootstrapped occupancy object. More details on how the unshared fraction is computed in \code{\link[=occupancy_to_unshared]{occupancy_to_unshared()}}.\cr
The \code{get_occupancy_unshared_bootstrap()} function attempts to reconstruct the volume of the unshared fraction from each bootstrapped occupancy object. At first, the volume of the union of hypervolumes is calculated for each hypervolume of the jth bootstrapped occupancy object as the ratio between the total number of random points and the number of random points of the ith hypervolume of the jth bootstrapped occupancy object, multiplied by the volume of the ith hypervolume of the jth bootstrapped occupancy object. This step results in a number of reconstructed volumes equal to the number of hypervolumes in the jth bootstrapped occupancy object. Reconstructed volumes are then compared among each other to ensure the consistency of the reconstruction. To do this, the distance among reconstructed volumes is calculated with the \code{dist()} function of the \code{stats} package. If at least one of the distances is greater than \code{tol} the computation is stopped and some suggestions are returned. The volume of the unshared fraction is then calculated as the ratio between the number of random points of the unshared fraction and the total number of random points, multiplied by the volume of the union of hypervolumes. \cr
When \code{relative = TRUE} relative instead of absolute volumes are returned. The relative volume is calculated as the ratio between the volume of the unshared fraction and the volume of the union of all the hypervolumes (or combination of hypervolumes when \code{method = "pairwise"}). The same approach described above is used to reconstruct the volume of the union of hypervolumes.
}

\value{
A \code{data.frame} with bootstrapped volumes or summary statistics of the unshared fraction.
}


\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}
\code{\link{occupancy_to_unshared}}
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

# get the unshared fraction
get_occupancy_unshared_bootstrap(hv_occupancy_boot_sex)

}
}