\name{hypervolume_n_occupancy}
\alias{hypervolume_n_occupancy}
\alias{hypervolume_n_occupancy_bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operations for groups of hypervolumes
}
\description{
Computes the occupancy of hyperspace by one or more groups of hypervolumes.
}
\usage{
hypervolume_n_occupancy(hv_list,
                        classification = NULL,
                        method = "subsample",
                        FUN = mean,
                        num.points.max = NULL,
                        verbose = TRUE,
                        distance.factor = 1,
                        check.hyperplane = FALSE,
                        box_density = 5000,
                        thin = FALSE,
                        quant.thin = 0.5,
                        seed = NULL,
                        print_log = FALSE)
                        
hypervolume_n_occupancy_bootstrap(path,
                                  name = NULL,
                                  classification = NULL,
                                  method = "subsample",
                                  FUN = mean,
                                  num.points.max = NULL,
                                  verbose = TRUE,
                                  distance.factor = 1,
                                  check.hyperplane = FALSE,
                                  box_density = 5000,
                                  thin = FALSE,
                                  quant.thin = 0.5,
                                  seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
An \code{HypervolumeList}.
}
  \item{classification}{
A vector assigning each hypervolume in the \code{HypervolumeList} to a group.
}
  \item{method}{
Can be \code{subsample} or \code{box}. See details.
}
  \item{FUN}{
A function to aggregate points within each group. Default to \code{mean}.
  }
  \item{num.points.max}{
Maximum number of random points to use for set operations. If \code{NULL} defaults to 10^(3+sqrt(n)) where n is the dimensionality of the input hypervolumes. Note that this default parameter value has been increased by a factor of 10 since the 1.2 release of this package.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
  \item{distance.factor}{
Numeric value; multiplicative factor applied to the critical distance for all inclusion tests (see below). Recommended to not change this parameter.
}
  \item{check.hyperplane}{
Check if data is hyperplanar.
}
  \item{box_density}{
Density of random points to fill the hyperbox when method is equal to \code{box}.
}
  \item{thin}{
Take a subsample of random points to get a more uniform distribution of random points. Intended to be used with \code{method = "subsample"}, but can be used with \code{method = "box"} too. Can be slow, especially in high dimensions. See details.
}
  \item{quant.thin}{
Set quantile for using when \code{thin = TRUE}. See details.
}
  \item{seed}{
Set seed for random number generation. Useful for having reproducible results and with the use of \code{\link[=find_optimal_occupancy_thin]{find_optimal_occupancy_thin()}}
}
  \item{print_log}{
Save a log file with the volume of each input hypervolume, recomputed volume and the ratio between the original and recomputed hypervolumes. It works for \code{hypervolume_n_occupancy()} only.
}
\item{path}{
A path to a directory of bootstrapped hypervolumes obtained with \cr
\code{hypervolume_n_resample()}.
}
\item{name}{
File name; The function writes hypervolumes to file in "./Objects/<name>".
}

}
\details{
Uses the inclusion test approach to count how many hypervolumes include each random point. Counts range from 0 (no hypervolumes contain a given random point), to the number of hypervolumes in a group (all the hypervolumes contain a given random point). A function \code{FUN}, usually \code{mean} or \code{sum}, is then applied. A hypervolume is then returned for each group and the occupancy stored in \code{ValueAtRandomPoints}. IMPORTANT: random points with 
\code{ValueAtRandomPoints} equal to 0 are not removed to ease downstream calculation. \cr
When \code{method = "subsample"} the computation is performed on a random sample from input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume and the point density calculated using the volumes of each input hypervolume divided by \code{num.points.max}. \cr
Because this algorithm is based on distances calculated between the distributions of random points, the critical distance (point density ^ (-1/n)) can be scaled by a user-specified factor to provide more or less liberal estimates (\code{distance_factor} greater than or less than 1). \cr
Two methods can be used for calculating the occupancy. The method \code{subsample} is based on a random sample of points from input hypervolumes. Each point is selected with a probability set to the inverse of the number of neighbour points calculated according to the critical distance. This method performs accurately when input hypervolumes have a low degree of overlap. The method \code{box} create a bounding box around the union of input hypervolumes. The bounding box is filled with points following a uniform distribution and with a density set with the argument \code{box_density}. A greater density usually provides more accurate results. The method \code{box} performs better than the method \code{subsample} in low dimensions, while in higher dimensions the method \code{box} become computationally inefficient as nearly all of the hyperbox sampling space will end up being empty and most of the points will be rejected. \cr
When \code{verbose = TRUE} the volume of each input hypervolume will be printed to screen togheter with the recomputed volume and the ratio between the original and recomputed hypervolumes. Mean absolute error (MAE) and root mean square error (RMSE) are also provided as overall measures of the goodness of fit. A log file will be saved in the working directory with the information about the volume of input hypervolumes, the recomputed volume and the ratio between the original and recomputed hypervolumes. \cr
When \code{thin = TRUE} an algorithm is applied to try to make the distribution of random points more uniform. Moderate departures from uniform distribution can in fact result from applying \code{hypervolume_n_occupancy()} on hypervolumes with a high overlap degree. At first, the algorithm in \code{thin} calculates the minimum distance from the neighboor points within the critical distance for each random point. A quantile (set with \code{quant.thin}) of these distances is taken and set as the threshold distance. Random points are then subset so that the distance of a point to another is greater than the threshold distance. \cr
The function \code{hypervolume_n_occupancy_bootstrap()} takes a path of bootstrapped hypervolumes generated with \code{hypervolume_n_resample()} as input. It creates a directory called Objects in the current working directory if a directory of that name doesn't already exist where storing occupancy objects. The function \code{hypervolume_n_occupancy_bootstrap()} returns the absolute path to the directory with bootstrapped hypervolumes. It automatically saves a log file with the volume of each input hypervolume, the recomputed volume and the ratio between the original and recomputed hypervolumes. The log file is used with \code{occupancy_bootstrap_gof()}.
}

\seealso{
\code{\link{find_optimal_occupancy_thin}}, \code{\link{occupancy_bootstrap_gof}}
}


\value{
\code{hypervolume_n_occupancy()} returns a \code{Hypervolume} or \code{HypervolumeList} whose number of hypervolumes equals the number of groups in \code{classification}. \code{hypervolume_n_occupancy_bootstrap()} returns a string containing an absolute path equivalent to ./Objects/<name>.
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
  hypervolume_gaussian(x[, c("bill_length_mm", "flipper_length_mm")],
                       samples.per.point=100, name = y), 
  x = penguins_no_na_split, 
  y = names(penguins_no_na_split))

hv_list <- hypervolume_join(hv_list)

# calculate occupancy without groups
hv_occupancy <- hypervolume_n_occupancy(hv_list)
plot(hv_occupancy, cex.random = 1)

# calculate occupancy with groups
hv_occupancy_list_sex <- hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), each = 3))

plot(hv_occupancy_list_sex, cex.random = 1, show.density = FALSE)


### hypervolume_n_occupancy_bootstrap  ###

# bootstrap the hypervolumes
hv_list_boot = hypervolume_n_resample(name = "example", hv_list)

# calculate occupancy on bootstrapped hypervolumes
hv_occupancy_boot_sex = hypervolume_n_occupancy_bootstrap(path = hv_list_boot,
                                    name = "example_occ",
                                    classification = rep(c("female", "male"), 3))

}
}