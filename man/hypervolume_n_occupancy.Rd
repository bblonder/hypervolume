\name{hypervolume_n_occupancy}
\alias{hypervolume_n_occupancy}
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
  box_density = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
An HypervolumeList.
}
  \item{classification}{
A vector assigning each Hypervolume in the HypervolumeList to a group.
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
Logical value; print diagnostic output if true.
}
  \item{distance.factor}{
Numeric value; multiplicative factor applied to the critical distance for all inclusion tests (see below). Recommended to not change this parameter.
}
  \item{check.hyperplane}{
Check if data is hyperplanar.
}
  \item{box_density}{
Density of random point to fill the hyperbox when method is equal to \code{box}.
}
}
\details{
Uses the inclusion test approach to count how many hypervolumes in each group includes random points. Counts range from 0 (no hypervolume contains a given random point), to the number of hypervolumes in a group (all the hypervolumes contains a given random point). A function \code{FUN}, usually \code{mean} or \code{sum}, is then applied. An hypervolume is then returned for each group and the occupancy stored in \code{@ValueAtRandomPoints}. IMPORTANT: random points with 
\code{@ValueAtRandomPoints} equal to 0 are not removed to ease downstream calculation.  

The computation is actually performed on a random sample from input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume, and the point density calculated using the volumes of each input hypervolume divided by \code{num.points.max}.

Because this algorithm is based on distances calculated between the distributions of random points, the critical distance (point density ^ (-1/n)) can be scaled by a user-specified factor to provide more or less liberal estimates (\code{distance_factor} greater than or less than 1).

Two methods can be used for calculating the occupancy. The method \code{subsample} is based on a random sample of points from input hypervolumes. Each point is selected with a probability set to the inverse of the number of neighbour points calculated according to the critical distance. This method performs accurately when input hypervolumes have a low degree of overlap. The method \code{box} create a bounding box around the union of input hypervolumes. The bounding box is filled with points following a uniform distribution and with a density set with the argument \code{box_density}. A greater density provides more accurate results. The method \code{box_density} performs better than the method \code{subsample} in low dimensions, while in higher dimensions \code{box_density} become computationally inefficient as nearly all of the hyperbox sampling space will end up being empty and most of the points will be rejected.

When \code{verbose = TRUE} the volume of each input hypervolume will be printed togheter with the recomputed volume and the ratio between the original and recomputed hypervolumes. Mean absolute error (MAE) and root mean square error (RMSE) will be also provided as overall measures of the goodness of fit.   

}




\value{
\code{hypervolume_n_occupancy} returns an HypervolumeList whose number of elements equals the number of groups in \code{classification}.
}

\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))


penguins_no_na_split = split(penguins_no_na, 
                        paste(penguins_no_na$species, penguins_no_na$sex, sep = "_"))


hv_list = lapply(penguins_no_na_split, function(x) 
              hypervolume_gaussian(x[, c("bill_length_mm","bill_depth_mm","flipper_length_mm")],
              samples.per.point=75))

names(hv_list) <- names(penguins_no_na_split)
hv_list <- hypervolume_join(hv_list)

hv_occupancy <- hypervolume_n_occupancy(hv_list)
plot(hv_occupancy, cex.random = 1)

hv_occupancy_list_sex <- hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), each = 3))

plot(hv_occupancy_list_sex, cex.random = 1, show.density = FALSE)

}
