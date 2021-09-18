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
  FUN = mean,
  num.points.max = NULL,
  verbose = TRUE,
  distance.factor = 1,
  check.hyperplane = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
An HypervolumeList.
}
  \item{classification}{
A vector assigning each Hypervolume in the HypervolumeList to a group.
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
}
\details{
Uses the inclusion test approach to count how many hypervolumes in each group includes random points. Counts range from 0 (no hypervolume contains a given random point), to the number of hypervolumes in a group (all the hypervolumes contains a given random point). A function \code{FUN}, usually \code{mean} or \code{sum}, is then applied. An hypervolume is then returned for each group and the occupancy stored in \code{@ValueAtRandomPoints}. IMPORTANT: random points with \code{@ValueAtRandomPoints} equal to 0 are not removed to ease downstream calculation.  

The computation is actually performed on a random sample from both input hypervolumes, constraining each to have the same point density given by the minimum of the point density of each input hypervolume, and the point density calculated using the volumes of each input hypervolume divided by \code{num.points.max}.

Because this algorithm is based on distances calculated between the distributions of random points, the critical distance (point density ^ (-1/n)) can be scaled by a user-specified factor to provide more or less liberal estimates (\code{distance_factor} greater than or less than 1).
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
              hypervolume_gaussian(x[, 
                c("bill_length_mm",
                "bill_depth_mm",
                "flipper_length_mm")],
              samples.per.point=75))

names(hv_list) <- names(penguins_no_na_split)
hv_list <- hypervolume_join(hv_list)

hv_occupancy <- hypervolume_n_occupancy(hv_list)
plot(hv_occupancy, cex.random = 1)

hv_occupancy_list_sex <- hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), each = 3))

plot(hv_occupancy_list_sex, cex.random = 1, show.density = FALSE)

}
