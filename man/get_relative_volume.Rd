\name{get_relative_volume}
\alias{get_relative_volume}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the relative volume
}

\description{
The function \code{get_relative_volume()} computes the relative volume from objects generated with the occupancy routine.
}

\usage{
get_relative_volume(hv_list, tol = 1e-10)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{Hypervolume} or \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()}, \code{hypervolume_n_occupancy_permute()}, \code{hypervolume_n_occupancy_test()}, \code{occupancy_to_union()}, \code{occupancy_to_intersection()}, \code{occupancy_to_unshared()}, or \code{occupancy_filter()}.
}
  \item{tol}{
Set the tolerance for reconstructing whole volume. See details.
 }
}

\details{
The relative volume is calculated as the ratio between hypervolumes of an \code{HypervolumeList} and the volume resulting from the union of hypervolumes in the same \code{HypervolumeList}. Relative volumes can be calculated only for \code{HypervolumeList} generated with functions \code{hypervolume_n_occupancy()}, \code{hypervolume_n_occupancy_test()}, \code{hypervolume_n_occupancy_permute()}, \code{occupancy_to_union()}, \code{occupancy_to_ushared()}, \code{occupancy_to_intersection()} or \code{occupancy_filter()}. \cr
The \code{get_relative_volume()} function attempts to reconstruct the volume of the union of hypervolumes from \code{hv_list}. At first, the volume of the union of hypervolumes is calculated for each hypervolume of \code{hv_list} as the the ratio between the total number of random points and the number of random points of the ith hypervolume of \code{hv_list}, multiplied by the volume of the ith hypervolume \code{hv_list}. This step results in a number of reconstructed volumes equal to the number of hypervolumes in \code{hv_list}. Reconstructed volumes are then compared to ensure the consistency of the reconstruction. To do this, the distance among reconstructed volumes is calculated with the \code{dist()} function of the \code{stats} package. If at least one of the distances is greater than \code{tol} the computation is stopped and some suggestions are returned.
}

\value{
A named numeric vector with the relative volume of each input hypervolume
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_permute}}, \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_union}},
\code{\link{occupancy_to_unshared}}, \code{\link{occupancy_to_intersection}},
\code{\link{occupancy_filter}}
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

# calculate occupancy based on sex
hv_occupancy_list_sex = hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), 3))


# get the relative volume
get_relative_volume(hv_occupancy_list_sex)
}
}