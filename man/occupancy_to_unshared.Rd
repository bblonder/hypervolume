\name{occupancy_to_unshared}
\alias{occupancy_to_unshared}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Unshared fraction from an occupancy object
}
\description{
The \code{occupancy_to_unshared()} function is used to get the unshared fraction of hypervolumes of an object generated with the occupancy routine.
}
\usage{
occupancy_to_unshared(hv_list, method = "all", tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()},  \code{hypervolume_n_occupancy_test()}, \code{occupancy_to_union()}, \code{occupancy_to_intersection()} or \code{occupancy_filter}.
}
  \item{method}{
If \code{all} compute the unshared fraction of each hypervolume in \code{hv_list}. If \code{pairwise} compute the unshared fraction for each pairwise combination of hypervolumes in \code{hv_list}.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume. See details.
 }
}
\details{
Unshared fraction is the fraction of the hypervolume not shared with other hypervolumes. It is calculated from occupancy objects only (check \code{See Also}). When \code{method = "all"} a \code{HypervolumeList} containing the unshared fraction of each hypervolume is returned. When \code{method = "pairwise"} an \code{HypervolumeList} containing the unshared fraction of the pairwise combination of hypervolumes is returned. Hypervolumes generated when \code{method = "pairwise"} include the unshared fraction of both hypervolumes under comparison. The first of the two hypervolumes is assigned with  \code{ValueAtRandomPoints} equal to 1 while, the second is assigned with \code{ValueAtRandomPoints} equal to -1. This is useful when used in combination with \code{occupancy_filter()} or \code{hypervolume_to_data_frame()}. \cr
The \code{occupancy_to_unshared()} function attempts to reconstruct the volume of the unshared fraction from the \code{hv_list} provided by the user. At first, the volume of the union of hypervolumes is calculated for each hypervolume in \code{hv_list} as the the ratio between the total number of random points and the number of random points of the ith hypervolume of \code{hv_list}, multiplied by the volume of the ith hypervolume of \code{hv_list}. This step results in a number of reconstructed volumes equal to the number of hypervolumes in the jth bootstrapped occupancy_object. Reconstructed volumes are then compared among each other to ensure the consistency of the reconstruction. To do this, the distance among reconstructed volumes is calculated with the \code{dist()} function of the \code{stats} package. If at least one of the distances is greater than \code{tol} the computation is stopped and some suggestions are returned. The volume of the unshared fraction is then calculated as the ratio between the number of random points of the unshared fraction and the total number of random points, multiplied by the volume of the union of hypervolumes.
}

\value{
A \code{\link{Hypervolume-class}} or \code{\link{HypervolumeList-class}} object.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_intersection}},
\code{\link{occupancy_to_union}}, \code{\link{occupancy_filter}}
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


# transform the list into an HypervolumeList
hv_list = hypervolume_join(hv_list)

# calculate occupancy based on sex
hv_occupancy_list_sex = hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), 3))
                          
# get hypervolumes with the unshared fraction
hv_occupancy_sex_unshared <- occupancy_to_unshared(hv_occupancy_list_sex)
plot(hv_occupancy_sex_unshared)

# get hypervolumes with the unshared fraction between each pairwise combination of hypervolumes
hv_occupancy_sex_unshared_pw <- occupancy_to_unshared(hv_occupancy_list_sex)

# plot the unshared fraction with ggplot2
require(ggplot2)

# extract data to plot
occupancy_sex_pw_df <- hypervolume_to_data_frame(hv_occupancy_sex_unshared_pw)

ggplot(occupancy_sex_pw_df, aes(bill_length_mm, flipper_length_mm, col = Name)) +
  geom_point() +
  theme_bw()
}
}