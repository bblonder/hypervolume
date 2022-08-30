\name{get_relative_volume}
\alias{get_relative_volume}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract the relative volume
}
\description{
Extract the relative volume from \code{hypervolume_n_occupancy()} \cr 
or \code{hypervolume_n_occupancy_test()}.
}
\usage{
get_relative_volume(hv_list, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{Hypervolume} or \code{HypervolumeList} object generated with \cr
  \code{hypervolume_n_occupancy()} or \code{hypervolume_n_occupancy_test()}.
}
  \item{tol}{
Set the tolerance for reconstructing whole volume.
 }
}
\details{
This function tries to reconstruct the volume of the union of hypervolumes from the object. For some reasons,
the reconstruction can fail. If the reconstruction fails, try to increase the tolerance or to check the accuracy of the \code{hypervolume_n_occupancy()} results.   
}
\value{
A named numeric vector with the relative volume of each input hypervolume
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