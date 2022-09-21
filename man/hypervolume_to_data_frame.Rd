\name{hypervolume_to_data_frame}
\alias{hypervolume_to_data_frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert hypervolumes to \code{data.frame}
}
\description{
Convert objects of class \code{Hypervolume} or \code{HypervolumeList} to a \code{data.frame}.
}
\usage{
hypervolume_to_data_frame(hv, remove_zeroes = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
A \code{Hypervolume} or \code{HypervolumeList}.
}
  \item{remove_zeroes}{
Remove zeroes from \code{ValuesAtRandomPoints}. See Details for further information. It works for \code{hypervolume_n_occupancy()}, \code{hypervolume_n_occupancy_test()}, \code{occupancy_to_union()}, \code{occupancy_to_intersection()} and \code{occupancy_to_unshared()}, otherwise ignored. 
}
}
\details{
Zero values are generated during the occupancy routine when a random point is included in some groups of hypervolumes but not in others. A tipical usage of \code{hypervolume_to_data_frame()} with objects generated with the occupancy routine should remove zeroes.
}
\value{
A \code{data.frame}.
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

# get the data.frame
hypervolume_to_data_frame(hv_list)
}
}
