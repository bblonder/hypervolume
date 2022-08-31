\name{occupancy_filter}
\alias{occupancy_filter}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Subset occupancy hypervolumes
}
\description{
The \code{occupancy_filter()} function is used to subset an hypervolume generated with \cr
\code{hypervolume_n_occupancy()} or \code{hypervolume_n_occupancy_test()}.
}
\usage{
  occupancy_filter(hv, operator = NULL, filter = NULL, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
  A \code{Hypervolume} or \code{HypervolumeList} object generated with  \cr
  \code{hypervolume_n_occupancy()} or \code{hypervolume_n_occupancy_test()}.
}
  \item{operator}{
Binary operator which allow the comparison.
 }
  \item{filter}{
Threshold value to perform the operation.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume.
 }
}
\details{
This function set the occupancy value based on the user-provided operation. Volume of the hypervolumes are changed accordingly.  
}
\value{
A \code{\link{Hypervolume-class}} or \code{\link{HypervolumeList-class}} object.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_test}}
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


# set to 0 values lower than 0.35
occupancy_filter(hv_occupancy_list_sex, operator = "<", filter = "0.35")
}
}