\name{occupancy_to_union}
\alias{occupancy_to_union}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Union of hypervolumes from an occupancy object
}
\description{
The \code{occupancy_to_union()} function is used to get the union of all the hypervolumes or of n_wise combinations of them generated with \code{hypervolume_n_occupancy()}.
}
\usage{
  occupancy_to_union(hv_list, method = "all", m = 2, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()} or \code{hypervolume_n_occupancy_test()}.
}
  \item{method}{
If \code{all} the function performs the union of all the hypervolumes in \code{hv_list}. If \code{n_wise} it compute the union of each n_wise combination of hypervolumes in \code{hv}.
 }
   \item{m}{
Number of elements to choose. Default to 2 (pairwise comparisons). This argument is ignored when \code{method} is set to \code{all}.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume.
 }
}
\details{
Unshared fraction is the fraction of the hypervolume not shared with other hypervolumes.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}, \cr
\code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_intersection}}, \cr
\code{\link{occupancy_to_unshared}} \code{\link{occupancy_filter}}
}


\value{
A \code{\link{Hypervolume-class}} or \code{\link{HypervolumeList-class}} object.
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
                          
# get the union of all the hypervolumes
hv_occupancy_sex_union <- occupancy_to_union(hv_occupancy_list_sex)
plot(hv_occupancy_sex_union)

}
}