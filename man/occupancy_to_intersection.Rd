\name{occupancy_to_intersection}
\alias{occupancy_to_intersection}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get the intersection of an occupancy object
}
\description{
The \code{occupancy_to_intersection()} function is used to get the intersection in an object generated with \code{hypervolume_n_occupancy()}.
}
\usage{
  occupancy_to_intersection(hv_list, method = "all", m = 2, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()}.
}
  \item{method}{
If \code{all} the function compares the target hypervolume with all the hypervolume in \code{hv}. If \code{n_wise} it computes the intersection of each n_wise combination of hypervolumes in \code{hv}.
 }
  \item{m}{
Number of elements to choose. Default to 2 (pairwise comparisons). This argument is ignored when \code{method} is set to \code{all}.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume.
 }
}

\value{
A \code{\link{Hypervolume-class}} or \code{\link{HypervolumeList-class}} object.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}, \cr \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_union}}, \cr
\code{\link{occupancy_to_unshared}} \code{\link{occupancy_filter}}
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
                    classification = rep(c("Adelie", "Chinstrap", "Gentoo"), 2))

# get the hypervolume of intersection
hv_occupancy_sex_intersection <- occupancy_to_intersection(hv_occupancy_list_sex)
plot(hv_occupancy_sex_intersection)

# get hypervolumes with the intersection among 3 or 2 combinations of hypervolumes
hv_occ_sex_intersection_3 <- occupancy_to_intersection(hv_occupancy_list_sex,
                                                       method = "n_wise",
                                                       m = 3)
hv_occ_intersection_2 <- occupancy_to_intersection(hv_occupancy_list_sex,
                                                   method = "n_wise",
                                                   m = 2)

}
}