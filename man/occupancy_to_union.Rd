\name{occupancy_to_union}
\alias{occupancy_to_union}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Union of hypervolumes from an occupancy object
}
\description{
The \code{occupancy_to_union()} function is used to get the union of hypervolumes of an object generated with \code{hypervolume_n_occupancy()}.
}
\usage{
occupancy_to_union(hv_list, method = "all", m = 2, tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
A \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()}, \code{hypervolume_n_occupancy_test()}, \code{occupancy_to_intersection()}, \code{occupancy_to_unshared()} or \code{occupancy_filter()}.
}
  \item{method}{
If \code{all} compute the union of all the hypervolumes in \code{hv_list}. If \code{n_wise} compute of the union for each n_wise combination of hypervolumes in \code{hv_list}.
 }
   \item{m}{
Number of elements to choose. Default to 2 (pairwise comparisons). This argument is ignored when \code{method} is set to \code{all}.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume. See details.
 }
}

\details{
The \code{occupancy_to_union()} function takes as input a \code{HypervolumeList} generated with an occupancy function (check \code{See Also}) and returns a \code{Hypervolume} or \code{HypervolumeList} depending on \code{method}. When \code{method = "all"} the \code{occupancy_to_union()} function returns a \code{Hypervolume} representing the union of all the hypervolumes in \code{hv_list}. When \code{method = "n_wise"} a \code{HypervolumeList} in which each hypervolume represent the union of a combination of the hypervolumes in \code{hv_list} is returned. The number of hypervolumes for each combination is set with the argument \code{m}. Argument \code{m} can not be higher than the number of hypervolumes in \code{hv_list} and lower than 2. \cr
The \code{occupancy_to_union()} function attempts to reconstruct the volume of the union from the \code{hv_list} provided by the user. For each hypervolume in \code{hv_list}, it calculates the volume of the union as the ratio between the total number of random points and the number of random points of the ith hypervolume of \code{hv_list}, multiplied by the volume of the ith hypervolume of \code{hv_list}. This step results in a number of reconstructed volumes equal to the number of hypervolumes in the jth bootstrapped occupancy_object. Reconstructed volumes are then compared among each other to ensure the consistency of the reconstruction. To do this, the distance among reconstructed volumes is calculated with the \code{dist()} function of the \code{stats} package. If at least one of the distances is greater than \code{tol} the computation is stopped and some suggestions are returned.   
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_intersection}}, \code{\link{occupancy_to_unshared}}, \code{\link{occupancy_filter}}
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