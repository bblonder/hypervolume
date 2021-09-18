\name{hypervolume_n_occupancy_test}
\alias{hypervolume_n_occupancy_test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Significance of random points occupancy
}
\description{
Calculates overlap for two hypervolumes.
}
\usage{
hypervolume_n_occupancy_test(observed,
  path,
  alternative = "two_sided",
  CI = 0.95,
  cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observed}{
An HypervolumeList generated from \code{\link{hypervolume_n_occupancy}}.
}
  \item{path}{
A path to a directory of permuted hypervolumes generated with \code{\link{hypervolume_n_occupancy_permute}}.
}
  \item{alternative}{
Alternative hypothesis, can be one of \code{two_sided}, \code{more}, \code{less} or \code{more_less}.
}
  \item{CI}{
Desired confidence interval proportion.
}
  \item{cores}{
Number of logical cores to use while generating permuted hypervolumes. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in cores.
}
}
\details{
The observed difference between ValueAtRandomPoints of two groups is compared against null expectations generated with \code{hypervolume_n_occupancy_permute}. 
}
\value{
An HypervolumeList with length equal to the number of pairwise combination of the observed HypervolumeList elements. ValueAtRandomPoints are calculated as the difference between the ValueAtRandomPoints of the first and the second group for each pairwise combinations. Only significant values are retained according to CI.
}

\examples{
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))


penguins_no_na_split = split(penguins_no_na, 
              paste(penguins_no_na$species, penguins_no_na$sex, sep = "_"))


hv_list = lapply(penguins_no_na_split, function(x) 
      hypervolume_gaussian(x[, c("bill_length_mm","bill_depth_mm","flipper_length_mm")],
      samples.per.point=100))

names(hv_list) <- names(penguins_no_na_split)
hv_list <- hypervolume_join(hv_list)


hv_occupancy_list_sex <- hypervolume_n_occupancy(hv_list, 
                          classification = rep(c("female", "male"), each = 3))

# takes 99 permutations on 1 core
hyper_permuted <- hypervolume_n_occupancy_permute("permute", 
          hv_occupancy_list_sex, hv_list, classification = rep(c("female", "male"), each = 3), 
          n = 99, cores = 1)

hypervolume_test <- hypervolume_n_occupancy_test(hv_occupancy_list_sex, hyper_permuted, 
                      alternative = "more")


}
}