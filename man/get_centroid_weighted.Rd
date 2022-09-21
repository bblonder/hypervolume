\name{get_centroid_weighted}
\alias{get_centroid_weighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get weighted centroid of hypervolume or hypervolume list
}
\description{
Returns the column weighted mean of the random points in each hypervolume. Useful for hypervolumes generated with \code{hypervolume_n_occupancy()}
or \code{hypervolume_n_occupancy_test()}.
}
\usage{
get_centroid_weighted(hv)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv}{
A \code{Hypervolume} or \code{HypervolumeList} object.
}
}
\details{
The function \code{get_centroid_weighted()} differs from \code{get_centroid()} because it uses occupancy values to weight random points for evaluating centroids position.
}

\value{
Either a vector or a matrix of column of centroid values along each axis. 
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_test}}
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

hv_list = hypervolume_join(hv_list)
hv_occupancy = hypervolume_n_occupancy(hv_list)

# unweighted centroids
get_centroid(hv_occupancy)

# weighted centroids
get_centroid_weighted(hv_occupancy)

}
}