\name{occupancy_to_unshared}
\alias{occupancy_to_unshared}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Unshared fraction from an occupancy object
}
\description{
The \code{occupancy_to_unshared()} function is used to get the unshared fraction of each hypervolumes in an object generated with \code{hypervolume_n_occupancy()}.
}
\usage{
  occupancy_to_unshared(hv_list, method = "all", tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv_list}{
  A \code{HypervolumeList} object generated with \code{hypervolume_n_occupancy()}.
}
  \item{method}{
If \code{all} the function compares the target hypervolume with all the hypervolume in \code{hv_list}. If \code{pairwise} it compute the unshared fraction for each pairwise combination of hypervolumes in \code{hv}.
 }
  \item{tol}{
Set the tolerance for reconstructing whole volume.
 }
}
\details{
Unshared fraction is the fraction of the hypervolume not shared with other hypervolumes.
}
\value{
A \code{\link{Hypervolume-class}} or \code{\link{HypervolumeList-class}} object.
}

\seealso{
\code{\link{hypervolume_n_occupancy}}, \code{\link{hypervolume_n_occupancy_bootstrap}}, \cr \code{\link{hypervolume_n_occupancy_test}}, \code{\link{occupancy_to_intersection}}, \cr
\code{\link{occupancy_to_union}} \code{\link{occupancy_filter}}
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