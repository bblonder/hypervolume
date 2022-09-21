\name{hypervolume_n_occupancy_permute}
\alias{hypervolume_n_occupancy_permute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolumes through permuting labels of n pairwise groups of hypervolumes
}
\description{
Permute labels of an \code{hypervolume_n_occupancy()} object and calculate \code{hypervolume_n_occupancy()} for the permuted objects. This function is meant for taking a sample of all permutations and does not guarantee that permutations are not repeated. Newly generated hypervolume objects are written to file. This function is to be used within the occupancy routine.
}
\usage{
hypervolume_n_occupancy_permute(name, 
                                hv_list1,
                                hv_list2,
                                verbose = TRUE,
                                n = 9,
                                cores = 1)
}
\arguments{
  \item{name}{
File name; The function writes hypervolumes to file in ./Objects/<name>
}
  \item{hv_list1}{
An hypervolume list generated with \code{hypervolume_n_occupancy()}.
}
  \item{hv_list2}{
The hypervolume list used to generate \code{hv_list1}.
}
  \item{verbose}{
Logical value; outputs progress bar in console.
}
  \item{n}{
Number of permutations to take.
}
  \item{cores}{
Number of logical cores to use while generating permuted hypervolumes. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in \code{cores}.
}

}
\details{
\code{hypervolume_n_occupancy_permute()} creates a directory called Objects in the current working directory if a directory of that name doesn't already exist. Within this directory, it creates a directory for each pairwise combinations of elements within \code{hv_list1}. Group labels are permuted and a new \code{HypervolumeList} is saved as a rds file for each pairwise combination. IMPORTANT: only group labels are permuted, random points are kept fixed and will be the same across all the permuted hypervolumes. 

\code{hypervolume_n_occupancy_test()} takes in a \code{hypervolume_n_occupuancy_permute()} filepath output.

It is also possible to access the hypervolumes by using readRDS to read the hypervolume objects in one by one.
}


\section{Warning}{
  \code{hypervolume_n_occupancy_permute()} requires a lot of disk space especially when building occupancy hypervolumes with \code{method = "box"}. Try with a small number of replications and check the folder Objects for memory usage before to proceed.
}

\value{
Returns a string containing an absolute path equivalent to ./Objects/<name>
}
\seealso{
\code{\link{hypervolume_n_occupancy}}
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
  

names(hv_list) <- names(penguins_no_na_split)
hv_list <- hypervolume_join(hv_list)


hv_occupancy_list_sex <- hypervolume_n_occupancy(hv_list, 
                              classification = rep(c("female", "male"), each = 3))

# takes 9 permutations on 1 core
hypervolume_n_occupancy_permute("permute", hv_occupancy_list_sex,
                                hv_list , n = 9, cores = 1)



}
}