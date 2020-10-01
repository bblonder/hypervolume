\name{permute}
\alias{permute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolumes through permuting data of two hypervolumes
}
\description{
Takes two data of two hypervolume objects (with same number of columns) and generates pairs of hypervolumes with the original sizes of the input hypervolumes but with permuted data. This function is meant for taking a sample of all permutations and does not guarantee that permutations are not repeated.
}
\usage{
permute <- function(name, hv1, hv2, n = 50, cores = 1, verbose = TRUE)
}
\arguments{
  \item{name}{
File name; The function writes hypervolumes to file in ./Objects/<name>
}
  \item{hv1}{
A n-dimensional hypervolume
}
  \item{hv2}{
A n-dimensional hypervolume
}
  \item{n}{
number of permutations to take
}
  \item{cores}{
Number of logical cores to use while generating permuted hypervolumes. If parallel backend already registered to doParallel, function will use that backend and ignore the argument in cores.
}
  \item{verbose}{
Logical value; If function is being run sequentially, outputs progress bar in console.
}
}
\details{
\code{permute} creates a directory called Objects in the current working directory if a directory of that name doesn't already exist. Returns an absolute path to directory with permuted hypervolumes. rds files are stored in separate subdirectories for each permutation. Use permute when generating null distribution of overlap statistics. \code{\link{overlap.test}} takes in a \code{permute} filepath output.

It is also possible to access the hypervolumes by using readRDS to read the hypervolume objects in one by one.
}

\value{
returns a string containing an absolute path equivalent to ./Objects/<name>
}
\examples{
\dontrun{
data("quercus")
# For this example consider taking two samples of size 150 from the data.
qsample1 = quercus[sample(1:nrow(quercus), 150),]
qsample2 = quercus[sample(1:nrow(quercus), 150),]

# Construct two hypervolumes from the samples
hv1 = hypervolume(qsample1[,2:3])
hv2 = hypervolume(qsample2[,2:3])

# Take 200 permutations of the 300 data points. Using more cores is faster.
perm_path = permute("Quercus_perm_150", hv1, hv2, n = 200, cores = 20)

# overlap.test takes perm_path as an input. Results show p value for the overlap statistics of hv1 and hv2 given null distribution generated from perm_path.
results = overlap.test(hv1, hv2, perm_path)
}
}