\name{hypervolume_permute}
\alias{hypervolume_permute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolumes through permuting data of two hypervolumes
}
\description{
Takes two data of two hypervolume objects (with the same column labels) and generates pairs of hypervolumes with the original sizes of the input hypervolumes but with permuted data (the rows of the original hypervolumes' data are combined and redistributed to the two new hypervolumes). This function is meant for taking a sample of all permutations and does not guarantee that permutations are not repeated. Newly generated hypervolume objects are written to file.
}
\usage{
hypervolume_permute(name, 
                    hv1, 
                    hv2, 
                    n = 50, 
                    cores = 1, 
                    verbose = TRUE)
}
\arguments{
  \item{name}{
File name; The function writes hypervolumes to file in ./Objects/<name>
}
  \item{hv1}{
A hypervolume object
}
  \item{hv2}{
A hypervolume object
}
  \item{n}{
number of permutations to take
}
  \item{cores}{
Number of logical cores to use while generating permuted hypervolumes. If parallel backend already registered to \code{doParallel}, function will use that backend and ignore the argument in cores.
}
  \item{verbose}{
Logical value; If function is being run sequentially, outputs progress bar in console.
}
}
\details{
\code{hypervolume_permute} creates a directory called Objects in the current working directory if a directory of that name doesn't already exist. Returns an absolute path to directory with permuted hypervolumes. rds files are stored in separate subdirectories for each permutation. Use \code{hypervolume_permute} when generating null distribution of overlap statistics. \code{\link{hypervolume_overlap_test}} takes in a \code{hypervolume_permute} filepath output.

It is also possible to access the hypervolumes by using readRDS to read the hypervolume objects in one by one.
}

\value{
returns a string containing an absolute path equivalent to ./Objects/<name>
}
\seealso{
\code{hypervolume_overlap_test}
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
perm_path = hypervolume_permute("Quercus_perm_150", hv1, hv2, n = 200, cores = 20)

# hypervolume_overlap_test takes perm_path as an input.
# Results include p value for the overlap statistics of hv1 and hv2 as well as
# null distribution generated from perm_path. The null distribution assumes data
# for hv1 and hv2 are drawn from the same distribution and permuting data will 
# not change the overlap statitsics.
results = hypervolume_overlap_test(hv1, hv2, perm_path)
}
}