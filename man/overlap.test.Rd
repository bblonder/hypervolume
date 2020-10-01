\name{overlap.test}
\alias{overlap.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Null distribution for overap statistics
}
\description{
Generates null distribution of four different overlap statistics under null hypothesis that two samples are drawn from the same population.
calculates p value for observed value of each statistic.
}
\usage{
overlap.test(hv1, hv2, path, alternative = "one-sided", bins = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hv1}{
A hypervolume object
}
  \item{hv2}{
A hypervolume object
}
  \item{path}{
a path to a directory containing permuted hypervolumes, bootstrapped hypervolumes, or a vector of two paths to bootstrapped hypervolumes
}
  \item{alternative}{
"one-sided" or "two-sided"
}
}
\details{
overlap.test can generate a null distribution from the output of permute, resample with method = "bootstrap", or a vector of two bootstrap outputs. See examples for how to use each case.

The four overlap statistics are Sorensen, Jaccard, frac_unique_1, frac_unique_2. See \code{\link{hypervolume_overlap_statistics}}
}
\value{
  \item{p_values}{
a list of p_values indexed by the name of the relevant statistic
}
  \item{plots}{
a list of ggplot objects indexed by the name of the relevant statistic. The observed value of each statistic is represented as a dot on the x axis.
}
  \item{distribution}{
a matrix of random points sampled from the null distribution
}
}

\examples{
\dontrun{
# We will use the data in "quercus" as our population in this example
data("quercus")
# Consider taking two samples of size 150 from the population and you want to figure out whether
# the samples are similar by seeing if they occupy the same area in feature space.
qsample1 = quercus[sample(1:nrow(quercus), 150),]
qsample2 = quercus[sample(1:nrow(quercus), 150),]

# Construct two hypervolumes from the samples
hv1 = hypervolume(qsample1[,2:3])
hv2 = hypervolume(qsample2[,2:3])

# Approach 1
# Take 200 permutations of the 300 data points. Using more cores is faster.
perm_path = permute("Quercus_perm_150", hv1, hv2, n = 200, cores = 20)

# overlap.test takes perm_path as an input. Results show p value for the overlap statistics of hv1 and hv2 given null distribution generated from perm_path.
results1 = overlap.test(hv1, hv2, perm_path)

# Approach 2
# Under our null hypothesis the samples come from the same population. Approximate the original population by combining the data and 
# simulate drawing 150 data points 50 times
hv_combine = hypervolume(rbind(qsample1[,2:3],qsample2[,2:3]))
bootstrap_path = resample("Quercus_boot_150", hv_combine, method = "bootstrap", n = 50, points_per_resample = 150, cores = 20)

# overlap.test splits the 50 resampled hypervolumes in half and gets overlap statistic for each of the 25*25 pairs to generate the null 
# distribution. This method allows us to approximate the null distribution using 625 data points while only generating 50 hypervolumes as 
# opposed to permute which uses 400 hypervolumes to generate 200 data points.
results2 = overlap.test(hv1, hv2, bootstrap_path)

# Approach 3
# Suppose we have a size 300 sample and a size 150 sample and we want to know whether they come from the same distribution
qsample3 = quercus[sample(1:nrow(quercus), 300),]
hv3 = hypervolume(qsample3[,2:3])

# Permute still works in this case, however we can also use bootstrap by combining the data and drawing size 150 then size 300 samples
hv_combine = hypervolume(rbind(qsample1[,2:3],qsample3[,2:3]))
b150_path = resample("Quercus_150", hv_combine, method = "bootstrap", n = 25, points_per_resample = 150, cores = 20)
b300_path = resample("Quercus_300", hv_combine, method = "bootstrap", n = 25, points_per_resample = 300, cores = 20)

# overlap.test generates overlap statistics for each of the 25*25 possible pairs of size 150 and size 300 hypervolumes.
results3 = overlap.test(hv1, hv2, c(b150_path, b300_path))
}
}
