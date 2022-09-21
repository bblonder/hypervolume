\name{find_optimal_occupancy_thin}
\alias{find_optimal_occupancy_thin}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Find optimal parameters to calculate occupancy
}
\description{
The \code{find_optimal_occupancy_thin()} function is used to find the optimal parameters for
\code{hypervolume_n_occupancy()}.
}
\usage{
find_optimal_occupancy_thin(..., 
                            verbose = TRUE, 
                            sequence = seq(0, 1, 0.1), 
                            n = 10, 
                            res_type = "raw")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
  Parameters to be used to run \code{hypervolume_n_occupancy()}.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
 }
  \item{sequence}{
Quantiles to be tested.
 }
  \item{n}{
Number of seeds to be tested.
 }
   \item{res_type}{
If \code{raw} print all the seeds and quantiles tested together with the resulting root mean square error (RMSE). If \code{summary} print RMSE mean and standard deviation for each quantile.
 }
 
}
\details{
The \code{find_optimal_occupancy_thin()} function searches for the optimal parameters for running \code{hypervolume_n_occupancy()}. It works by testing different quantiles and \code{n} seeds for random number generation (the same set of n seeds is tested for each quantile). RMSE is returned as the measure of the goodness of fit and results are ordered by increasing RMSE when \code{res_type = "raw"}. Quantile equal to 0 correspond to no thin. The obtained parameters can be used to feed arguments \code{quant.thin} and \code{seed} within the function \code{\link[=hypervolume_n_occupancy]{hypervolume_n_occupancy()}}.
}
\value{
A \code{data.frame}.
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
  hypervolume_gaussian(x[, c("bill_length_mm","bill_depth_mm","flipper_length_mm")],
                       samples.per.point=100, name = y), 
  x = penguins_no_na_split, 
  y = names(penguins_no_na_split))


# transform the list into an HypervolumeList
hv_list = hypervolume_join(hv_list)

# find optimal parameters
opt_par = find_optimal_occupancy_thin(hv_list, 
                                       classification = rep(c("female", "male"), 3),
                                       n = 20)

head(opt_par)

unoptimized_hv_occ = hypervolume_n_occupancy(hv_list, 
                        classification = rep(c("female", "male"), 3))

optimized_hv_occ = hypervolume_n_occupancy(hv_list, 
                        classification = rep(c("female", "male"), 3), 
                        quant.thin = opt_par[1, 2], seed = opt_par[1, 1])

}
}