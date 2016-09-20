\name{hypervolume_svm}
\alias{hypervolume_svm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hypervolume construction via one-class support vector machine (SVM) learning model
}
\description{
Constructs a hypervolume by building a one-class support vector machine that classifies data points as 'in' and other locations as 'out'. This is accomplished by 1) transforming the input data into a high-dimensional nonlinear space in which the data points can be optimally separated from background by a single hyperplane, 2) back-transforming the hyperplane into the original space, 3) delineating an adaptive grid of random points near the original data points, and 4) using the SVM to predict if each of these points is in or out.
}
\usage{
hypervolume_svm(data, 
	name = NULL, verbose = TRUE, 
	output.density = 10^ncol(data), 
	expectation.num.shifts = 1, expectation.bin.widths = 2*estimate_bandwidth(data), 
	svm.nu = 0.01, svm.gamma = 0.5, svm.chunksize = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
  \item{output.density}{
Point density of output hypervolume.
}
  \item{expectation.num.shifts}{
The integer number of empty grid cells away from each filled grid cell that are retained in the adaptive grid. Should not need to be changed.
}
  \item{expectation.bin.widths}{
The spacing between grid points for the adaptive grid. Should not need to be changed.
}
  \item{svm.nu}{
A SVM parameter determining an upper bound on the fraction of training errors and a lower bound of the fraction of support vectors. Lower values result in tighter wrapping of the shape to the data (see section 2.2. of https://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf). 

}
  \item{svm.gamma}{
A SVM parameter defining the inverse radius of influence of a single point. Low values yield large influences (smooth less complex wraps around the data) and high values yield small influences (tighter but potentially noiser wraps around the data) (see http://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html).
}
  \item{svm.chunksize}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
}

\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}

\seealso{
\code{\link{hypervolume_quantile_threshold}}, \code{\link{expectation_adaptive_box}}
}
\examples{
data(iris)
hv = hypervolume_svm(data=subset(iris, Species=="setosa")[,1:3],name='setosa')
summary(hv)
}
