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
hypervolume_svm(data, name = NULL, 
                  samples.per.point = ceiling((10^(3 + sqrt(ncol(data))))/nrow(data)), 
                  svm.nu = 0.01, svm.gamma = 0.5, 
                  scale.factor = 1,
                  chunk.size = 1000, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix or data frame, where m is the number of observations and n is the dimensionality.
}
  \item{name}{
A string to assign to the hypervolume for later output and plotting. Defaults to the name of the variable if NULL.
}
  \item{samples.per.point}{
Number of random points to be evaluated per data point in \code{data}.
}
  \item{svm.nu}{
A SVM parameter determining an upper bound on the fraction of training errors and a lower bound of the fraction of support vectors. Lower values result in tighter wrapping of the shape to the data (see section 2.2. of https://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.pdf). 

}
  \item{svm.gamma}{
A SVM parameter defining the inverse radius of influence of a single point. Low values yield large influences (smooth less complex wraps around the data) and high values yield small influences (tighter but potentially noiser wraps around the data) (see http://scikit-learn.org/stable/auto_examples/svm/plot_rbf_parameters.html).
}
  \item{scale.factor}{
A multiplicative factor used to determine the boundaries of the hyperelliptical sampling region. Larger values yield larger boundaries and can prevent clipping. Should not need to be changed in almost any situation.
}
  \item{chunk.size}{
Number of random points to process per internal step. Larger values may have better performance on machines with large amounts of free memory. Changing this parameter does not change the output of the function; only how this output is internally assembled.
}
  \item{verbose}{
Logical value; print diagnostic output if \code{TRUE}.
}
}

\value{
A \code{\link{Hypervolume-class}} object corresponding to the inferred hypervolume.
}

\seealso{
\code{\link{hypervolume_threshold}}
}
\examples{
data(iris)
hv = hypervolume_svm(data=subset(iris, Species=="setosa")[,1:2],name='setosa')
summary(hv)
}
