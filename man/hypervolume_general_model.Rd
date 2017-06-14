\name{hypervolume_general_model}
\alias{hypervolume_general_model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates hypervolume by sampling from arbitrary model object.
}
\description{
Uses rejection sampling to obtain predicted values of a model object at uniformly random points within a range box, then converts output to a hypervolume.
}
\usage{
hypervolume_general_model(model, name = NULL, verbose = TRUE, 
  data = NULL, range.box = NULL, num.samples = ceiling(10^(3 + sqrt(ncol(range.box)))), 
  chunk.size = 10000, min.value = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
Any model object which can be used within a \code{predict(model, newdata, ...) call.}
}
  \item{name}{
Name of the output hypervolume
}
  \item{verbose}{
If \code{TRUE}, prints diagnostic output.
}
  \item{data}{
If not \code{NULL}, used to specify \code{range.box=padded_range(data)}.
}
  \item{range.box}{
A 2 x n matrix, where n is the number of dimensions of the data, and the first row corresponds to a lower limit and the second row to an upper limit. Each column is thus the low and high limits of the range box along each axis. Can be generated via \code{\link{padded_range}}.
}
  \item{num.samples}{
Number of samples to draw from the range box.
}
  \item{chunk.size}{
Number of samples to process in each \code{predict} call. Changing this value may affect the speed of function return but not the returned values.
}
  \item{min.value}{
If \code{TRUE}, discards sampled values below this threshold. Effectively used to set hypervolume boundaries.
}
  \item{\dots}{
Other arguments to be passed to \code{predict}, e.g. \code{type='response'}.
}
}

\value{
A \code{Hypervolume}-class object corresponding to retained values within the hyperbox of interest.
}
\examples{
data(iris)
iris[,"Species"] <- iris[,"Species"] == "setosa"
m_glm = glm(Species~.,data=iris)

hv_general_glm = hypervolume_general_model(m_glm, 
  range.box=padded_range(iris[,1:4]),type='response')
plot(hv_general_glm)	
	
}
