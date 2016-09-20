\name{expectation_maximal}
\alias{expectation_maximal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximal expectation
}
\description{
Creates a hypervolume from a set of points reflecting the maximal expectation. 
}
\usage{
expectation_maximal(input, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A dataset to be used as input to the \code{hypervolume} function
}
  \item{\dots}{
Arguments to the \code{hypervolume} function
}
}
\details{
This function is effectively an alias for the \code{hypervolume} function. You must decide what the maximal expectation is yourself!
}
\value{
A Hypervolume object.
}