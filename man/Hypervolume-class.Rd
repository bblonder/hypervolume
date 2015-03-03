\name{Hypervolume-class}
\Rdversion{1.1}
\docType{class}
\alias{Hypervolume-class}

\title{Class hypervolume}
\description{
A class used to store a stochastic description of a hypervolume.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Hypervolume", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"}. A string naming the hypervolume, used in plotting. }
    \item{\code{Data}:}{Object of class \code{"matrix"}. If available, the raw data used to construct the hypervolume. Defaults to a one-row NaN vector for hypervolumes returned by set operations.   }
    \item{\code{Dimensionality}:}{Object of class \code{"numeric"}. The dimensionality of the hypervolume. }
    \item{\code{Volume}:}{Object of class \code{"numeric"}. The volume of the hypervolume, in units of the product of all dimensions. }
    \item{\code{PointDensity}:}{Object of class \code{"numeric"}. The number density of the uniformly sampled random points characterizing the hypervolume. }
    \item{\code{Bandwidth}:}{Object of class \code{"numeric"}. If available, the bandwidth vector used to construct the hypervolume. Defaults to a one-row NaN vector for hypervolumes returned by set operations. }
    \item{\code{DisjunctFactor}:}{Object of class \code{"numeric"}. The ratio of the inferred volume to the volume of a hypervolume constructed from the same data with disjunct data points (i.e. no kernels overlap). Varies from zero to one. High values suggest that bandwidth should be increased. }
    \item{\code{RepsPerPoint}:}{Object of class \code{"numeric"}. If available, the number of random points used per observation to construct the hypervolume. Defaults to NaN for hypervolumes returned by set operations. }
    \item{\code{QuantileThresholdDesired}:}{Object of class \code{"numeric"}. If available, the quantile requested by the user and used to construct the hypervolume. Defaults to NaN for hypervolumes returned by set operations. }
    \item{\code{QuantileThresholdObtained}:}{Object of class \code{"numeric"}. If available, the quantile obtained by the hypervolume algorithm. Defaults to NaN for hypervolumes returned by set operations. }
    \item{\code{RandomUniformPointsThresholded}:}{Object of class \code{"matrix"} A set of uniformly random points guaranteed to be in the hypervolume. }
    \item{\code{ProbabilityDensityAtRandomUniformPoints}:}{Object of class \code{"numeric"} A vector of integers proportional to the probability density at each uniformly random point in the hypervolume. Defaults to a 1-valued vector for hypervolumes returned by set operations because set operations are well defined for volumes and not for probability density functions. }
  }
}
\section{Methods}{
Summary and plot methods are available.
}
