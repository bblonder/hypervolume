\name{Hypervolume-class}
\Rdversion{1.1}
\docType{class}
\alias{Hypervolume-class}
\title{Class \code{"Hypervolume"}}
\description{Primary storage class for stochastic descriptions of hypervolumes}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Hypervolume", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{Name}:}{Object of class \code{"character"} ~~ the name of the hypervolume}
    \item{\code{Method}:}{Object of class \code{"character"} ~~ the method used to construct this hypervolume}
    \item{\code{Data}:}{Object of class \code{"matrix"} ~~ May be empty if the hypervolume is not associated with data (e.g. convex expectation, set operations)}
    \item{\code{Dimensionality}:}{Object of class \code{"numeric"} ~~ Dimensionality of the hypervolume}
    \item{\code{Volume}:}{Object of class \code{"numeric"} ~~ Volume of the hypervolume}
    \item{\code{PointDensity}:}{Object of class \code{"numeric"} ~~ Number of random points per unit volume}
    \item{\code{Parameters}:}{Object of class \code{"list"} ~~ List of parameters that will depend on the method used to construct the hypervolume}
    \item{\code{RandomPoints}:}{Object of class \code{"matrix"} ~~ A matrix of uniformly random points distributed within the hypervolume}
    \item{\code{ValueAtRandomPoints}:}{Object of class \code{"numeric"} ~~ A vector of positive numbers representing the probabilty density at each random point in \code{@RandomPoints}}
  }
}