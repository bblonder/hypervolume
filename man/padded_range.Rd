\name{padded_range}
\alias{padded_range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generates axis-wise range limits with padding
}
\description{
For each data axis, finds the minimum and maximum values. Then pads this range by a multiplicative factor of the range interval, and pads again by an additive amount.
}
\usage{
padded_range(data, multiply.interval.amount = 0, add.amount = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A m x n matrix whose range limits should be found.
}
  \item{multiply.interval.amount}{
A non-negative factor used to multiply the range interval. Can have either dimensionality 1 or n.
}
  \item{add.amount}{
A non-negative factor used to add to the range limits. Can have either dimensionality 1 or n.
}
}

\value{
A 2 x n matrix, whose first row is the low value along each axis and whose second row is the high value along each axis.
}
\examples{
data(iris)
iris_rangebox_nopadding = padded_range(iris[,1:3])
iris_rangebox_padding = padded_range(iris[,1:3], multiply.interval.amount=0.5, add.amount=0.1)
}