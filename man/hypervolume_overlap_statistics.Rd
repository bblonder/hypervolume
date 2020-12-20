\name{hypervolume_overlap_statistics}
\alias{hypervolume_overlap_statistics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Overlap statistics for set operations (Sorensen, Jaccard, etc.)
}
\description{
Calculates overlap metrics for two hypervolumes
}
\usage{
hypervolume_overlap_statistics(hvlist)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hvlist}{
A set of hypervolumes calculated from \code{\link{hypervolume_set}}
}
}
\value{
A set of multiple metrics
\item{\code{jaccard} }{Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)}
\item{\code{sorensen} }{Sorensen similarity (twice the volume of intersection of 1 and 2 divided by volume of 1 plus volume of 2)}
\item{\code{frac_unique_1} }{Unique fraction 1 (volume of unique component of 1 divided by volume of 1))}
\item{\code{frac_unique_2} }{Unique fraction 2 (volume of unique component of 2 divided by volume of 2))}
}

\examples{
\dontrun{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
penguins_chinstrap = penguins_no_na[penguins_no_na$species=="Chinstrap",
                      c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv1 = hypervolume_box(penguins_adelie,name='Adelie')
hv2 = hypervolume_box(penguins_chinstrap,name='Chinstrap')

hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)

hypervolume_overlap_statistics(hv_set)
}
}