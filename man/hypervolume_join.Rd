\name{hypervolume_join}
\alias{hypervolume_join}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Concatenate hypervolumes
}
\description{
Combines multiple hypervolumes or hypervolume lists into a single HypervolumeList suitable for analysis or plotting.
}
\usage{
hypervolume_join(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
One or more objects of class \code{Hypervolume} or \code{HypervolumeList}, or a \code{list()} of \code{Hypervolume} objects.
}
}

\value{
A \code{HypervolumeList} containing all hypervolumes in all arguments.
}

\examples{
data(penguins,package='palmerpenguins')
penguins_no_na = as.data.frame(na.omit(penguins))
penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                    c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
penguins_chinstrap = penguins_no_na[penguins_no_na$species=="Chinstrap",
                      c("bill_length_mm","bill_depth_mm","flipper_length_mm")]

hv1 = hypervolume_box(penguins_adelie,name='Adelie')
hv2 = hypervolume_box(penguins_chinstrap,name='Chinstrap')

hvs_joined = hypervolume_join(hv1, hv2)
}
