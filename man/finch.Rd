\name{finch}
\alias{finch}
\docType{data}
\title{
Data and demo for Darwin's finches
}
\description{
Data for nine morphological traits for five species of Darwin's finches occurring on Isabela Island. Demonstration analysis of hypervolume overlaps between all pairs of species.
}
\usage{data(finch)}
\format{
  A data frame with 146 observations on the following 10 variables.
  \describe{
    \item{\code{Species}}{a factor with levels \code{Geospiza fortis fortis} \code{Geospiza fortis platyrhyncha} \code{Geospiza fuliginosa parvula} \code{Geospiza heliobates } \code{Geospiza prosthemelas prosthemelas}}
    \item{\code{BodyL}}{a numeric vector}
    \item{\code{WingL}}{a numeric vector}
    \item{\code{TailL}}{a numeric vector}
    \item{\code{BeakW}}{a numeric vector}
    \item{\code{BeakH}}{a numeric vector}
    \item{\code{LBeakL}}{a numeric vector}
    \item{\code{UBeakL}}{a numeric vector}
    \item{\code{N.UBkL}}{a numeric vector}
    \item{\code{TarsusL}}{a numeric vector}
  }
}
\source{
Data are from Snodgrass RE and Heller E (1904) Papers from the Hopkins-Stanford Galapagos Expedition, 1898-99. XVI. Birds. Proceedings of the Washington Academy of Sciences 5: 231-372.
}
\references{
Blonder, B., Lamanna, C., Violle, C., Enquist, B. The n-dimensional hypervolume. Global Ecology and Biogeography (in press)
}
\examples{
demo('finch', package='hypervolume')
}