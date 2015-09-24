\name{morphSnodgrassHeller}
\alias{morphSnodgrassHeller}
\docType{data}
\title{
Morphological data for Darwin's finches
}
\description{
Data for nine morphological traits for species of Darwin's finches occurring on the Galapagos Islands.

Note that the underlying morphological dataset has been augmented and improved since version 1.3.1 to include more species and islands. Results are not comparable to version 1.3.0 and below. To duplicate results in the Blonder et al. (2014) paper please install an older version of the package.
}
\usage{data("morphSnodgrassHeller")}
\format{
  A data frame with 549 observations on the following 20 variables.
  \describe{
    \item{\code{Source}}{a factor with levels \code{Snodgrass & Heller (1904)}}
    \item{\code{IslandID}}{a factor with levels \code{Balt_SS} \code{Drwn_Clp} \code{Esp_Hd} \code{Flor_Chrl} \code{Frn_Nrb} \code{Gnov_Twr} \code{Isa_Alb} \code{Mrch_Bndl} \code{Pnt_Abng} \code{Pnz_Dnc} \code{SCris_Chat} \code{SCru_Inde} \code{SFe_Brngt} \code{Snti_Jams} \code{Wlf_Wnm}}
    \item{\code{TaxonOrig}}{a factor with levels \code{Certhidea cinerascens bifasciata} \code{Certhidea cinerascens cinerascens} \code{Certhidea olivacea becki} \code{Certhidea olivacea fusca} \code{Certhidea olivacea luteola} \code{Certhidea olivacea mentalis} \code{Certhidea olivacea olivacea} \code{Geospiza affinis } \code{Geospiza conirostris conirostris} \code{Geospiza conirostris propinqua} \code{Geospiza crassirostris } \code{Geospiza fortis dubia} \code{Geospiza fortis fortis} \code{Geospiza fortis fratercula} \code{Geospiza fortis platyrhyncha} \code{Geospiza fuliginosa acutirostris} \code{Geospiza fuliginosa difficilis} \code{Geospiza fuliginosa fuliginosa} \code{Geospiza fuliginosa minor} \code{Geospiza fuliginosa parvula} \code{Geospiza habeli } \code{Geospiza heliobates } \code{Geospiza paupera } \code{Geospiza prosthemelas prosthemelas} \code{Geospiza prosthemelas salvini} \code{Geospiza psittacula psittacula} \code{Geospiza scandens abingdoni} \code{Geospiza scandens fatigata} \code{Geospiza scandens rothschildi} \code{Geospiza scandens scandens} \code{Geospiza septentrionalis } \code{Geospiza strenua }}
    \item{\code{GenusL69}}{a factor with levels \code{Camarhynchus} \code{Certhidea} \code{Geospiza} \code{Platyspiza}}
    \item{\code{SpeciesL69}}{a factor with levels \code{conirostris} \code{crassirostris} \code{difficilis} \code{fortis} \code{fuliginosa} \code{heliobates} \code{magnirostris} \code{olivacea} \code{parvulus} \code{pauper} \code{psittacula} \code{scandens}}
    \item{\code{SubspL69}}{a factor with levels \code{} \code{abingdoni} \code{affinis} \code{becki} \code{bifasciatus} \code{cinerascens} \code{conirostris} \code{darwini} \code{fusca} \code{habeli} \code{intermedia} \code{luteola} \code{mentalis} \code{olivacea} \code{parvulus} \code{propinqua} \code{psittacula} \code{rothschildi} \code{salvini} \code{scandens} \code{septentrionalis} \code{strenua}}
    \item{\code{SpeciesID}}{a factor with levels \code{Cam.hel} \code{Cam.par} \code{Cam.pau} \code{Cam.psi} \code{Cer.oli} \code{Geo.con} \code{Geo.dif} \code{Geo.for} \code{Geo.ful} \code{Geo.mag} \code{Geo.sca} \code{Pla.cra}}
    \item{\code{SubspID}}{a factor with levels \code{Cam.hel} \code{Cam.par.par} \code{Cam.par.sal} \code{Cam.pau} \code{Cam.psi.aff} \code{Cam.psi.hab} \code{Cam.psi.psi} \code{Cer.oli.bec} \code{Cer.oli.bif} \code{Cer.oli.cin} \code{Cer.oli.fus} \code{Cer.oli.lut} \code{Cer.oli.men} \code{Cer.oli.oli} \code{Geo.con.con} \code{Geo.con.dar} \code{Geo.con.pro} \code{Geo.dif.sep} \code{Geo.for} \code{Geo.ful} \code{Geo.mag.str} \code{Geo.sca.abi} \code{Geo.sca.int} \code{Geo.sca.rot} \code{Geo.sca.sca} \code{Pla.cra}}
    \item{\code{Sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{Plumage}}{a logical vector}
    \item{\code{BodyL}}{a numeric vector}
    \item{\code{WingL}}{a numeric vector}
    \item{\code{TailL}}{a numeric vector}
    \item{\code{BeakW}}{a numeric vector}
    \item{\code{BeakH}}{a numeric vector}
    \item{\code{LBeakL}}{a numeric vector}
    \item{\code{UBeakL}}{a numeric vector}
    \item{\code{N.UBkL}}{a factor with levels \code{} \code{10} \code{10.3} \code{10.5} \code{10.7} \code{11} \code{11.3} \code{11.5} \code{11.7} \code{12} \code{12.3} \code{12.5} \code{12.7} \code{13} \code{13.3} \code{13.5} \code{13.7} \code{14} \code{14.3} \code{14.5} \code{14.7} \code{15} \code{15.3} \code{15.5} \code{15.7} \code{16} \code{16.3} \code{16.5} \code{16.7} \code{17} \code{17.5} \code{6.5} \code{6.7} \code{7} \code{7.3} \code{7.5} \code{7.7} \code{8} \code{8..3} \code{8.3} \code{8.5} \code{8.7} \code{9} \code{9.3} \code{9.5} \code{9.7}}
    \item{\code{TarsusL}}{a numeric vector}
    \item{\code{MToeL}}{a logical vector}
  }
}
\source{
Snodgrass RE and Heller E (1904) Papers from the Hopkins-Stanford Galapagos Expedition, 1898-99. XVI. Birds. Proceedings of the Washington Academy of Sciences 5: 231-372. 

Downloaded from http://datadryad.org/resource/doi:10.5061/dryad.152
}

\examples{
data(morphSnodgrassHeller)
finch_isabela <- morphSnodgrassHeller[morphSnodgrassHeller$IslandID=="Isa_Alb",]
}
\keyword{datasets}
