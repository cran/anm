\name{temp.era}
\alias{temp.era}
\title{Daily winter common EOF.}
\description{
  Common EOFs for daily December-February temperature. 
}
\usage{data(temp.era)}
\format{A list}
\value{
\item{EOF}{EOF patterns.}
\item{W}{Eigen values.} 
\item{PC}{Principal components of common PCA.} 
\item{n.fld}{Number of different predictors (see \link[clim.pact]{mixFields}).} 
\item{tot.var}{ Sum of all W squared.} 
\item{id.t}{Time labels for the fields (see \link[clim.pact]{catFields}) - used in \link[clim.pact]{DS}. }
\item{id.x}{Spatial labels for the fields (see \link[clim.pact]{mixFields}) - used in \link[clim.pact]{plotEOF}.}
\item{mon}{Month (1-12) [season (1-4) for daily data] to extract.}
\item{id.lon}{Spatial labels for the fields (see \link[clim.pact]{mixFields}) - used in \link[clim.pact]{plotEOF}.}
\item{id.lat}{Spatial labels for the fields (see \link[clim.pact]{mixFields}) - used in \link[clim.pact]{plotEOF}.}
\item{region}{Describes the region analysed.}
\item{tim}{Time information (usually redundant).}
\item{lon}{Longitudes associated with EOF patterns.} 
\item{lat}{Latitudes associated with EOF patterns.} 
\item{var.eof}{Fractional variances associated with EOF patterns.} 
\item{yy}{years.} 
\item{mm}{months.}
\item{dd}{days.}
\item{v.name}{Name of element.} 
\item{c.mon}{Month-season information.}
\item{f.name}{File name of original data.} 
}
\source{Rasmus E. Benestad \url{rasmus.benestad@met.no}.}
\references{
Reference to methodology: R.E. Benestad (2001), "A comparison between two empirical downscaling strategies", Int. J. Climatology, vol 210, pp.1645-1668. [DOI 10.1002/joc.703]. 
}
\keyword{datasets}


