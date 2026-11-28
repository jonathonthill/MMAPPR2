% Manually synchronized with R/mmapprData.R.
\name{MmapprData-functions}
\alias{MmapprData-functions}
\alias{param}
\alias{param<-}
\alias{snpDistance}
\alias{peaks}
\alias{candidates}
\alias{param,MmapprData-method}
\alias{snpDistance,MmapprData-method}
\alias{peaks,MmapprData-method}
\alias{candidates,MmapprData-method}
\alias{param<-,MmapprData-method}
\title{MmapprData Getters and Setter}
\usage{
\S4method{param}{MmapprData}(obj)
\S4method{snpDistance}{MmapprData}(obj)
\S4method{peaks}{MmapprData}(obj)
\S4method{candidates}{MmapprData}(obj)
\S4method{param}{MmapprData}(obj) <- value
}
\arguments{
\item{obj}{A \code{MmapprData} object.}
\item{value}{A replacement \code{MmapprParam} object for \code{param<-}.}
}
\value{The requested slot; the replacement method returns the updated object.}
\description{Access the major MMAPPR2 result containers. Pipeline-generated result slots are read-only through the public API.}
\seealso{\code{\link{mmapprData}}}
