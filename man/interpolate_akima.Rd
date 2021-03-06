% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate_akima.R
\name{interpolate_akima}
\alias{interpolate_akima}
\title{Interpolate with package akima}
\usage{
interpolate_akima(x, linear = TRUE)
}
\arguments{
\item{x}{\code{\linkS4class{RasterLayer}}. Single layer}

\item{...}{Arguments passed on to \code{akima::interp}
\describe{
  \item{linear}{logical -- indicating wether linear or spline
    interpolation should be used.}
}}
}
\value{
\code{\linkS4class{RasterLayer}}
}
\description{
Wrapper function for \code{\link[akima]{interp}}.
}
\examples{
dem <- raster(system.file("external/dtm.grd", package="rdemtools"))
plot(dem)

xy <- xyFromCell(dem, sampleRegular(dem, 10, cells = TRUE)[, 1])
voids <- fake_voids(xy, 100, dem)
plot(voids)

dem[voids] <- NA
plot(dem)

dem <- interpolate_akima(dem)
plot(dem)
}
