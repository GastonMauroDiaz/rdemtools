% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_voids.R
\name{filter_voids}
\alias{filter_voids}
\title{Filter objects by size}
\usage{
filter_voids(x, thr_size_in_pixels = 50, smallers_stay = TRUE)
}
\arguments{
\item{x}{\code{\linkS4class{RasterLayer}}. Binary image.}

\item{thr_size_in_pixels}{numeric. Threshold value, an object size in pixels.}

\item{smallers_stay}{logical}
}
\value{
\code{\linkS4class{RasterLayer}}
}
\description{
Filter objects by size
}
\examples{
r <- raster(ncol = 200, nrow = 100)
extent(r) <- extent(0, 200, 0, 100)
projection(r) <- projection(r) <- CRS("+init=epsg:32718")

###reference
set.seed(11)
reference <- fake_dem(r, n_random_data = 60, z_range = c(0, 600))
plot(reference)


p <- sampleRegular(reference, 10, sp = TRUE)
bigVoidsMask <- fake_voids(p, void_size = 20, reference)
p <- sampleRandom(reference, 10, sp = TRUE)
set.seed(1)
smallVoidsMask <- fake_voids(p, void_size = 3, reference, is_circular = TRUE)
voidsMask <- bigVoidsMask + smallVoidsMask
voidsMask <- voidsMask != 0

plot(voidsMask)

m <- filter_voids(voidsMask)
plot(m)
}
