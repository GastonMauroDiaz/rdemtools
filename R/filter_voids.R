#' Filter objects by size
#'
#' @param x \code{\linkS4class{RasterLayer}}. Binary image. Voids equal 0, data equal 1.
#' @param thr_size_in_pixels numeric. Threshold value, an object size in pixels.
#' @param smallers_stay logical
#'
#' @return \code{\linkS4class{RasterLayer}}
#' @example /inst/examples/filter_voids_example.R
#' @export filter_voids
filter_voids <- function(x, thr_size_in_pixels = 50, smallers_stay = TRUE) {
  foo <- EBImage::bwlabel(raster::as.matrix(x))
  area <- EBImage::computeFeatures.shape(foo)
  area <- as.data.frame(area)
  foo <- raster::raster(foo)
  extent(foo) <- extent(x) # projection?
  if (smallers_stay) {
    index <- row.names(area)[area$s.area < thr_size_in_pixels]
  } else {
    index <- row.names(area)[area$s.area > thr_size_in_pixels]
  }
  index <- as.numeric(index)
  y <- data.frame(index, rep(1, length(index)))
  out <- raster::subs(foo, y)
  return(!is.na(out))
}

