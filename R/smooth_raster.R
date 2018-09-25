#' @title Smooth Raster Layer
#'
#' @description Wrapper function for \code{\link[fields]{image.smooth}}.
#'
#' @param x \code{\linkS4class{RasterLayer}}. Single Layer.
#' @inheritDotParams fields::image.smooth theta
#'
#' @return \code{\linkS4class{RasterLayer}}
#'
#' @example /inst/examples/smooth_dem_example.R
#' @export smooth_dem
smooth_dem <- function(x, theta) {
  xy <- data.frame(xyFromCell(x, seq(length=ncell(x))))
  v <- getValues(x)
  out <- fields::as.image(v, x=xy, nrow=ncol(x), ncol=nrow(x), na.rm=FALSE)
  rm(xy, v)
  out <- fields::image.smooth(out, theta=theta)
  out <- .xyzImageToRasterLayer(out)
  projection(out) <- projection(x)
  extent(out) <- extent(x)
  return(out)
}
