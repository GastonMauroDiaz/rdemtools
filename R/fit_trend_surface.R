#' @title Fits a trend surface by least-squares to a RasterLayer.
#'
#' @details Fits a trend surface by least-squares to an object of class
#'   RasterLayer. This function is powered by \code{\link[spatial]{surf.ls}}.
#'
#' @param x \code{RasterObject}.
#'
#' @param sampleProportion numeric. Fraction. Approximate proportion of pixels
#'   that will be passed on to \code{\link[spatial]{surf.ls}} as xyz coordinates.
#'   \code{\link[raster]{sampleRandom}} is used with the argument \code{size} as
#'   \code{ncell(x) * sampleProportion}.
#'
#' @inheritParams spatial::surf.ls
#' @inheritParams raster::sampleRandom
#'
#' @return \code{RasterObject}.
#'
#' @example /inst/examples/fit_trend_surface_example.R
#'
#' @export fit_trend_surface
fit_trend_surface <- function (x, size = ncell(x), np = 3) {
    for (i in seq(length=nlayers(x))) {
      if (nlayers(x) > 1) {
        tmp <- sampleRandom(raster::subset(x, i), size,
          sp = TRUE)
      } else {
        tmp <- sampleRandom(x, size, sp = TRUE)
      }

      tmp <- cbind(tmp@coords, tmp@data)

      fit <- spatial::surf.ls(x = tmp[, 1], y = tmp[, 2], z = tmp[, 3], np)
      xl <- xmin(x)
      xu <- xmax(x)
      yl <- ymin(x)
      yu <- ymax(x)

      out <- spatial::trmat(fit, xl, xu, yl, yu, ncol(x))
      out <- raster(out)
      out <- resample(out, x)
      if (nlayers(x) > 1) {
        x[[i]] <- out
      } else {
        x <- out
      }
    }
    return(x)
}

