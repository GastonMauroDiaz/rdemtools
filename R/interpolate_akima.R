#' @title Interpolate with package akima
#'
#' @description Wrapper function for \code{\link[akima]{interp}}.
#'
#' @param x \code{\linkS4class{RasterLayer}}. Single layer
#' @inheritDotParams akima::interp linear
#'
#' @return \code{\linkS4class{RasterLayer}}
#' @export interpolate_akima
#'
#' @example /inst/examples/interpolate_akima_example.R
interpolate_akima <- function(x, linear = TRUE) {

  out  <- x
  meta <- is.na(x)

  ma <- matrix(values(meta), ncol = ncol(x), nrow = nrow(x), byrow = TRUE)
  ma <- EBImage::bwlabel(ma)
  for (i in seq(length = max(ma))) {
    r <- raster(ma)
    extent(r) <- extent(x)
    r <- r == i
    r <- r[r == 1, drop = FALSE]

    xCropped <- raster::crop(x, extent(r) * 2)

    xy <- coordinates(xCropped)
    z <- values(xCropped)

    xy <- xy[!is.na(z),]
    z <- z[!is.na(z)]

    xo <- seq(xmin(xCropped), xmax(xCropped), length=ncol(xCropped))
    yo <- seq(ymin(xCropped), ymax(xCropped), length=nrow(xCropped))
    r <- akima::interp(xy[,1], xy[,2], z, xo=xo, yo=yo, linear)

    r <- raster(r)
    extent(r) <- extent(xCropped)
    r <- extend(r, x)
    projection(r) <- projection(out)
    out <- cover(out, r)
  }
  out <- stack(out, meta)
  return(out)
}

