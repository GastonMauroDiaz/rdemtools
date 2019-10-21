#' @title  Voids are filling by fitting the fill source to the elevations
#'   surrounded the void
#'
#' @description \code{fitFillSource} returns a DEM with its voids filled with a
#'   source locally fitted with \code{\link{fit_trend_surface}}.
#'
#' @param x \code{\linkS4class{Raster}}. The dem with voids should be the first
#'   layer, the filling source or sources should be the other layer of layers.
#' @param show_progress_bar logical.
#' @inheritParams fit_trend_surface
#'
#' @details If the length of the shorter side of the bounding box of the void if
#'   less than \code{10}, the void is not filled.
#'
#'   To estimate the trend surface between DEM and fill source, the area of the
#'   bounding box is quadruplicate by side duplication (the void remain
#'   centered). The trend surface is used to correct the differences between DEM
#'   and fill source. If there are more than one filled source available for the
#'   void (i.e, \code{x} has more than two layers and them have data for the
#'   given void), the algorithm uses the median of the corrected difference to
#'   estimate which is the better choice to filling the void. To obtain the
#'   median, the algorithm only uses the pixels surrounding the void.
#'
#' @return \code{\linkS4class{RasterStack}}. First layer has the dem, second
#'   layer is a metadata layer that indicate which fill source was uses in each
#'   void.
#'
#' @references Grohman, G., Kroenung, G., Strebeck, J., 2006. Filling SRTM
#'   voids: the delta surface fill method. Photogramm. Eng. Remote Sensing 72,
#'   213-216.
#'
#' @example /inst/examples/fill_voids_g06_example.R
#'
#' @export fill_voids_g06
fill_voids_g06 <- function (x, np = 3, show_progress_bar = TRUE) {

    out  <- raster::subset(x, 1)
    meta <- raster(x)
    values(meta) <- NA
    kern <- EBImage::makeBrush(3, shape = "box")

    m <- raster::as.matrix(raster::subset(is.na(out), 1))
    m <- EBImage::bwlabel(m)

    if (show_progress_bar) pb <- pbCreate(max(m), progress="text")

    for (i in seq(length = max(m))){
      #print(i)
      #if (i == 369) browser()
      
      if (show_progress_bar) pbStep(pb, i)
      r <- raster(m)
      extent(r) <- extent(x)

      r[r != i] <- 0
      r <- r / i
      rc <- r[r == 1, drop = FALSE]
      xc <- crop(x, extent(rc) * 2)
      extended_rc <- extend(rc, extent(rc) * 2)
      extended_rc <- is.na(extended_rc)

      indices <- seq(from = 2, to = nlayers(xc))
      diff <- raster::subset(xc, subset=indices) -
              raster::subset(xc, 1)
      calc_diff <- function(index) {
        raster::subset(xc, subset=index) - raster::subset(xc, 1)
      }
      diff <- Map(calc_diff, indices)
      
      #indices <- unlist(Map(function(x) all(is.na(x[])), diff))

      no_of <- function(x){
          x[extended_rc] <- NA
          freq(is.na(x), value = 0)
        }
      no_of <- Map(no_of, diff)
      no_of <- unlist(no_of)

      indices <- no_of > 20 
      indices <- seq_along(diff)[indices]
      
      unlock <- FALSE
      diff <- stack(diff)
      if (any(indices) & ncol(diff) > 20 & nrow(diff) > 20) {
        unlock <- TRUE
        diff <- raster::subset(diff, subset=indices)
        xc <- raster::subset(xc, subset=indices + 1)
      }

      if (unlock) {

        trend <- fit_trend_surface(diff, size = ncell(diff), np)

        correctedDiff <- diff + trend

        b <- raster::as.matrix(is.na(raster::subset(xc, 1)))
        b <- EBImage::dilate(b, kern) - b
        b <- EBImage::bwlabel(b)
        b <- raster(b)
        extent(b) <- extent(xc)
        foo <- freq(b)
        index <- foo[which.max(foo[-1, 2]) + 1, 1]
        b <- b / index
        b[b != 1] <- 0
        b[b == 0] <- NA

        d <- correctedDiff * b

        fun <- function(x) {stats::median(x, na.rm = TRUE)}

        if (nlayers(diff) > 1) {
          d <- apply(values(d), FUN = fun, 2)
        } else {
          d <- fun(values(d))
        }

        index <- which.min(d)

        rc <- extend(rc, meta)
        meta <- cover(meta, rc * index)

        r <- raster::subset(trend, index)
        r <- raster::subset(xc, index + 1) - r

        rc <- crop(rc, r)
        r[is.na(rc)] <- NA
        r <- extend(r, out)

        out <- cover(out, r)
      }
    }
    if (show_progress_bar) pbClose(pb)
    out <- stack(out, meta)
    return(out)
}
