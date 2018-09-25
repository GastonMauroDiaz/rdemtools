.karkeeetal <- function(hNgb, alpha, theta, weights, rx, ry) {
  d1 <- sqrt(sum(rx^2, ry^2))
  d <- c(d1, ry, d1, rx, 0, rx, d1, ry, d1)
  beta <- c(7/4, 0, 1/4, 3/2, NA, 1/2, 5/4, 1, 3/4) * pi
  deltaH <- c()
  for (i in seq(length=length(hNgb))) {
    deltaH[i] <- cos(alpha[i] - beta[i]) * tan(theta[i]) * d[i]  # eq. 3
  }
  h <- hNgb + deltaH  # eq. 4

  if (sum(is.na(h)) != length(h)) {
    weightsTemp <- weights[!is.na(h)]
    weight  <- 1 / (sum(1 / weightsTemp) / length(weightsTemp) + 1)  # eq. 6
    out <- sum(h * weights, na.rm=TRUE) / sum(weightsTemp)
    out <- c(out, weight)
  } else {
    out <- c(NA, NA)
  }

  return(out)
}

#' @title Void filling using derivatives
#'
#' @description Void filling using trigonometric, elevation of pixels
#'   surrounding the voids, and slope and aspect from the fill source.
#'
#' @param x \code{\linkS4class{RasterLayer}}. DEM with voids.It must be in a
#'   cartographic projection with the same unit that the elevation data (meters
#'   is strongly recommended).
#' @param asp \code{\linkS4class{RasterLayer}}. Aspect in radians derivate from the fill source.
#' @param slp \code{\linkS4class{RasterLayer}}. Slope in radians derivate from the fill source.
#' @param show_progress_bar logical.
#'
#' @details This function implement the algorithm developed by Karkee et al.
#'   (2008) (see references). The filling process iteratively grows from the
#'   surrounding pixels to the center of the void.
#'
#' @example /inst/examples/fill_voids_k08_example.R
#'
#' @export fill_voids_k08
#'
#' @references Karkee, M., Steward, B., Aziz, S., 2008. Improving quality of
#'   public domain digital elevation models through data fusion. Biosyst. Eng.
#'   101, 293-305. doi:10.1016/j.biosystemseng.2008.09.010
fill_voids_k08 <- function (x, slp, asp, show_progress_bar = TRUE) {

    rx <- res(x)[1]
    ry <- res(x)[2]
    compareRaster(asp, slp)
    compareRaster(x, asp)

    weight <- raster(x)
    values(weight) <- 1
    weight[is.na(x)]  <- NA
    kern = EBImage::makeBrush(3, shape='box')

    ies <- seq(length=ncell(x))[is.na(values(x))]
    if (show_progress_bar) pb <- pbCreate(length(ies), progress="text")
    count <- 0

    on <- TRUE
    while(on) {
      tmp <- raster::as.matrix(is.na(weight))
      tmp <- tmp - EBImage::erode(tmp, kern)
      tmp <- raster(tmp)
      extent(tmp) <- extent(weight)

      if(freq(tmp, v=1)) {
        ies <- seq(length=ncell(x))[values(tmp) == 1]
        for(i in ies){
          count <- count + 1
          xy <- xyFromCell(x, i)
          xy <- matrix(c(c(-1,0,1,-1,0,1,-1,0,1) * res(x)[1] + xy[1],
                    c(1,1,1,0,0,0,-1,-1,-1) * res(x)[2] + xy[2]), ncol = 2)
          p <- sp::SpatialPoints(xy, proj4string=CRS(projection(x)))
          hNgb <- extract(x, p)
          alpha <- extract(asp, p)
          theta <- extract(slp, p)
          weights <- extract(weight, p)
          out <- .karkeeetal(hNgb, alpha, theta, weights, rx, ry)
          x[i] <- out[1]
          weight[i] <- out[2]
          if (show_progress_bar) pbStep(pb, count)
        }
      } else {on <- FALSE}
    }
    if (show_progress_bar) pbClose(pb)
    return(x)
}
