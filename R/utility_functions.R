#' Get a fake DEM
#'
#' @param r \code{\linkS4class{RasterLayer}}
#' @param true_dem \code{\linkS4class{RasterLayer}}
#' @param n_random_data
#' @param z_range
#'
#' @return \code{\linkS4class{RasterLayer}}
#' @export fake_dem
#'
#' @example /inst/examples/fake_dem_example.R
fake_dem <- function(r = NULL, true_dem = NULL, n_random_data = ncell(r) * 0.01,
                       z_range = c(0, 500)) {

  if (!requireNamespace("automap", quietly = TRUE)) {
    stop("Package \"automap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("gstat", quietly = TRUE)) {
    stop("Package \"automap\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(r)) {
    # gerenate the frame for the fake DEM
    size <- 100
    r <- raster(ncol = size, nrow = size)
    projection(r) <- CRS("+init=epsg:32718")
    extent(r) <- extent(0, size, 0, size)
  }

  if (is.null(true_dem)) {
    # extract a variogram from a real DEM
    path <- system.file("external/dtm.grd", package="rdemtools")
    true_dem <- raster(path)

  }

  r[] <- 1
  names(r) <- "z"
  names(true_dem) <- "z"

  input_data <- sampleRandom(r, n_random_data, sp = TRUE)
  input_data@data$z <-
    min(z_range) + (runif(n_random_data, 0, 1) * (max(z_range) - min(z_range)))

  new_data <- as(r, "SpatialPixels")

  data_variogram <-
    sampleRandom(true_dem, 40, sp = TRUE)

  krigingResult <-
    suppressWarnings(automap::autoKrige(z ~ 1, input_data, new_data, data_variogram))

  raster(krigingResult$krige_output)
}


.getCircleEdgeXY <-
  function(center = c(0, 0),
           diameter = 1,
           npoints = 100) {
    radio = diameter / 2
    t <- seq(0, 2 * pi, length.out = npoints)
    x <- center[1] + radio * cos(t)
    y <- center[2] + radio * sin(t)
    return(data.frame(x, y))
  }


#' Get a fake mask of voids
#'
#' @param center dataframe or matrix. Coordinates.
#' @param void_size numeric
#' @param is_circular logical.
#' @inheritDotParams raster::rasterize y
#'
#' @return \code{\linkS4class{RasterLayer}}
#' @export fake_voids
#'
#' @example /inst/examples/fake_voids_example.R
fake_voids <- function(center, void_size, y, is_circular = FALSE){

  makeCircle <- function(x, i) {
    center <- as.numeric(x[, 1:2])
    void_size <- as.numeric(x[, 3])
    circle <- .getCircleEdgeXY(center = center, diameter = void_size)
    p1 <- Polygon(circle)
    Polygons(list(p1), i)
  }

  makePolygon <- function(x, i) {

    center <- as.numeric(x[, 1:2])
    diameter <- as.numeric(x[, 3])
    foo <- function(x) .getCircleEdgeXY(center, x)
    concentricCircles <- Map(foo, (diameter/2):diameter)
    concentricCircles <- data.frame(concentricCircles)
    concentricCircles <-
      data.frame(x = unlist(concentricCircles[, seq(1, ncol(concentricCircles), 2)]),
                 y = unlist(concentricCircles[, seq(2, ncol(concentricCircles), 2)]))

    index <- sample(1:nrow(concentricCircles), nrow(concentricCircles) * 0.01)
    concentricCircles <- concentricCircles[index,]

    index <- chull(concentricCircles)
    index <- c(index, index[1])

    p1 <- Polygon(concentricCircles[index,])
    Polygons(list(p1), i)
  }

  ds <- data.frame(sp::coordinates(center), void_size)
  ds <- Map(function(x) ds[x,], 1:nrow(ds)) # to list

  if (is_circular) {
    p2 <- Map(makeCircle, ds, 1:length(ds))
  } else {
    p2 <- Map(makePolygon, ds, 1:length(ds))
  }

  p <- SpatialPolygons(p2, proj4string=CRS(projection(y)))
  !is.na(rasterize(p, y))
}



.xyzImageToRasterLayer <- function(x) {
  xy <- expand.grid(x$x, x$y)
  r <- raster(nrows=length(x$y), ncols=length(x$x),
              xmn=min(x$x), xmx=max(x$x), ymn=min(x$y), ymx=max(x$y), crs=NA)
  cell <- cellFromXY(r, xy)
  r[cell] <- as.numeric(x$z)
  return(r)
}

.filesWithoutAux <- function(ext=".tif") {
  files <- dir(pattern=ext)
  tmp <- grep(paste0(ext, "."), files)
  if (length(tmp) > 0) files <- files[-tmp]
  return(files)
}

.makeF8single <- function(a, ...) { # single layer output

  function(x, filename = "", ...) {
    # code from raster-package vignette
    # "Writing functions for large raster files"
    # function referred as f8
    out <- raster(x)
    big <- ! canProcessInMemory(out, 3)
    filename <- trim(filename)
    if (big & filename == '') {
      filename <- rasterTmpFile()
    }
    if (filename != '') {
      out <- writeStart(out, filename, ...)
      todisk <- TRUE
    } else {
      vv <- matrix(ncol=nrow(out), nrow=ncol(out))
      todisk <- FALSE
    }

    bs <- blockSize(x)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...) # new code
        out <- writeValues(out, v, bs$row[i])
        pbStep(pb, i)
      }
      out <- writeStop(out)
    } else {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...) # new code
        cols <- bs$row[i]:(bs$row[i]+bs$nrows[i]-1)
        vv[,cols] <- matrix(v, nrow=out@ncols)
        pbStep(pb, i)
      }
      out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)
    return(out)
  }
}

.makeF8multi <- function(a, ...) { # multi layer output

  function(x, filename = "", ...) {
    # based on code from raster-package vignette
    # "Writing functions for large raster files"
    # function referred as f8
    out <- brick(x) #
    big <- ! canProcessInMemory(out, 3)
    filename <- trim(filename)
    if (big & filename == '') {
      filename <- rasterTmpFile()
    }
    if (filename != '') {
      out <- writeStart(out, filename, ...)
      todisk <- TRUE
    } else {
      vv <- array(dim = c(ncol(out), nrow(out), nlayers(out)))
      todisk <- FALSE
    }

    bs <- blockSize(x)
    pb <- pbCreate(bs$n, ...)

    if (todisk) {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...)
        out <- writeValues(out, v, bs$row[i])
        pbStep(pb, i)
      }
      out <- writeStop(out)
    } else {
      for (i in 1:bs$n) {
        v <- getValues(x, row=bs$row[i], nrows=bs$nrows[i] )
        v <- a(v, ...)
        cols <- bs$row[i]:(bs$row[i]+bs$nrows[i]-1)
        vv[,cols,] <- array(v, dim=c(bs$nrows[i],nrow=out@ncols, nlayers(x)))
        pbStep(pb, i)
      }
      out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)
    return(out)
  }
}
