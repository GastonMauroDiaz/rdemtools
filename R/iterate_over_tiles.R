#' @title Apply a function in tiles
#'
#' @description Apply a function in tiles.
#'
#' @param x \code{\linkS4class{Raster}}
#' @param fun function
#' @param tile_size numeric. Tile size in pixels.
#' @param overlapping_percentage
#' @param show_progress_bar
#' @param filename character. To create a file during the process.
#' @inheritParams raster::focal
#'
#'
#' @details Useful if you are running out of \emph{RAM} while running spatial sensitive
#'   process.
#'
#' @return \code{\linkS4class{Raster}}
#' @export iterate_over_tiles
#'
#' @examples /inst/examples/iterate_over_tiles_example.R
iterate_over_tiles <- function(x, fun, tile_size=floor(sqrt(blockSize(x)$nrows[1] * ncol(x))),
                      overlapping_percentage = NULL, show_progress_bar = TRUE,
                      filename="", ...) {

  if (floor(tile_size) != tile_size) stop("Argument \"tile_size\" should be an integer")

  ncols <- rep(tile_size, floor(ncol(x)/tile_size))
  ncols <- c(ncols, ncol(x) - sum(ncols))
  col <- seq(1, ncol(x), tile_size)
  nrows<- rep(tile_size, floor(nrow(x)/tile_size))
  nrows <- c(nrows, nrow(x) - sum(nrows))
  row <- seq(1, nrow(x), tile_size)

  r <- raster(x)
  res(r) <- res(x) * tile_size
  if (!compareRaster(x, r, extent=TRUE, rowcol=FALSE, res=FALSE, stopiffalse=FALSE)) {
    a <- xmax(r) + res(r)[1]
    b <- ymin(r) - res(r)[1]
    e <- extent(xmin(r), a, b, ymax(r))
    r <- extend(r, e)
  }
  r <- crop(r, x, snap="out")
  values(r) <- 1:ncell(r)

  if (show_progress_bar) pb <- pbCreate(ncell(r), progress="text")

  rowCol <- matrix(c(sort(rep(1:nrow(r), ncol(r))), rep(1:ncol(r), nrow(r))), ncol=2)


	filename <- trim(filename)
  if(class(filename) !=  "character")
    stop("Argument \"filename\" should be an object of class \"character\"")

  if(nchar(filename) == 0) filename <- rasterTmpFile()

  writeRaster(x, filename, datatype="FLT4S", ...)
  x <- raster(filename)

  numberOfSquares <- r@data@max
  for (i in seq(length=numberOfSquares)) {
    r2 <- r[i, drop = FALSE]

    if (is.null(overlapping_percentage)) {
      x2 <- crop(x, r2)
    } else {
      stopifnot(class(overlapping_percentage) == "numeric")
      x2 <- crop(x, extent(r2) * (1 + overlapping_percentage/100))
    }

    if (!canProcessInMemory(x2, 3))
      stop(paste0("Decrease argument \"tile_size\", current value is ", tile_size, "."))

    x2 <- fun(x2, ...)

    if (!is.null(overlapping_percentage)) x2 <- crop(x2, r2)

    cell <- cellFromRowCol(x,row[rowCol[i, 1]], col[rowCol[i, 2]])

    for (w in seq(length=nrows[rowCol[i, 1]])) {
      if (w != 1) {
        cell <- cell + ncol(x)
        cells  <- c(cells, seq(cell, cell + ncols[rowCol[i, 2]] - 1))
      } else {
        cells  <- seq(cell, cell + ncols[rowCol[i, 2]] - 1)
      }
    }
    x <- update(x, v=values(x2), cell=cells)
    if (show_progress_bar) pbStep(pb, i)
  }
  if (show_progress_bar) pbClose(pb)
  return(x)
}
