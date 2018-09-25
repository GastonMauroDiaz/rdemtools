#' Mask artifacts
#'
#' @param x \code{linkS4class{RasterLayer}}
#' @param kernel_sizes numeric. Size of the kernels, could be one or more. See details.
#' @param times_sd numeric. Thresold value. See details.
#'
#' @return x \code{linkS4class{RasterLayer}}. Binary image.
#' @export mask_artifacts
#'
#' @details Is a faster approach that I pick up from Hengl and Reuter (2009).
#'   The algorithm compute the average elevation with a kernel and compute the
#'   difference between the average elevation and the elevation at the core
#'   cell. A threshold is used to classify this difference.
#'
#'   If users uses several kernels, the masks are aggregated.
#'
#' @references Hengl, T., Reuter, H.I., 2009. GEOMORPHOMETRY Concepts, Software,
#'   Applications, Developmen. ed. Elsevier, Amsterdam. (p. 100)
#'
#' @example /inst/examples/mask_artifacts_example.R
#'
mask_artifacts <- function(x, kernel_sizes = 5, times_sd = 3) {
  for (i in seq(length = length(kernel_sizes))) {
    kernel_size <- kernel_sizes[i]
    f <- focal(x, matrix(1, kernel_size, kernel_size), fun = mean, na.rm = TRUE)
    s <- focal(x, matrix(1, kernel_size, kernel_size), fun = sd, na.rm = TRUE)
    d <- x - f
    s <- s * times_sd
    if(i == 1) {out <- abs(d) > s} else {out <- out + (abs(d) > s)}
  }
  out[out > 0] <- 1
  return(out)
}
