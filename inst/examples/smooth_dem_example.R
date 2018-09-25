r <- raster()
r[] <- seq(1, ncell(r)) + rnorm(ncell(r), 20000, 10000)
plot(r)
r.smoothed <- smooth_dem(r, 10)
plot(r.smoothed)


