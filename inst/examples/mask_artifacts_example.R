r <- raster(ncol = 200, nrow = 100)
projection(r) <- NA
extent(r) <- extent(0, 200, 0, 100)

set.seed(11)
dem <- fake_dem(r, n_random_data = ncell(r) * 0.01, z_range = c(0, 500))

mask_artifacts(m)

plot(m)
