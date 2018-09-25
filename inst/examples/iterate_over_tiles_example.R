r <- raster(ncol = 200, nrow = 100)
projection(r) <- NA
extent(r) <- extent(0, 200, 0, 100)

set.seed(11)
dem <- fake_dem(r, n_random_data = ncell(r) * 0.01, z_range = c(0, 500))

m <-
  iterate_over_tiles(
    dem,
    mask_artifacts,
    tile_size = 40,
    overlapping_percentage = 20,
    show_progress_bar = TRUE
  )

plot(m)

