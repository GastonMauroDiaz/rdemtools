# gerenate the frame for the fake DEM
r <- raster(ncol = 200, nrow = 100)
projection(r) <- NA
extent(r) <- extent(0, 200, 0, 100)

# load a DEM to use as the true DEM
path <- system.file("external/dtm.grd", package="rdems")
dtm <- raster(path)

fdem <- fake_dem(r, dtm)
plot(fdem)
