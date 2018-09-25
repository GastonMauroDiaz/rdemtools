dem <- raster(system.file("external/dtm.grd", package="rdemtools"))
plot(dem)

xy <- xyFromCell(dem, sampleRegular(dem, 10, cells = TRUE)[, 1])
voids <- fake_voids(xy, 100, dem)
plot(voids)

dem[voids] <- NA
plot(dem)

dem <- interpolate_akima(dem)
plot(dem)
