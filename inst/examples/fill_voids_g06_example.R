r <- raster(ncol = 200, nrow = 100)
extent(r) <- extent(0, 200, 0, 100)
projection(r) <- projection(r) <- CRS("+init=epsg:32718")

###reference
set.seed(11)
reference <- fake_dem(r, n_random_data = 60, z_range = c(0, 600))
plot(reference)

set.seed(18)
fillSource <- fake_dem(r, z_range = c(0, 20))
fillSource <- fillSource - mean(fillSource[])
fillSource <- fillSource + reference
fillSource <- smooth_dem(fillSource, theta = 6)
plot(fillSource)

upperHalfCells <- 1:(ncell(r)/2)
fillSource1 <- fillSource
fillSource1[upperHalfCells] <- (fillSource * 2)[upperHalfCells]
fillSource2 <- fillSource * 2
fillSource2[upperHalfCells] <- fillSource[upperHalfCells]

p <- sampleRegular(reference, 10, sp = TRUE)
voidsMask <- fake_voids(p, void_size = 20, reference)

dem <- reference
dem[voidsMask] <- NA
plot(dem)

demF1 <- fill_voids_g06(stack(dem, fillSource1, fillSource2), show_progress_bar = TRUE)
plot(subset(demF1, 1))
plot(subset(demF1, 2))

plot(subset(demF1, 1) - reference)

hs <-
  hillShade(terrain(subset(demF1, 1), "slope"),
            terrain(subset(demF1, 1), "aspect"))
plot(hs , col = grey((0:255) / 255), legend = FALSE)
