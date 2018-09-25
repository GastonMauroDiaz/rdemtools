fdem <- fake_dem(r, dtm)
plot(fdem)

center <- sampleRandom(fdem, 10, xy = TRUE, sp = TRUE)
voids_mask <- fake_voids(center, 20, fdem)
plot(voids_mask)
