set.seed(11)
reference <- fake_dem(n_random_data = 60, z_range = c(0, 600))
plot(reference)

set.seed(18)
fillSource <- fake_dem(z_range = c(0, 20))
fillSource <- fillSource - mean(fillSource[])
fillSource <- fillSource + reference
fillSource <- smooth_dem(fillSource, theta = 6)
plot(fillSource)

p <- sampleRegular(reference, 10, sp = TRUE)
voidsMask <- fake_voids(p, void_size = 20, reference)

dem <- reference
dem[voidsMask] <- NA
plot(dem)

slp <- terrain(fillSource, "slope")
asp <- terrain(fillSource, "aspect")
demF2 <- fill_voids_k08(dem, slp, asp, show_progress_bar = TRUE)

plot(demF2)

plot(demF2 - reference)

hs <- hillShade(terrain(demF2, "slope"), terrain(demF2, "aspect"))
plot(hs , col= grey((0:255)/255), legend = FALSE)
