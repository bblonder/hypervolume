# Covers proposal #7 (hoist the kd-tree build out of the sampler's while loop,
# and stop round-tripping through data.frame).
#
# The safety argument for hoisting is that ellipsoid_inverse_weight is a pure
# function of (samples, centers, scales) and `centers` never changes inside the
# loop. These tests assert that purity, then pin the volume accounting on cases
# where the answer is exact rather than Monte Carlo.

test_that("nball_volume matches known values", {
  expect_equal(hypervolume:::nball_volume(1, 1), 2, tolerance = 1e-12)
  expect_equal(hypervolume:::nball_volume(2, 1), pi, tolerance = 1e-12)
  expect_equal(hypervolume:::nball_volume(3, 1), 4 / 3 * pi, tolerance = 1e-12)
  # radius scaling
  expect_equal(
    hypervolume:::nball_volume(3, 2),
    8 * hypervolume:::nball_volume(3, 1),
    tolerance = 1e-12
  )
})

test_that("ellipsoid_volume is the unit ball scaled by prod(scales)", {
  scales <- c(2, 3, 0.5)
  expect_equal(
    hypervolume:::ellipsoid_volume(scales),
    hypervolume:::nball_volume(3, 1) * prod(scales),
    tolerance = 1e-12
  )
})

test_that("sample_ellipsoid draws inside the ellipsoid and centres correctly", {
  centre <- c(10, -4)
  scales <- c(2, 0.5)
  pts <- with_seed(HV_SEED, hypervolume:::sample_ellipsoid(centre, 5000L, scales))

  expect_equal(dim(pts), c(5000L, 2L))

  # Every draw must satisfy sum(((x - c)/s)^2) <= 1.
  normed <- sweep(sweep(pts, 2, centre, "-"), 2, scales, "/")
  expect_true(all(rowSums(normed^2) <= 1 + 1e-9))

  # And the cloud should be centred on `centre`.
  expect_equal(colMeans(pts), centre, tolerance = 0.05)
})

test_that("ellipsoid_inverse_weight counts containing ellipsoids", {
  centres <- rbind(c(0, 0), c(10, 10))
  colnames(centres) <- c("X1", "X2")
  scales <- c(1, 1)

  samples <- rbind(
    c(0, 0),      # inside centre 1 only
    c(10, 10),    # inside centre 2 only
    c(0.5, 0.5),  # inside centre 1 only
    c(50, 50)     # inside neither
  )
  colnames(samples) <- c("X1", "X2")

  got <- quietly(hypervolume:::ellipsoid_inverse_weight(
    samples, centers = centres, scales = scales, verbose = FALSE
  ))

  expect_equal(as.integer(got), c(1L, 1L, 1L, 0L))
})

test_that("ellipsoid_inverse_weight counts overlaps", {
  # Two centres 1 apart with scale 1: the midpoint sits in both unit balls.
  centres <- rbind(c(0, 0), c(1, 0))
  colnames(centres) <- c("X1", "X2")
  samples <- rbind(c(0.5, 0), c(-0.9, 0), c(1.9, 0))
  colnames(samples) <- c("X1", "X2")

  got <- quietly(hypervolume:::ellipsoid_inverse_weight(
    samples, centers = centres, scales = c(1, 1), verbose = FALSE
  ))

  expect_equal(as.integer(got), c(2L, 1L, 1L))
})

test_that("ellipsoid_inverse_weight is pure and repeatable", {
  # This is the invariant that makes proposal #7 safe: the result depends only
  # on its arguments, so the tree over `centers` can be built once and reused
  # across iterations of the sampler's while loop.
  centres <- toy_data(30L, 2L, seed = HV_SEED)
  samples <- toy_data(120L, 2L, seed = HV_SEED + 21L)
  scales <- c(0.8, 1.2)

  a <- quietly(hypervolume:::ellipsoid_inverse_weight(samples, centres, scales, FALSE))
  b <- quietly(hypervolume:::ellipsoid_inverse_weight(samples, centres, scales, FALSE))
  c3 <- quietly(hypervolume:::ellipsoid_inverse_weight(samples, centres, scales, FALSE))

  expect_equal(as.integer(a), as.integer(b))
  expect_equal(as.integer(a), as.integer(c3))
})

test_that("ellipsoid_inverse_weight respects anisotropic scales", {
  centres <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  samples <- rbind(c(1.5, 0), c(0, 1.5))
  colnames(samples) <- c("X1", "X2")

  # Wide in X1, narrow in X2: the first sample is in, the second is out.
  got <- quietly(hypervolume:::ellipsoid_inverse_weight(
    samples, centers = centres, scales = c(2, 1), verbose = FALSE
  ))
  expect_equal(as.integer(got), c(1L, 0L))
})

# ---- sample_model_ellipsoid volume accounting -------------------------------

accept_everything <- function(x) rep(1, nrow(x))

test_that("sampled volume of one ellipsoid is exactly the ellipsoid volume", {
  # With a single centre, every draw lies in exactly one ellipsoid, so the
  # inverse weights are all 1 and the Monte Carlo correction is exactly 1.
  # The reported volume is therefore analytic, not noisy.
  data <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  scales <- c(1, 1)

  res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
    predict_function = accept_everything,
    data = data, scales = scales, min.value = 0,
    samples.per.point = 100, chunk.size = 500, verbose = FALSE
  )))

  expect_equal(res$volume, hypervolume:::ellipsoid_volume(scales), tolerance = 1e-8)
  expect_equal(ncol(res$samples), 3L)  # d columns plus "value"
  expect_equal(nrow(res$samples), 100L)
})

test_that("two coincident centres do not double-count volume", {
  # Both ellipsoids are identical, so every sample has inverse weight 2 and the
  # union volume must still be one ellipsoid. This is the overlap correction
  # working, and it is the part of the sampler proposal #7 must not disturb.
  data <- rbind(c(0, 0), c(0, 0))
  colnames(data) <- c("X1", "X2")
  scales <- c(1, 1)

  res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
    predict_function = accept_everything,
    data = data, scales = scales, min.value = 0,
    samples.per.point = 100, chunk.size = 500, verbose = FALSE
  )))

  expect_equal(res$volume, hypervolume:::ellipsoid_volume(scales), tolerance = 1e-8)
})

test_that("two far-apart centres give twice the ellipsoid volume", {
  data <- rbind(c(0, 0), c(100, 100))
  colnames(data) <- c("X1", "X2")
  scales <- c(1, 1)

  res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
    predict_function = accept_everything,
    data = data, scales = scales, min.value = 0,
    samples.per.point = 100, chunk.size = 500, verbose = FALSE
  )))

  expect_equal(res$volume, 2 * hypervolume:::ellipsoid_volume(scales), tolerance = 1e-8)
})

test_that("rejecting half the samples roughly halves the volume", {
  # Accept only points with X1 > 0. The ellipsoid is symmetric about 0, so the
  # retained fraction should be ~1/2 and the volume ~half the ellipsoid.
  data <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  scales <- c(1, 1)
  half_plane <- function(x) as.numeric(x[, 1] > 0)

  res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
    predict_function = half_plane,
    data = data, scales = scales, min.value = 0,
    samples.per.point = 400, chunk.size = 2000, verbose = FALSE
  )))

  expect_volume_close(res$volume, hypervolume:::ellipsoid_volume(scales) / 2, tol = 0.10)
  expect_true(all(res$samples[, "X1"] > 0))
})

test_that("sample_model_ellipsoid drops incomplete cases and names its columns", {
  data <- rbind(c(0, 0), c(NA, 1), c(0.2, 0.1))
  colnames(data) <- c("X1", "X2")

  res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
    predict_function = accept_everything,
    data = data, scales = c(1, 1), min.value = 0,
    samples.per.point = 50, chunk.size = 400, verbose = FALSE
  )))

  expect_equal(colnames(res$samples), c("X1", "X2", "value"))
  expect_false(anyNA(res$samples))
})

test_that("sample_model_rejection honours its range and acceptance threshold", {
  # The other sampler mentioned in proposal #7 (data.frame round trip in the
  # rescaling loop). Uses lm() purely because it comes with a predict method.
  # X2 must not be collinear with X1, or lm() is rank-deficient and every
  # predict() call warns.
  df <- data.frame(
    X1 = seq(0, 1, length.out = 40),
    X2 = seq(0, 1, length.out = 40)^2
  )
  df$y <- df$X1
  model <- stats::lm(y ~ X1 + X2, data = df)

  rng <- rbind(c(0, 0), c(1, 1))
  dimnames(rng) <- list(NULL, c("X1", "X2"))

  samples <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_rejection(
    model = model, range = rng, N.samples = 200,
    chunk.size = 500, verbose = FALSE, min.value = 0.5
  )))

  expect_true(is.matrix(samples))
  expect_equal(colnames(samples), c("X1", "X2", "value"))
  # every draw must sit inside the requested range
  expect_true(all(samples[, "X1"] >= 0 & samples[, "X1"] <= 1))
  expect_true(all(samples[, "X2"] >= 0 & samples[, "X2"] <= 1))
  # and the cutoff must leave exactly N.samples above the threshold
  expect_equal(sum(samples[, "value"] > 0.5), 200L)
})

test_that("sample_model_ellipsoid returns the requested number of samples", {
  data <- toy_data(5L, 2L)
  for (spp in c(20, 75, 150)) {
    res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
      predict_function = accept_everything,
      data = data, scales = c(1, 1), min.value = 0,
      samples.per.point = spp, chunk.size = 400, verbose = FALSE
    )))
    expect_equal(nrow(res$samples), ceiling(spp * nrow(data)), info = paste("spp =", spp))
  }
})
