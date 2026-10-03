# End-to-end cover for proposal #2 and #3 acting together: hypervolume_gaussian
# is the caller that drives calculate_density, sample_model_ellipsoid and
# hypervolume_threshold(num.thresholds = 1000).
#
# Kept deliberately small (few data points, low samples.per.point) so the whole
# file runs in seconds. The heavier checks are behind skip_if_not_slow().

fixed_bw <- function(data, value) {
  estimate_bandwidth(data, method = "fixed", value = value)
}

small_gaussian_hv <- function(n = 12L, d = 2L, value = 0.6,
                              samples.per.point = 40, seed = HV_SEED) {
  data <- toy_data(n, d, seed = seed)
  with_seed(seed, quietly(hypervolume_gaussian(
    data,
    name = "test",
    kde.bandwidth = fixed_bw(data, value),
    samples.per.point = samples.per.point,
    chunk.size = 300,
    verbose = FALSE
  )))
}

test_that("hypervolume_gaussian returns a well-formed Hypervolume", {
  hv <- small_gaussian_hv()

  expect_valid_hypervolume(hv)
  expect_equal(hv@Dimensionality, 2)
  expect_equal(hv@Method, "Gaussian kernel density estimate")
  expect_gt(hv@Volume, 0)
  expect_true(all(hv@ValueAtRandomPoints > 0))
  expect_equal(colnames(hv@RandomPoints), c("X1", "X2"))
})

test_that("hypervolume_gaussian records its parameters", {
  hv <- small_gaussian_hv()

  expect_true(all(c(
    "kde.bandwidth", "kde.method", "samples.per.point",
    "sd.count", "quantile.requested", "quantile.requested.type"
  ) %in% names(hv@Parameters)))
  expect_equal(hv@Parameters$kde.method, "fixed")
  expect_equal(hv@Parameters$quantile.requested.type, "probability")
})

test_that("hypervolume_gaussian is reproducible under a fixed seed", {
  a <- small_gaussian_hv(seed = HV_SEED)
  b <- small_gaussian_hv(seed = HV_SEED)

  expect_equal(a@Volume, b@Volume, tolerance = 1e-12)
  expect_equal(nrow(a@RandomPoints), nrow(b@RandomPoints))
  expect_equal(a@ValueAtRandomPoints, b@ValueAtRandomPoints, tolerance = 1e-12)
})

test_that("hypervolume_gaussian volume scales as c^d under rescaling", {
  # THE KEY INVARIANT FOR PROPOSAL #2. Multiplying the data and the bandwidth by
  # c maps every internal quantity by a known power of c: sampled points by c,
  # densities by c^-d, volume by c^d. Because acceptance is deterministic
  # (min.value = 0 and the density is strictly positive) and thresholding is
  # quantile-based, the same RNG stream picks the same points in both runs, so
  # this should hold to floating-point precision, not just Monte Carlo
  # precision. Any rewrite of calculate_density that gets the normalising
  # constant wrong will break it.
  d <- 2L
  base_data <- toy_data(10L, d, seed = HV_SEED)
  cc <- 3

  hv1 <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    base_data, kde.bandwidth = fixed_bw(base_data, 0.6),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))
  scaled <- base_data * cc
  hv2 <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    scaled, kde.bandwidth = fixed_bw(scaled, 0.6 * cc),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))

  expect_equal(hv2@Volume / hv1@Volume, cc^d, tolerance = 1e-6)
  expect_equal(nrow(hv2@RandomPoints), nrow(hv1@RandomPoints))
})

test_that("hypervolume_gaussian volume is translation invariant", {
  base_data <- toy_data(10L, 2L, seed = HV_SEED)
  shift <- c(100, -50)

  hv1 <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    base_data, kde.bandwidth = fixed_bw(base_data, 0.6),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))
  moved <- sweep(base_data, 2, shift, "+")
  hv2 <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    moved, kde.bandwidth = fixed_bw(moved, 0.6),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))

  expect_equal(hv2@Volume, hv1@Volume, tolerance = 1e-6)
  expect_equal(
    unname(get_centroid(hv2) - get_centroid(hv1)),
    unname(shift),
    tolerance = 0.2
  )
})

test_that("a wider bandwidth yields a larger volume", {
  data <- toy_data(10L, 2L, seed = HV_SEED)

  narrow <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    data, kde.bandwidth = fixed_bw(data, 0.4),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))
  wide <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    data, kde.bandwidth = fixed_bw(data, 0.9),
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))

  expect_gt(wide@Volume, narrow@Volume)
})

test_that("uniform weights and default weights agree", {
  # hypervolume_gaussian fills weight with rep(1/n, n) when it is NULL, so
  # passing that explicitly must be a no-op. Proposal #2 changes how weight
  # enters the sum, so this pins the equivalence.
  data <- toy_data(10L, 2L, seed = HV_SEED)
  bw <- fixed_bw(data, 0.6)

  default_w <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    data, kde.bandwidth = bw, samples.per.point = 40,
    chunk.size = 300, verbose = FALSE
  )))
  explicit_w <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    data, weight = rep(1 / nrow(data), nrow(data)), kde.bandwidth = bw,
    samples.per.point = 40, chunk.size = 300, verbose = FALSE
  )))

  expect_equal(explicit_w@Volume, default_w@Volume, tolerance = 1e-10)
  expect_equal(
    explicit_w@ValueAtRandomPoints,
    default_w@ValueAtRandomPoints,
    tolerance = 1e-10
  )
})

test_that("unnormalised weights are renormalised with a warning", {
  data <- toy_data(8L, 2L, seed = HV_SEED)
  bw <- fixed_bw(data, 0.6)

  expect_warning(
    with_seed(HV_SEED, utils::capture.output(suppressMessages(
      hypervolume_gaussian(
        data, weight = rep(1, nrow(data)), kde.bandwidth = bw,
        samples.per.point = 20, chunk.size = 200, verbose = FALSE
      )
    ))),
    "sum of the weights"
  )
})

test_that("hypervolume_gaussian validates its inputs", {
  data <- toy_data(8L, 2L, seed = HV_SEED)
  bw <- fixed_bw(data, 0.6)

  # bandwidth of the wrong length
  expect_error(
    quietly(hypervolume_gaussian(
      data, kde.bandwidth = fixed_bw(data, 0.6)[1], verbose = FALSE
    )),
    "same dimensionality"
  )

  # zero bandwidth
  expect_error(
    quietly(hypervolume_gaussian(
      data, kde.bandwidth = estimate_bandwidth(data, method = "fixed", value = 0),
      verbose = FALSE
    )),
    "non-zero"
  )

  # bandwidth without the method attribute set by estimate_bandwidth()
  raw_bw <- c(0.6, 0.6)
  expect_error(
    quietly(hypervolume_gaussian(data, kde.bandwidth = raw_bw, verbose = FALSE)),
    "estimate_bandwidth"
  )

  # weight of the wrong length
  expect_error(
    quietly(hypervolume_gaussian(
      data, weight = rep(1 / 3, 3), kde.bandwidth = bw, verbose = FALSE
    )),
    "length of the weights"
  )
})

test_that("a low sd.count warns", {
  data <- toy_data(8L, 2L, seed = HV_SEED)
  bw <- fixed_bw(data, 0.6)

  expect_warning(
    with_seed(HV_SEED, utils::capture.output(suppressMessages(
      hypervolume_gaussian(
        data, kde.bandwidth = bw, sd.count = 2,
        samples.per.point = 20, chunk.size = 200, verbose = FALSE
      )
    ))),
    "sd.count"
  )
})

test_that("hypervolume_gaussian works in 1 and 3 dimensions", {
  for (d in c(1L, 3L)) {
    data <- toy_data(10L, d, seed = HV_SEED + d)
    hv <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
      data, kde.bandwidth = fixed_bw(data, 0.6),
      samples.per.point = 30, chunk.size = 300, verbose = FALSE
    )))
    expect_valid_hypervolume(hv)
    expect_equal(hv@Dimensionality, d, info = paste("d =", d))
  }
})

test_that("get_volume and get_centroid agree with the slots", {
  hv <- small_gaussian_hv()

  expect_equal(unname(get_volume(hv)), hv@Volume)
  expect_equal(
    unname(get_centroid(hv)),
    unname(colMeans(hv@RandomPoints)),
    tolerance = 1e-10
  )
})

test_that("a recognisable 2D Gaussian recovers a plausible volume", {
  # Loose sanity check against theory rather than a tight assertion: for a
  # standard bivariate normal the 95%-probability region has area
  # -2*pi*log(0.05) ~= 18.8. With a small sample and a fixed bandwidth we only
  # ask for the right order of magnitude.
  skip_if_not_slow()

  data <- with_seed(HV_SEED, matrix(stats::rnorm(400), ncol = 2))
  colnames(data) <- c("X1", "X2")

  hv <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
    data, kde.bandwidth = estimate_bandwidth(data, method = "silverman-1d"),
    samples.per.point = 200, verbose = FALSE
  )))

  expect_gt(hv@Volume, 5)
  expect_lt(hv@Volume, 60)
})
