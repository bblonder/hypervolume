# Covers proposal #3 (replace the 1000-iteration threshold loop with sorted
# cumulative sums) and the volume_calculation helper it sits next to.
#
# The rewrite must reproduce the Statistics data frame exactly and return the
# same selected hypervolume. Everything below is phrased against that contract.

test_that("volume_calculation reproduces a hand-worked example", {
  counts <- c(1L, 1L, 2L, 2L, 2L, 3L)
  # tabulate -> [2, 3, 1] for values 1, 2, 3
  # weighted  -> 2/1, 3/2, 1/3
  weighted <- c(2 / 1, 3 / 2, 1 / 3)
  cumulated <- cumsum(weighted) / sum(weighted)

  res <- quietly(hypervolume:::volume_calculation(
    counts, point_density = 2, quantile = 0.5, verbose = FALSE
  ))

  qindex <- head(which(cumulated >= 0.5), 1)
  expect_equal(res$index_out, qindex)
  expect_equal(res$final_volume, sum(weighted[qindex:length(weighted)]) / 2)
})

test_that("volume_calculation scales inversely with point density", {
  counts <- c(1L, 2L, 2L, 3L, 3L, 3L)

  a <- quietly(hypervolume:::volume_calculation(counts, 1, 0.5, verbose = FALSE))
  b <- quietly(hypervolume:::volume_calculation(counts, 4, 0.5, verbose = FALSE))

  expect_equal(a$final_volume / b$final_volume, 4, tolerance = 1e-10)
})

test_that("volume_calculation is monotone decreasing in the requested quantile", {
  counts <- rep(1:6, times = 6:1)
  vols <- vapply(c(0.1, 0.3, 0.5, 0.7, 0.9), function(q) {
    quietly(hypervolume:::volume_calculation(counts, 1, q, verbose = FALSE))$final_volume
  }, numeric(1))

  expect_true(all(diff(vols) <= 1e-12))
})

# ---- hypervolume_threshold --------------------------------------------------

make_threshold_hv <- function(n = 400L, seed = HV_SEED) {
  pts <- box_points(n, c(0, 0), c(1, 1), seed = seed)
  vals <- with_seed(seed + 1L, stats::runif(n))
  methods::new(
    "Hypervolume",
    Name = "thresh", Method = "test fixture",
    Data = pts, Dimensionality = 2,
    Volume = 1, PointDensity = n / 1, Parameters = list(),
    RandomPoints = pts, ValueAtRandomPoints = vals
  )
}

test_that("hypervolume_threshold returns one hypervolume per threshold", {
  hv <- make_threshold_hv()
  res <- quietly(hypervolume_threshold(
    hv, num.thresholds = 12, plot = FALSE, verbose = FALSE
  ))

  expect_named(res, c("HypervolumesThresholded", "Statistics"))
  expect_s4_class(res$HypervolumesThresholded, "HypervolumeList")
  expect_length(res$HypervolumesThresholded@HVList, 12L)
  expect_equal(nrow(res$Statistics), 12L)
  expect_named(
    res$Statistics,
    c("threshold", "volume", "volume_quantile", "probability_quantile")
  )
})

test_that("threshold statistics equal a direct cumulative-sum calculation", {
  # This is exactly the identity proposal #3 relies on: every row of
  # Statistics is derivable from the sorted ValueAtRandomPoints without
  # materialising a single sub-hypervolume.
  hv <- make_threshold_hv(n = 300L)
  res <- quietly(hypervolume_threshold(
    hv, num.thresholds = 25, plot = FALSE, verbose = FALSE
  ))

  v <- hv@ValueAtRandomPoints
  thresholds <- res$Statistics$threshold

  expected_volume <- vapply(thresholds, function(t) sum(v > t) / hv@PointDensity, numeric(1))
  expected_psum   <- vapply(thresholds, function(t) sum(v[v > t]), numeric(1))

  expect_equal(res$Statistics$volume, expected_volume, tolerance = 1e-12)
  expect_equal(
    res$Statistics$probability_quantile,
    expected_psum / max(expected_psum),
    tolerance = 1e-12
  )
  expect_equal(
    res$Statistics$volume_quantile,
    expected_volume / max(expected_volume),
    tolerance = 1e-12
  )
})

test_that("threshold volumes and quantiles are monotone non-increasing", {
  hv <- make_threshold_hv()
  res <- quietly(hypervolume_threshold(
    hv, num.thresholds = 30, plot = FALSE, verbose = FALSE
  ))

  expect_true(all(diff(res$Statistics$volume) <= 1e-12))
  expect_true(all(diff(res$Statistics$volume_quantile) <= 1e-12))
  expect_true(all(diff(res$Statistics$probability_quantile) <= 1e-12))
  expect_lte(max(res$Statistics$volume_quantile), 1 + 1e-12)
})

test_that("each thresholded hypervolume holds exactly the points above its threshold", {
  hv <- make_threshold_hv(n = 200L)
  res <- quietly(hypervolume_threshold(
    hv, num.thresholds = 8, plot = FALSE, verbose = FALSE
  ))

  for (i in seq_along(res$Statistics$threshold)) {
    t_i <- res$Statistics$threshold[i]
    keep <- hv@ValueAtRandomPoints > t_i
    got <- res$HypervolumesThresholded@HVList[[i]]

    expect_equal(nrow(got@RandomPoints), sum(keep), info = paste("i =", i))
    expect_equal(
      unname(got@RandomPoints),
      unname(hv@RandomPoints[keep, , drop = FALSE]),
      info = paste("i =", i)
    )
  }
})

test_that("uniform.density controls whether values are flattened to 1", {
  hv <- make_threshold_hv(n = 150L)

  flat <- quietly(hypervolume_threshold(
    hv, num.thresholds = 5, uniform.density = TRUE, plot = FALSE, verbose = FALSE
  ))
  kept <- quietly(hypervolume_threshold(
    hv, num.thresholds = 5, uniform.density = FALSE, plot = FALSE, verbose = FALSE
  ))

  flat2 <- flat$HypervolumesThresholded@HVList[[2]]
  expect_true(all(flat2@ValueAtRandomPoints == 1))

  # PINS A KNOWN DEFECT. hypervolume_threshold.R does
  #     hv_new@ValueAtRandomPoints <- rep(1, length(hv_new@ValueAtRandomPoints))
  # but at that line hv_new is still a copy of hv with only @RandomPoints
  # subset, so length() is the ORIGINAL n, not the retained count. Every
  # uniform.density = TRUE result therefore carries more values than points.
  #
  # Asserting the buggy length rather than the correct one, so the suite stays
  # green and proposal #3 cannot change it by accident. The correct assertion
  # is nrow(flat2@RandomPoints); swap to that when the slot is fixed.
  expect_length(flat2@ValueAtRandomPoints, nrow(hv@RandomPoints))
  expect_gt(length(flat2@ValueAtRandomPoints), nrow(flat2@RandomPoints))

  t2 <- kept$Statistics$threshold[2]
  expect_equal(
    kept$HypervolumesThresholded@HVList[[2]]@ValueAtRandomPoints,
    hv@ValueAtRandomPoints[hv@ValueAtRandomPoints > t2],
    tolerance = 1e-12
  )
})

test_that("requesting a quantile returns a single hypervolume at the right index", {
  # This is the path hypervolume_gaussian uses (num.thresholds = 1000,
  # quantile.requested.type = "probability"), and the only one whose output
  # survives. Proposal #3 must keep the selected index identical.
  hv <- make_threshold_hv(n = 300L)

  for (type in c("volume", "probability")) {
    res <- quietly(hypervolume_threshold(
      hv, num.thresholds = 50, quantile.requested = 0.9,
      quantile.requested.type = type, plot = FALSE, verbose = FALSE
    ))

    expect_s4_class(res$HypervolumesThresholded, "Hypervolume")

    col <- if (type == "probability") "probability_quantile" else "volume_quantile"
    idx <- head(which(res$Statistics[[col]] <= 0.9), 1)
    t_idx <- res$Statistics$threshold[idx]

    expect_equal(
      nrow(res$HypervolumesThresholded@RandomPoints),
      sum(hv@ValueAtRandomPoints > t_idx),
      info = type
    )
  }
})

test_that("hypervolume_threshold rejects an unknown quantile type", {
  hv <- make_threshold_hv(n = 50L)
  expect_error(
    quietly(hypervolume_threshold(
      hv, num.thresholds = 5, quantile.requested = 0.5,
      quantile.requested.type = "nonsense", plot = FALSE, verbose = FALSE
    )),
    "Unrecognized"
  )
})

test_that("explicit thresholds are used verbatim", {
  hv <- make_threshold_hv(n = 120L)
  tvals <- c(0.1, 0.4, 0.8)

  res <- quietly(hypervolume_threshold(
    hv, thresholds = tvals, plot = FALSE, verbose = FALSE
  ))

  expect_equal(res$Statistics$threshold, tvals)
})
