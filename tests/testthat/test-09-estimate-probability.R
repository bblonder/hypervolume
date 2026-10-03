# Covers proposal #4 (replace the per-test-point pdist::pdist call with one
# blocked matrix computation).
#
# The contract is an inverse-distance-weighted mean of ValueAtRandomPoints over
# the retained hypervolume points, with an optional hard zero outside the
# characteristic distance. These tests reconstruct that arithmetic directly so
# a blocked rewrite has an exact target.

prob_fixture <- function(n = 500L, seed = HV_SEED) {
  pts <- box_points(n, c(0, 0), c(1, 1), seed = seed)
  # A value field that varies across the box, so weighting actually matters.
  vals <- 0.2 + 0.6 * pts[, 1]
  methods::new(
    "Hypervolume",
    Name = "prob", Method = "test fixture",
    Data = pts, Dimensionality = 2,
    Volume = 1, PointDensity = n / 1, Parameters = list(),
    RandomPoints = pts, ValueAtRandomPoints = vals
  )
}

# Reference implementation of the per-point arithmetic, written out longhand.
ref_probability <- function(hv, points, weight.exponent = -1,
                            set.edges.zero = TRUE,
                            edges.zero.distance.factor = 1) {
  cutoff <- (nrow(hv@RandomPoints) / hv@Volume)^(-1 / ncol(hv@RandomPoints))
  vapply(seq_len(nrow(points)), function(i) {
    d <- sqrt(colSums((t(hv@RandomPoints) - points[i, ])^2))
    w <- d^weight.exponent
    res <- sum(hv@ValueAtRandomPoints * w) / sum(w)
    if (set.edges.zero && min(d) > cutoff * edges.zero.distance.factor) {
      res <- 0
    }
    res
  }, numeric(1))
}

test_that("estimate_probability matches the longhand reference", {
  hv <- prob_fixture()
  pts <- rbind(c(0.25, 0.5), c(0.5, 0.5), c(0.75, 0.5))
  colnames(pts) <- c("X1", "X2")

  observed <- quietly(hypervolume_estimate_probability(
    hv, pts, reduction.factor = 1, verbose = FALSE
  ))
  expected <- ref_probability(hv, pts)

  expect_equal(as.numeric(observed), expected, tolerance = 1e-8)
})

test_that("estimate_probability tracks the underlying value field", {
  hv <- prob_fixture()
  pts <- rbind(c(0.15, 0.5), c(0.5, 0.5), c(0.85, 0.5))
  colnames(pts) <- c("X1", "X2")

  p <- as.numeric(quietly(hypervolume_estimate_probability(hv, pts, verbose = FALSE)))

  # ValueAtRandomPoints increases with X1, so the estimates should too.
  # Inverse-distance weighting over the whole cloud shrinks everything hard
  # toward the global mean, so only assert the ordering of the extremes --
  # asserting every adjacent pair would be flaky on 500 random points.
  expect_lt(p[1], p[3])
  expect_true(all(p >= 0))
})

test_that("points far outside are zeroed when set.edges.zero is on", {
  hv <- prob_fixture()
  far <- matrix(c(50, 50), nrow = 1, dimnames = list(NULL, c("X1", "X2")))

  zeroed <- quietly(hypervolume_estimate_probability(
    hv, far, set.edges.zero = TRUE, verbose = FALSE
  ))
  kept <- quietly(hypervolume_estimate_probability(
    hv, far, set.edges.zero = FALSE, verbose = FALSE
  ))

  expect_equal(as.numeric(zeroed), 0)
  expect_gt(as.numeric(kept), 0)
})

test_that("weight.exponent changes the weighting as expected", {
  hv <- prob_fixture()
  pts <- matrix(c(0.3, 0.5), nrow = 1, dimnames = list(NULL, c("X1", "X2")))

  for (we in c(-1, -2, -0.5)) {
    observed <- quietly(hypervolume_estimate_probability(
      hv, pts, weight.exponent = we, verbose = FALSE
    ))
    expect_equal(
      as.numeric(observed),
      ref_probability(hv, pts, weight.exponent = we),
      tolerance = 1e-8,
      info = paste("weight.exponent =", we)
    )
  }
})

test_that("NA test points are dropped and recorded in the NAIdx attribute", {
  hv <- prob_fixture(n = 200L)
  pts <- rbind(c(0.3, 0.5), c(NA, 0.5), c(0.7, 0.5))
  colnames(pts) <- c("X1", "X2")

  p <- quietly(hypervolume_estimate_probability(hv, pts, verbose = FALSE))

  expect_length(as.numeric(p), 2L)
  expect_equal(as.integer(attr(p, "NAIdx")), 2L)
})

test_that("estimate_probability rejects mismatched dimensionality", {
  hv <- prob_fixture(n = 100L)
  pts <- matrix(c(0.5, 0.5, 0.5), nrow = 1)

  expect_error(
    quietly(hypervolume_estimate_probability(hv, pts, verbose = FALSE)),
    "not the same"
  )
})

test_that("estimate_probability rejects an out-of-range reduction factor", {
  hv <- prob_fixture(n = 100L)
  pts <- matrix(c(0.5, 0.5), nrow = 1)

  expect_error(
    quietly(hypervolume_estimate_probability(hv, pts, reduction.factor = 0, verbose = FALSE)),
    "not in"
  )
  expect_error(
    quietly(hypervolume_estimate_probability(hv, pts, reduction.factor = 2, verbose = FALSE)),
    "not in"
  )
})

test_that("inclusion test fast path agrees with a ball query", {
  hv <- prob_fixture(n = 800L)

  # Query at points the hypervolume already contains, so the ball query finds
  # them at distance 0 regardless of how the random cloud fell. Picking
  # arbitrary interior coordinates instead would be flaky: the characteristic
  # distance is ~n^(-1/2), so an interior point has a real chance of having no
  # neighbour inside it.
  inside_pts <- hv@RandomPoints[c(1L, 137L, 500L), , drop = FALSE]
  far <- matrix(c(5, 5), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  pts <- rbind(inside_pts, far)

  inside <- suppressWarnings(quietly(hypervolume_inclusion_test(
    hv, pts, fast.or.accurate = "fast", verbose = FALSE
  )))

  expect_type(as.vector(inside), "logical")
  expect_length(as.vector(inside), 4L)
  expect_true(all(inside[1:3]))
  expect_false(inside[4])
})

test_that("inclusion test accurate path returns logicals", {
  hv <- prob_fixture(n = 400L)

  # Two hazards on the accurate path (which forwards to
  # hypervolume_estimate_probability with set.edges.zero = TRUE):
  # 1. Characteristic distance is only 400^(-1/2) = 0.05, so an arbitrary
  #    interior point has a ~4% chance of having no neighbour that close,
  #    which would hard-zero the estimate. Query a point the cloud contains.
  # 2. Querying an exact fixture point makes one pdist distance zero;
  #    weight = distance^(-1) is then Inf and sum(prob*Inf)/sum(Inf) is NaN
  #    -- see the "known-defects" section of TESTING.md. Offset by an
  #    epsilon large enough to survive pdist's ~1e-6 catastrophic
  #    cancellation floor but well below cutoff_dist (=0.05 here).
  inside <- hv@RandomPoints[1L, , drop = FALSE] + 1e-4
  far <- matrix(c(20, 20), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  pts <- rbind(inside, far)

  res <- quietly(hypervolume_inclusion_test(
    hv, pts, fast.or.accurate = "accurate",
    accurate.method.threshold = 0.1, verbose = FALSE
  ))

  expect_type(as.vector(res), "logical")
  expect_true(res[1])
  expect_false(res[2])
})

test_that("inclusion test rejects an unknown method", {
  hv <- prob_fixture(n = 100L)
  pts <- matrix(c(0.5, 0.5), nrow = 1, dimnames = list(NULL, c("X1", "X2")))

  expect_error(
    quietly(hypervolume_inclusion_test(hv, pts, fast.or.accurate = "nonsense", verbose = FALSE)),
    "Unsupported"
  )
})
