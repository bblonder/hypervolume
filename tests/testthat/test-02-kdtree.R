# Covers proposal #1 (the O(n^2) workarray clear in build_recursively) and
# proposal #5 (flat point storage / bucketed leaves) and proposal #6 (OpenMP).
#
# All three change how the tree is built or traversed, so what matters is that
# query RESULTS are invariant. These tests compare against brute force on
# inputs chosen to stress the split logic: duplicates, collinear points, one
# point, two points (the heapsort len==2 bug noted in the source comments),
# and sizes that straddle the `m log m < n` branch switch.

test_that("kdtree_build accepts matrices and data frames alike", {
  d <- toy_data(20L, 3L)

  t1 <- quietly(hypervolume:::kdtree_build(d, verbose = FALSE))
  t2 <- quietly(hypervolume:::kdtree_build(as.data.frame(d), verbose = FALSE))

  # Rcpp XPtr comes back as a base externalptr, not an S4 object.
  expect_type(t1, "externalptr")
  expect_type(t2, "externalptr")
})

test_that("kdtree_build rejects non-tabular input", {
  expect_error(
    hypervolume:::kdtree_build(1:10, verbose = FALSE),
    "not a matrix or data frame"
  )
})

test_that("ball query counts match brute force", {
  data <- toy_data(60L, 3L, seed = HV_SEED)
  pts  <- toy_data(25L, 3L, seed = HV_SEED + 7L)

  for (radius in c(0.25, 0.75, 1.5, 4.0)) {
    observed <- quietly(hypervolume:::evalfspherical(
      data = data, radius = radius, points = pts, verbose = FALSE
    ))
    expected <- ref_ball_counts(data, pts, radius)
    expect_equal(as.integer(observed), expected, info = paste("radius =", radius))
  }
})

test_that("ball query counts match brute force across dimensionalities", {
  for (d in c(1L, 2L, 4L, 6L)) {
    data <- toy_data(40L, d, seed = HV_SEED + d)
    pts  <- toy_data(15L, d, seed = HV_SEED + 500L + d)
    radius <- 1.25

    observed <- quietly(hypervolume:::evalfspherical(
      data = data, radius = radius, points = pts, verbose = FALSE
    ))
    expect_equal(
      as.integer(observed),
      ref_ball_counts(data, pts, radius),
      info = paste("d =", d)
    )
  }
})

test_that("ball query is correct at the branch switch in build_recursively", {
  # build_recursively picks heapsort vs. linear scan on `m log m < npoints`.
  # These sizes put nodes on both sides of that boundary.
  for (n in c(2L, 3L, 4L, 8L, 17L, 33L, 64L, 129L)) {
    data <- toy_data(n, 2L, seed = HV_SEED + n)
    pts  <- toy_data(6L, 2L, seed = HV_SEED + 900L + n)
    radius <- 1.0

    observed <- quietly(hypervolume:::evalfspherical(
      data = data, radius = radius, points = pts, verbose = FALSE
    ))
    expect_equal(
      as.integer(observed),
      ref_ball_counts(data, pts, radius),
      info = paste("n =", n)
    )
  }
})

test_that("ball query handles duplicate and collinear points", {
  # Duplicates force zero-width splits; collinear points make every split in
  # one dimension degenerate.
  dup <- matrix(rep(c(0, 0), each = 12), ncol = 2)
  colnames(dup) <- c("X1", "X2")
  collinear <- cbind(X1 = seq(0, 1, length.out = 16), X2 = 0)
  data <- rbind(dup, collinear)

  pts <- rbind(c(0, 0), c(0.5, 0), c(2, 2))
  colnames(pts) <- c("X1", "X2")

  for (radius in c(0.1, 0.6, 3.0)) {
    observed <- quietly(hypervolume:::evalfspherical(
      data = data, radius = radius, points = pts, verbose = FALSE
    ))
    expect_equal(
      as.integer(observed),
      ref_ball_counts(data, pts, radius),
      info = paste("radius =", radius)
    )
  }
})

test_that("ball query returns zero when nothing is in range", {
  data <- toy_data(20L, 2L)
  far  <- matrix(c(1e6, 1e6), nrow = 1, dimnames = list(NULL, c("X1", "X2")))

  observed <- quietly(hypervolume:::evalfspherical(
    data = data, radius = 0.5, points = far, verbose = FALSE
  ))
  expect_equal(as.integer(observed), 0L)
})

test_that("ball query id variant returns the same membership as the counts", {
  data <- toy_data(40L, 2L, seed = HV_SEED)
  pts  <- toy_data(10L, 2L, seed = HV_SEED + 3L)
  radius <- 1.0

  res <- quietly(hypervolume:::evalfspherical(
    data = data, radius = radius, points = pts,
    getid.nearestneighbor = TRUE, verbose = FALSE
  ))

  expect_type(res, "list")
  expect_length(res, 2L)

  counts <- as.integer(res[[1]])
  ids    <- res[[2]]
  expected_ids <- ref_ball_ids(data, pts, radius)

  expect_equal(counts, ref_ball_counts(data, pts, radius))
  expect_length(ids, nrow(pts))

  for (i in seq_along(ids)) {
    this <- ids[[i]]
    if (identical(as.integer(this), -1L)) {
      # documented sentinel for "nothing in range"
      expect_equal(length(expected_ids[[i]]), 0L, info = paste("i =", i))
    } else {
      # C++ returns 0-based indices
      expect_setequal(sort(as.integer(this) + 1L), expected_ids[[i]])
    }
  }
})

test_that("dimension mismatch between tree and query points is an error", {
  data <- toy_data(15L, 3L)
  tree <- quietly(hypervolume:::kdtree_build(data, verbose = FALSE))
  pts  <- toy_data(5L, 2L)

  expect_error(
    hypervolume:::kdtree_ball_query_multiple(
      tree, t(pts), nrow(pts), ncol(pts), 1.0, FALSE
    ),
    "dimensionality"
  )
})

test_that("range query matches brute force", {
  data <- toy_data(50L, 3L, seed = HV_SEED)
  pts  <- toy_data(12L, 3L, seed = HV_SEED + 11L)
  bandwidth <- c(0.4, 0.8, 1.2)

  observed <- quietly(hypervolume:::evalfrectangular(
    data = data, bandwidth = bandwidth, points = pts, verbose = FALSE
  ))
  expect_equal(
    as.integer(observed),
    ref_box_counts(data, pts, bandwidth),
    tolerance = 0
  )
})

test_that("repmat tiles as the rectangular evaluator expects", {
  X <- matrix(1:6, nrow = 2)
  out <- hypervolume:::repmat(X, 3, 2)

  expect_equal(dim(out), c(6L, 6L))
  expect_equal(out[1:2, 1:3], X)
  expect_equal(out[3:4, 1:3], X)
  expect_equal(out[1:2, 4:6], X)
})

test_that("tree build cost is not quadratic in n", {
  # Acceptance test for proposal #1 (KDTree.h workarray clear).
  # 4x the points at O(n log n) costs ~5x the time; at O(n^2) it costs
  # ~16x. Assert the ratio stays below 8, which fails on the unfixed
  # source (~15x at these sizes) and passes with the fix (~5x).
  #
  # Uses microbenchmark rather than a single proc.time() reading so we
  # can keep n small (fixed builds at n=5k are ~3ms -- individual reads
  # are noisy at that scale, but the median of 15 microbenchmark
  # iterations is stable to a few percent). Slow-gated because 20k-point
  # builds unfixed are hundreds of ms even in the median.
  skip_if_not_slow()
  skip_on_cran()
  skip_if_not_installed("microbenchmark")

  bench_build <- function(n) {
    d <- toy_data(n, 3L, seed = HV_SEED + n)
    mb <- microbenchmark::microbenchmark(
      quietly(hypervolume:::kdtree_build(d, verbose = FALSE)),
      times = 15L
    )
    median(mb$time) / 1e6  # ms
  }

  small <- bench_build(5000L)
  large <- bench_build(20000L)

  expect_lt(large / small, 8)
})

test_that("ball-query results are independent of OpenMP thread count", {
  # Acceptance test for proposal #6 (OpenMP over kdtree_ball_query_multiple
  # and friends). Skipped until proposal #6 lands and the build sets
  # SHLIB_OPENMP_CXXFLAGS -- see skip_if_not_omp() in helper-hv.R.
  #
  # What this guards: proposal #6's #pragma omp parallel for over query points
  # must write into a pre-sized result vector, not push_back into an
  # unreserved one. A race in the count/id append path would show up here as
  # counts drifting between thread configurations or ids landing in the wrong
  # query-point slot. Sizes are chosen large enough that a 4-thread run
  # actually splits work (n=4000 data, 200 queries).
  skip_if_not_omp()
  skip_on_cran()

  data  <- toy_data(4000L, 4L, seed = HV_SEED)
  query <- toy_data(200L, 4L, seed = HV_SEED + 11L)

  run_at_threads <- function(nt) {
    hypervolume:::omp_set_threads(nt)
    on.exit(hypervolume:::omp_set_threads(1L), add = TRUE)
    list(
      counts_075 = as.numeric(quietly(
        hypervolume:::evalfspherical(data, 0.75, query, verbose = FALSE)
      )),
      counts_150 = as.numeric(quietly(
        hypervolume:::evalfspherical(data, 1.50, query, verbose = FALSE)
      ))
    )
  }

  serial   <- run_at_threads(1L)
  parallel <- run_at_threads(4L)

  # Exact equality: counts are integer-valued and the algorithm is
  # deterministic per query point regardless of thread count.
  expect_identical(parallel$counts_075, serial$counts_075)
  expect_identical(parallel$counts_150, serial$counts_150)
})
