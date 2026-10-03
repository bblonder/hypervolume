# Covers proposal #2 (replace the per-data-point dmvnorm loop with a single
# blocked matrix product).
#
# The rewrite changes summation order and the normalising constant path, so
# these tests pin: the value, the SHAPE of the value, the weight semantics
# including the suspected double normalisation, and invariance to chunking.

test_that("calculate_density matches an independent mixture-of-Gaussians", {
  means <- toy_data(15L, 3L, seed = HV_SEED)
  s     <- toy_data(9L, 3L, seed = HV_SEED + 4L)
  bw    <- c(0.5, 0.7, 0.9)
  w     <- rep(1 / nrow(means), nrow(means))

  observed <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w, verbose = FALSE
  )
  expected <- ref_density(s, means, bw, w)

  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("calculate_density returns a bare numeric vector, not a matrix", {
  # A matrix-product rewrite returns an n x 1 matrix unless it is unwrapped.
  # Downstream, sample_model_ellipsoid does `as.numeric(predicted_values) >
  # min.value` and then indexes rows with the result, so a stray matrix would
  # survive that line and fail confusingly later.
  means <- toy_data(8L, 2L)
  s     <- toy_data(5L, 2L, seed = HV_SEED + 6L)
  bw    <- c(0.6, 0.6)

  out <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw,
    weight = rep(1 / 8, 8), verbose = FALSE
  )

  expect_true(is.numeric(out))
  expect_null(dim(out))
  expect_length(out, nrow(s))
})

test_that("calculate_density output is strictly positive and finite", {
  means <- toy_data(12L, 2L)
  s     <- toy_data(20L, 2L, seed = HV_SEED + 8L)
  bw    <- c(0.4, 0.4)

  out <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw,
    weight = rep(1 / 12, 12), verbose = FALSE
  )

  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
})

test_that("calculate_density is invariant to chunksize", {
  # chunksize only drives the progress bar. It must not touch the result --
  # and the chunk bookkeeping is one of the things proposal #2 deletes.
  means <- toy_data(20L, 2L)
  s     <- toy_data(7L, 2L, seed = HV_SEED + 9L)
  bw    <- c(0.5, 0.5)
  w     <- rep(1 / 20, 20)

  base <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w,
    chunksize = 10, verbose = FALSE
  )

  for (cs in c(1, 3, 20, 1000)) {
    other <- hypervolume:::calculate_density(
      s = s, means = means, kernel_sd = bw, weight = w,
      chunksize = cs, verbose = FALSE
    )
    expect_equal(other, base, tolerance = 1e-12, info = paste("chunksize =", cs))
  }
})

test_that("calculate_density is invariant to the order of the data points", {
  means <- toy_data(18L, 2L)
  s     <- toy_data(6L, 2L, seed = HV_SEED + 12L)
  bw    <- c(0.55, 0.65)
  w     <- rep(1 / 18, 18)

  base <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w, verbose = FALSE
  )
  perm <- with_seed(HV_SEED + 13L, sample.int(nrow(means)))
  shuffled <- hypervolume:::calculate_density(
    s = s, means = means[perm, , drop = FALSE], kernel_sd = bw,
    weight = w[perm], verbose = FALSE
  )

  # Tolerance rather than exactness: any reordering of a floating-point sum
  # moves the last bits, and the BLAS rewrite reorders it further.
  expect_equal(shuffled, base, tolerance = 1e-12)
})

test_that("calculate_density respects non-uniform weights", {
  means <- rbind(c(0, 0), c(4, 4))
  colnames(means) <- c("X1", "X2")
  s <- rbind(c(0, 0), c(4, 4))
  colnames(s) <- c("X1", "X2")
  bw <- c(0.5, 0.5)

  # All weight on the first kernel: the point at (0,0) should dominate.
  w <- c(1, 0)
  out <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w, verbose = FALSE
  )
  expect_gt(out[1], out[2])
  expect_equal(out, ref_density(s, means, bw, w), tolerance = 1e-12)

  # Mirror image.
  w2 <- c(0, 1)
  out2 <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = w2, verbose = FALSE
  )
  expect_lt(out2[1], out2[2])
  expect_equal(out2, rev(out), tolerance = 1e-10)
})

test_that("calculate_density with weight=NULL returns numeric(0)", {
  # DOCUMENTS A TRAP, not desired behaviour. The signature defaults weight to
  # NULL, but the body does `weight[i] * dmvnorm(...)`, and NULL[i] is NULL, so
  # the accumulator collapses to numeric(0) and the function silently returns
  # nothing. hypervolume_gaussian never hits this because it always fills
  # weight first. Any rewrite of this function should either give weight a real
  # default of rep(1/n, n) or stop() on NULL -- at which point this test should
  # be updated to assert the new, better contract.
  means <- toy_data(6L, 2L)
  s     <- toy_data(4L, 2L, seed = HV_SEED + 14L)

  out <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = c(0.5, 0.5),
    weight = NULL, verbose = FALSE
  )
  expect_length(out, 0L)
})

test_that("calculate_density integrates to 1 as a proper probability density", {
  # After the 1/m normalisation fix: with weight summing to 1, the returned
  # mixture density k(x) = sum_i w_i * dnorm(x; mean_i, sd) integrates to 1
  # over R^d, matching Blonder 2014 Box 2's alpha = 1/(m * prod(h)).
  # Previously the body divided by n_points a second time and this test
  # pinned the buggy 1/n integral -- see PROPOSAL2_VERIFICATION.md sec 1
  # and TESTING.md's known-defects list (removed with this commit).
  grid <- matrix(seq(-12, 12, by = 0.002), ncol = 1)
  colnames(grid) <- "X1"

  for (n in c(1L, 2L, 5L)) {
    means <- matrix(seq(-1, 1, length.out = n), ncol = 1)
    colnames(means) <- "X1"
    w <- rep(1 / n, n)

    dens <- hypervolume:::calculate_density(
      s = grid, means = means, kernel_sd = 1.0, weight = w, verbose = FALSE
    )
    integral <- sum(dens) * 0.002

    expect_equal(integral, 1, tolerance = 1e-4, info = paste("n =", n))
  }
})

test_that("calculate_density scales correctly with bandwidth", {
  # Peak height of a d-dimensional diagonal Gaussian goes as 1 / prod(sd).
  means <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("X1", "X2")))
  s     <- matrix(c(0, 0), nrow = 1, dimnames = list(NULL, c("X1", "X2")))

  narrow <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = c(0.5, 0.5), weight = 1, verbose = FALSE
  )
  wide <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = c(1.0, 1.0), weight = 1, verbose = FALSE
  )

  expect_equal(narrow / wide, 4, tolerance = 1e-8)
})

test_that("calculate_density handles a single data point and a single test point", {
  means <- matrix(c(1, 2, 3), nrow = 1, dimnames = list(NULL, c("X1", "X2", "X3")))
  s     <- matrix(c(1, 2, 3), nrow = 1, dimnames = list(NULL, c("X1", "X2", "X3")))
  bw    <- c(0.5, 0.5, 0.5)

  out <- hypervolume:::calculate_density(
    s = s, means = means, kernel_sd = bw, weight = 1, verbose = FALSE
  )

  expect_length(out, 1L)
  # ignore_attr: R's 1-row-matrix column-subset (s[, 1] on a matrix whose
  # only dimnames are column names) inherits the column name, and dnorm
  # propagates it, so ref_density ends up with a stray "X1" name.
  # calculate_density itself correctly returns a bare numeric.
  expect_equal(out, ref_density(s, means, bw, 1), tolerance = 1e-12, ignore_attr = TRUE)
})
