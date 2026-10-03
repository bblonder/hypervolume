# Covers proposal #9 (fastPdist2 numerical stability).
#
# fastPdist2 uses the |a|^2 + |b|^2 - 2a.b expansion, which is fast but can
# produce small negative values for near-coincident points, and sqrt() of those
# is NaN. These tests pin the current behaviour so the clamp can be added
# without changing anything else.

test_that("fastPdist2 matches a naive distance calculation", {
  A <- toy_data(12L, 3L, seed = HV_SEED)
  B <- toy_data(9L, 3L, seed = HV_SEED + 1L)

  observed <- hypervolume:::fastPdist2(A, B)
  expected <- ref_pdist(A, B)

  expect_equal(dim(observed), c(nrow(A), nrow(B)))
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-8)
})

test_that("fastPdist2 handles single-row inputs on both sides", {
  A <- toy_data(1L, 4L, seed = HV_SEED)
  B <- toy_data(6L, 4L, seed = HV_SEED + 2L)

  expect_equal(dim(hypervolume:::fastPdist2(A, B)), c(1L, 6L))
  expect_equal(dim(hypervolume:::fastPdist2(B, A)), c(6L, 1L))
  expect_equal(
    as.numeric(hypervolume:::fastPdist2(A, B)),
    as.numeric(t(hypervolume:::fastPdist2(B, A))),
    tolerance = 1e-8
  )
})

test_that("fastPdist2 is symmetric and zero on the diagonal for A vs A", {
  A <- toy_data(10L, 2L, seed = HV_SEED)
  D <- hypervolume:::fastPdist2(A, A)

  expect_equal(as.numeric(D), as.numeric(t(D)), tolerance = 1e-8)
  # The diagonal is the exact case where catastrophic cancellation bites:
  # |a|^2 + |a|^2 - 2a.a can land either side of zero. Proposal #9's clamp
  # guarantees the sqrt is real, so no NaN is tolerable.
  expect_false(anyNA(diag(D)))
  expect_true(all(abs(diag(D)) < 1e-6))
})

test_that("fastPdist2 returns real distances for near-coincident points", {
  # Acceptance test for proposal #9 (clamp negative squared distances
  # before sqrt). Points separated by far less than sqrt(.Machine$double.eps)
  # relative to their magnitude: at |a| ~ 1e6 the squared norms are ~3e12,
  # where one ulp is ~1e-3, but the true squared distance is ~4e-18. So
  # |a|^2 + |b|^2 - 2a.b is pure rounding noise of either sign; without the
  # clamp, sqrt() of the negative half returns NaN and hypervolume_distance
  # (type = "minimum") silently hides it via na.rm = TRUE.
  base <- matrix(rep(1e6, 3), nrow = 1)
  A <- rbind(base, base + 1e-9, base - 1e-9)
  colnames(A) <- paste0("X", 1:3)

  D <- hypervolume:::fastPdist2(A, A)

  expect_false(anyNA(D))
  expect_true(all(D >= 0))
})

test_that("fastPdist2 preserves scale across dimensionalities", {
  for (d in c(1L, 2L, 5L, 8L)) {
    A <- toy_data(7L, d, seed = HV_SEED + d)
    B <- toy_data(5L, d, seed = HV_SEED + 100L + d)
    expect_equal(
      as.numeric(hypervolume:::fastPdist2(A, B)),
      as.numeric(ref_pdist(A, B)),
      tolerance = 1e-8,
      info = paste("d =", d)
    )
  }
})
