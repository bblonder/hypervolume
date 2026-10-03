# Covers proposal #8 (stop round-tripping point clouds through data.frame and
# unique.data.frame in hypervolume_set) and the evalfspherical integration that
# proposals #1/#5/#6 sit underneath.
#
# Ground truth comes from axis-aligned boxes, whose intersection and union
# volumes are known analytically.

unit_box_hv <- function(lower, upper, n = 4000L, seed = HV_SEED, name = "box") {
  make_box_hv(lower, upper, n = n, seed = seed, name = name)
}

test_that("hypervolume_set recovers a known box intersection", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  expect_s4_class(res, "HypervolumeList")
  expect_named(
    res@HVList,
    c("HV1", "HV2", "Intersection", "Union", "Unique_1", "Unique_2")
  )

  # true intersection is [0.5, 1] x [0, 1] = 0.5
  expect_volume_close(res@HVList$Intersection@Volume, 0.5, tol = 0.15)
  # true union is 1 + 1 - 0.5 = 1.5
  expect_volume_close(res@HVList$Union@Volume, 1.5, tol = 0.15)
})

test_that("union volume is exactly V1 + V2 - intersection", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  expect_equal(
    res@HVList$Union@Volume,
    hv1@Volume + hv2@Volume - res@HVList$Intersection@Volume,
    tolerance = 1e-10
  )
  expect_equal(
    res@HVList$Unique_1@Volume,
    hv1@Volume - res@HVList$Intersection@Volume,
    tolerance = 1e-10
  )
  expect_equal(
    res@HVList$Unique_2@Volume,
    hv2@Volume - res@HVList$Intersection@Volume,
    tolerance = 1e-10
  )
})

test_that("disjoint hypervolumes have (near) zero intersection", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(10, 10), c(11, 11), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  expect_lt(res@HVList$Intersection@Volume, 0.02)
  expect_equal(res@HVList$Union@Volume, 2, tolerance = 0.05)
})

test_that("self-intersection recovers (almost) the whole volume", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED + 99L, name = "A copy")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  # Not exactly 1: points are matched within a characteristic distance, so a
  # thin shell at the boundary is lost.
  expect_volume_close(res@HVList$Intersection@Volume, 1, tol = 0.15)
})

test_that("set output point clouds are numeric matrices with the right shape", {
  # hypervolume_set currently builds these via
  # matrix(as.matrix(as.data.frame(x)), ncol = dim). Proposal #8 removes the
  # round trip; the resulting object must be unchanged in type and dimnames.
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  for (nm in c("Intersection", "Union", "Unique_1", "Unique_2")) {
    pts <- res@HVList[[nm]]@RandomPoints
    expect_true(is.matrix(pts), info = nm)
    expect_type(pts, "double")
    expect_equal(ncol(pts), 2L, info = nm)
    expect_equal(colnames(pts), c("X1", "X2"), info = nm)
    expect_false(anyNA(pts), info = nm)
    expect_equal(
      length(res@HVList[[nm]]@ValueAtRandomPoints),
      nrow(pts),
      info = nm
    )
  }
})

test_that("set output point clouds contain no duplicate rows", {
  # The unique() calls proposal #8 questions. If they are removed, this test
  # tells you whether they were load-bearing.
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))

  for (nm in c("Intersection", "Union")) {
    pts <- res@HVList[[nm]]@RandomPoints
    expect_equal(sum(duplicated(pts)), 0L, info = nm)
  }
})

test_that("hypervolume_set returns NULL when check.memory is left on", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 200L, seed = HV_SEED)
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), n = 200L, seed = HV_SEED + 1L)

  expect_null(quietly(hypervolume_set(hv1, hv2, check.memory = TRUE, verbose = FALSE)))
})

test_that("hypervolume_set rejects mismatched dimensionality", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 200L)
  hv2 <- unit_box_hv(c(0, 0, 0), c(1, 1, 1), n = 200L)

  expect_error(
    quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE)),
    "not the same"
  )
})

test_that("an empty hypervolume short-circuits to the documented result", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 500L)
  empty <- hv1
  empty@Name <- "empty"
  empty@Volume <- 0
  empty@PointDensity <- 0
  empty@RandomPoints <- matrix(
    numeric(0), nrow = 0, ncol = 2,
    dimnames = list(NULL, c("X1", "X2"))
  )
  empty@ValueAtRandomPoints <- numeric(0)

  # testthat's expect_warning(expr, regexp) returns the caught warning
  # condition when regexp is supplied, not the expression value -- so we
  # can't reuse the result. Capture the value separately.
  expect_warning(
    quietly(hypervolume_set(empty, hv1, check.memory = FALSE, verbose = FALSE)),
    "no random points"
  )
  res <- suppressWarnings(
    quietly(hypervolume_set(empty, hv1, check.memory = FALSE, verbose = FALSE))
  )

  expect_equal(res@HVList$Intersection@Volume, 0)
  expect_equal(res@HVList$Union@Volume, hv1@Volume)
  expect_equal(res@HVList$Unique_2@Volume, hv1@Volume)
})

test_that("overlap statistics follow from the set volumes", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED, name = "A")
  hv2 <- unit_box_hv(c(0.5, 0), c(1.5, 1), seed = HV_SEED + 1L, name = "B")

  res <- quietly(hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE))
  stats <- hypervolume_overlap_statistics(res)

  expect_named(stats, c("jaccard", "sorensen", "frac_unique_1", "frac_unique_2"))
  expect_equal(
    unname(stats[["jaccard"]]),
    res@HVList$Intersection@Volume / res@HVList$Union@Volume,
    tolerance = 1e-12
  )
  expect_equal(
    unname(stats[["sorensen"]]),
    2 * res@HVList$Intersection@Volume / (hv1@Volume + hv2@Volume),
    tolerance = 1e-12
  )
  # For these two boxes: jaccard ~ 0.5/1.5 = 1/3
  expect_gt(stats[["jaccard"]], 0.2)
  expect_lt(stats[["jaccard"]], 0.45)
})

test_that("hypervolume_distance centroid matches a direct calculation", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), seed = HV_SEED)
  hv2 <- unit_box_hv(c(4, 0), c(5, 1), seed = HV_SEED + 1L)

  d <- hypervolume_distance(hv1, hv2, type = "centroid")
  expected <- sqrt(sum((colMeans(hv1@RandomPoints) - colMeans(hv2@RandomPoints))^2))

  expect_equal(d, expected, tolerance = 1e-10)
  expect_equal(d, 4, tolerance = 0.1)
})

test_that("hypervolume_distance minimum uses fastPdist2 consistently", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 300L, seed = HV_SEED)
  hv2 <- unit_box_hv(c(4, 0), c(5, 1), n = 300L, seed = HV_SEED + 1L)

  d <- quietly(hypervolume_distance(
    hv1, hv2, type = "minimum", num.points.max = 200, check.memory = FALSE
  ))

  expect_true(is.finite(d))
  expect_gt(d, 2.8)
  expect_lt(d, 4.1)
  # the minimum cross distance can never exceed the centroid distance
  expect_lte(d, hypervolume_distance(hv1, hv2, type = "centroid") + 1e-9)
})

test_that("hypervolume_distance rejects an unknown type", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 100L)
  expect_error(
    hypervolume_distance(hv1, hv1, type = "nonsense"),
    "unrecognized value"
  )
})

test_that("hypervolume_join concatenates into a HypervolumeList", {
  hv1 <- unit_box_hv(c(0, 0), c(1, 1), n = 100L, name = "A")
  hv2 <- unit_box_hv(c(2, 2), c(3, 3), n = 100L, name = "B")

  joined <- hypervolume_join(hv1, hv2)
  expect_s4_class(joined, "HypervolumeList")
  expect_length(joined@HVList, 2L)
  expect_equal(unname(get_volume(joined)), c(hv1@Volume, hv2@Volume))
})
