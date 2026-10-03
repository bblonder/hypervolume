# Shared fixtures and independent reference implementations.
#
# Everything here is deliberately written WITHOUT reusing package internals, so
# that a test failure means the package changed, not that the reference changed
# with it. Where a reference re-derives a quantity the package also computes,
# it does so by the most naive route available (explicit loops, dist(), etc.).

# ---- determinism ------------------------------------------------------------

HV_SEED <- 20260925L

# Run expr under a fixed seed and restore the caller's RNG state afterwards, so
# test order can never affect a later test's random draws.
with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(".Random.seed", envir = globalenv())),
      add = TRUE
    )
  }
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  force(expr)
}

# Tests that take more than a second or two. Off by default; opt in with
#   Sys.setenv(HYPERVOLUME_TEST_SLOW = "true")
skip_if_not_slow <- function() {
  if (!identical(tolower(Sys.getenv("HYPERVOLUME_TEST_SLOW", "false")), "true")) {
    skip("slow test; set HYPERVOLUME_TEST_SLOW=true to run")
  }
  invisible(TRUE)
}

# Gate for tests that exercise the OpenMP-parallelised query paths introduced
# by proposal #6. Skipped when the package was not built with OpenMP, or when
# the runtime helper for changing thread count is not present (both of those
# ship together with proposal #6). Same env-var idiom as the slow gate: set
#   Sys.setenv(HYPERVOLUME_TEST_OMP = "true")
# in the CI job that builds with SHLIB_OPENMP_CXXFLAGS.
skip_if_not_omp <- function() {
  if (!identical(tolower(Sys.getenv("HYPERVOLUME_TEST_OMP", "false")), "true")) {
    skip("OpenMP test; set HYPERVOLUME_TEST_OMP=true after proposal #6 lands")
  }
  if (!exists("omp_set_threads", envir = asNamespace("hypervolume"), inherits = FALSE)) {
    skip("hypervolume:::omp_set_threads() not exported; add it with proposal #6")
  }
  invisible(TRUE)
}

# Quieten the cat()/message() chatter that most of these functions emit.
# Deliberately does NOT suppress warnings: several tests assert on them, and
# swallowing them here would make those assertions silently unreachable.
quietly <- function(expr) {
  suppressMessages(utils::capture.output(res <- force(expr)))
  res
}

# ---- toy datasets -----------------------------------------------------------

# Small, well-separated, fully deterministic. Named dimensions throughout,
# because several code paths behave differently when dimnames are NULL.
toy_data <- function(n = 25L, d = 2L, seed = HV_SEED) {
  m <- with_seed(seed, matrix(stats::rnorm(n * d), nrow = n, ncol = d))
  dimnames(m) <- list(NULL, paste0("X", seq_len(d)))
  m
}

# Points drawn uniformly from an axis-aligned box, so the true volume is known
# exactly: prod(upper - lower).
box_points <- function(n, lower, upper, seed = HV_SEED) {
  d <- length(lower)
  stopifnot(length(upper) == d)
  m <- with_seed(seed, matrix(stats::runif(n * d), nrow = n, ncol = d))
  for (j in seq_len(d)) {
    m[, j] <- lower[j] + m[, j] * (upper[j] - lower[j])
  }
  dimnames(m) <- list(NULL, paste0("X", seq_len(d)))
  m
}

# A Hypervolume object over a known box. Bypasses the estimators entirely so
# that set-operation tests have an analytic ground truth to compare against.
make_box_hv <- function(lower, upper, n = 3000L, seed = HV_SEED,
                        name = "box", value = 1) {
  pts <- box_points(n, lower, upper, seed = seed)
  vol <- prod(upper - lower)
  methods::new(
    "Hypervolume",
    Name = name,
    Method = "test fixture",
    Data = pts,
    Dimensionality = length(lower),
    Volume = vol,
    PointDensity = nrow(pts) / vol,
    Parameters = list(),
    RandomPoints = pts,
    ValueAtRandomPoints = rep(value, nrow(pts))
  )
}

# ---- independent reference implementations ----------------------------------

# Euclidean cross-distance matrix, the slow and obvious way.
ref_pdist <- function(A, B) {
  A <- as.matrix(A); B <- as.matrix(B)
  out <- matrix(NA_real_, nrow = nrow(A), ncol = nrow(B))
  for (i in seq_len(nrow(A))) {
    for (j in seq_len(nrow(B))) {
      out[i, j] <- sqrt(sum((A[i, ] - B[j, ])^2))
    }
  }
  out
}

# For each row of `points`, how many rows of `data` lie within `radius`.
# The kd-tree uses `distance_squared <= radius^2`, i.e. a CLOSED ball.
ref_ball_counts <- function(data, points, radius) {
  d <- ref_pdist(points, data)
  as.integer(rowSums(d <= radius))
}

# Indices (1-based) of the data rows inside each ball.
ref_ball_ids <- function(data, points, radius) {
  d <- ref_pdist(points, data)
  lapply(seq_len(nrow(d)), function(i) which(d[i, ] <= radius))
}

# For each row of `points`, how many rows of `data` fall in the axis-aligned
# box of half-width `bandwidth`. range_query uses a CLOSED box.
ref_box_counts <- function(data, points, bandwidth) {
  data <- as.matrix(data); points <- as.matrix(points)
  vapply(seq_len(nrow(points)), function(i) {
    inside <- rep(TRUE, nrow(data))
    for (j in seq_len(ncol(data))) {
      lo <- points[i, j] - bandwidth[j]
      hi <- points[i, j] + bandwidth[j]
      inside <- inside & data[, j] >= lo & data[, j] <= hi
    }
    sum(inside)
  }, integer(1))
}

# Mixture-of-diagonal-Gaussians density, built from univariate dnorm only.
# Deliberately does NOT use mvtnorm, so it is independent of the package's
# dmvnorm call. Historically had a `divide_by_n` argument that mirrored
# the package's trailing `/ n_points`; that trailing division was a
# double-normalisation bug and has been removed from calculate_density.
# The argument is retained (default FALSE) so any test that was pinning
# the buggy behaviour can still opt in explicitly if needed.
ref_density <- function(s, means, kernel_sd, weight, divide_by_n = FALSE) {
  s <- as.matrix(s); means <- as.matrix(means)
  d <- ncol(means)
  n <- nrow(means)
  out <- numeric(nrow(s))
  for (i in seq_len(n)) {
    dens <- rep(1, nrow(s))
    for (j in seq_len(d)) {
      dens <- dens * stats::dnorm(s[, j], mean = means[i, j], sd = kernel_sd[j])
    }
    out <- out + weight[i] * dens
  }
  if (divide_by_n) out <- out / n
  out
}

# ---- assertion helpers ------------------------------------------------------

# Volume comparisons on Monte Carlo estimates need proportional tolerance.
expect_volume_close <- function(observed, expected, tol = 0.10) {
  rel <- abs(observed - expected) / expected
  expect_lt(rel, tol)
}

expect_valid_hypervolume <- function(hv) {
  expect_s4_class(hv, "Hypervolume")
  expect_true(is.matrix(hv@RandomPoints))
  expect_equal(ncol(hv@RandomPoints), hv@Dimensionality)
  expect_equal(length(hv@ValueAtRandomPoints), nrow(hv@RandomPoints))
  expect_true(is.finite(hv@Volume))
  expect_gte(hv@Volume, 0)
  expect_false(anyNA(hv@RandomPoints))
  invisible(hv)
}
