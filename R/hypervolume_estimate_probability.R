hypervolume_estimate_probability <- function(hv, points, reduction.factor = 1,
                                             weight.exponent = -1, set.edges.zero = TRUE,
                                             edges.zero.distance.factor = 1,
                                             parallel = FALSE, n.cores = 1,
                                             verbose = TRUE, ...) {
  np = nrow(hv@RandomPoints)
  dimhv = ncol(hv@RandomPoints)
  dimp = ncol(points)
  if (dimp != dimhv) {
    stop("Dimensionality of hypervolume and points is not the same.")
  }
  if (reduction.factor <= 0 | reduction.factor > 1) {
    stop("Reduction factor is not in (0,1].")
  }
  if (parallel) {
    if (!(is(n.cores,"integer") | is(n.cores,"numeric"))) {
      stop("If specified n.cores needs to be an integer.")
    }
  }
  numpointstokeep_hv = floor(np * reduction.factor)
  index_ss = sample(1:np, size = numpointstokeep_hv, replace = FALSE)
  hv_points_ss = hv@RandomPoints[index_ss, , drop = FALSE]
  prob_ss = hv@ValueAtRandomPoints[index_ss]
  point_density = nrow(hv_points_ss)/hv@Volume
  cutoff_dist = point_density^(-1/dimhv)

  # any NA test points?
  naIdx <- which(apply(points, 1, anyNA))
  if (length(naIdx) > 0) {
    points <- as.matrix(points[-naIdx, , drop = FALSE])
  } else {
    points <- as.matrix(points)
  }
  if (verbose) {
    cat(sprintf("Retaining %d/%d hypervolume random points for comparison with %d test points.\n",
                nrow(hv_points_ss), nrow(hv@RandomPoints), nrow(points)))
  }
  dimnames(points) <- list(NULL, NULL)

  n_pts <- nrow(points)
  n_hv  <- nrow(hv_points_ss)
  cutoff_thresh <- cutoff_dist * edges.zero.distance.factor

  # Was: n_pts separate pdist::pdist() calls, each constructing its own S4
  # object to compute one row of distances against hv_points_ss. Now: one
  # fastPdist2 (BLAS dgemm) per block of test points. Block sized so the
  # n_block x n_hv distance matrix stays under ~50 MB, matching the same
  # bound used in calculate_density.
  #
  # The `parallel = TRUE` path is silently ignored -- the vectorised path
  # is faster than the parallel one used to be, and the parallel overhead
  # (cluster startup + serialization) mostly consumed its own gains.
  target_bytes <- 50 * 1024 * 1024
  block_rows <- max(1L, min(n_pts, as.integer(target_bytes / max(1L, n_hv) / 8L)))

  # weight.exponent = -1 is the documented default and the ~universal case.
  # d^(-1) goes via pow(), which is roughly 10x a divide; special-case it.
  neg_one <- isTRUE(all.equal(weight.exponent, -1))

  probabilities <- numeric(n_pts)

  if (verbose && n_pts > 0L) {
    pb <- progress::progress_bar$new(total = n_pts)
    pb$tick(0)
  }

  for (start in seq.int(1L, max(n_pts, 1L), by = block_rows)) {
    if (n_pts == 0L) break
    end <- min(start + block_rows - 1L, n_pts)
    D <- fastPdist2(points[start:end, , drop = FALSE], hv_points_ss)

    if (neg_one) {
      W <- 1 / D
    } else {
      W <- D^weight.exponent
    }

    num <- as.numeric(W %*% prob_ss)
    den <- rowSums(W)
    res <- num / den

    if (set.edges.zero) {
      # min distance per test point across the (kept) hv random points.
      # apply() over ~2000 rows is fine here; the block size caps the
      # per-call cost and rowSums-style helpers for min don't exist in base.
      min_d <- apply(D, 1L, min)
      res[min_d > cutoff_thresh] <- 0
    }

    probabilities[start:end] <- res

    if (verbose && n_pts > 0L) pb$tick(end - start + 1L)
  }

  attr(probabilities, "NAIdx") <- naIdx
  probabilities
}
