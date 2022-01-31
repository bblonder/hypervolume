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
    if (!class(n.cores) %in% c("integer", "numeric")) {
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
  probabilities <- rep(NA, nrow(points))
  if (parallel) {
    cat("Processing points in parallel...\n")
    pbapply::pboptions(style = 3, nout = 20, type = "timer")
    cls <- makeCluster(n.cores)
    clusterExport(cls, c("points", "hv_points_ss", "weight.exponent",
                         "set.edges.zero", "cutoff_dist", "prob_ss",
                         "edges.zero.distance.factor"),
                  envir = environment())
    probabilities <- pbapply::pbsapply(1:nrow(points), function(i) {
      distances <- pdist::pdist(points[i, , drop = FALSE], 
                                hv_points_ss)@dist
      weights <- distances^weight.exponent
      result <- sum(prob_ss * weights)/sum(weights)
      if (set.edges.zero == TRUE) {
        if (min(distances) > cutoff_dist * edges.zero.distance.factor) {
          result <- 0
        }
      }
      return(result)
    }, cl = cls)
    parallel::stopCluster(cls)
  } else {
    pb <- progress::progress_bar$new(total = nrow(points))
    if (verbose == TRUE) {
      pb$tick(0)
    }
    probabilities <- sapply(1:nrow(points), function(i) {
      if (verbose == TRUE & (i%%10 == 0)) {
        if (!pb$finished == TRUE) {
          pb$update(i/nrow(points))
        }
      }
      distances <- pdist::pdist(points[i, , drop = FALSE], hv_points_ss)@dist
      weights <- distances^weight.exponent
      result <- sum(prob_ss * weights)/sum(weights)
      if (set.edges.zero == TRUE) {
        if (min(distances) > cutoff_dist * edges.zero.distance.factor) {
          result <- 0
        }
      }
      return(result)
    })
  }
  attr(probabilities, "NAIdx") <- naIdx
  return(probabilities)
}
