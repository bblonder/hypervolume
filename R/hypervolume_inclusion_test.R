hypervolume_inclusion_test <- function(hv, points, reduction.factor = 1, 
                                       fast.or.accurate = "fast",
                                       fast.method.distance.factor = 1, 
                                       accurate.method.threshold = quantile(hv@ValueAtRandomPoints, 0.5), 
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
  numpointstokeep_hv = floor(np * reduction.factor)
  if (reduction.factor < 1) {
    hv_points_ss = hv@RandomPoints[sample(1:np, size = numpointstokeep_hv), 
    ]
  } else {
    hv_points_ss = hv@RandomPoints
  }
  # any NA test points?
  naIdx <- which(apply(points, 1, anyNA))
  if (length(naIdx) > 0) {
    points <- as.matrix(points[-naIdx, , drop = FALSE])
  } else {
    points <- as.matrix(points)
  }
  if (fast.or.accurate == "fast") {
    warning("Results may have a high error rate. Consider setting fast.or.accurate='accurate'.")
    if (verbose == TRUE) {
      cat(sprintf("Retaining %d/%d hypervolume random points for comparison with %d test points.\n", 
                  nrow(hv_points_ss), nrow(hv@RandomPoints), nrow(points)))
    }
    point_density = nrow(hv_points_ss)/hv@Volume
    cutoff_dist = point_density^(-1/dimhv) * fast.method.distance.factor
    if (nrow(hv_points_ss) > 0) {
      points_in_hv_all = evalfspherical(hv_points_ss, 
                                        cutoff_dist, points, verbose = verbose)
      points_in = (points_in_hv_all > 0)
      attr(points_in, "NAIdx") <- naIdx
      return(points_in)
    } else {
      warning("No points in hypervolume - increase reduction.factor or use non-empty input hypervolume! Returning NULL.")
      return(NULL)
    }
  } else if (fast.or.accurate == "accurate") {
    probabilities <- hypervolume_estimate_probability(hv, points, 
                                                      reduction.factor = reduction.factor, 
                                                      verbose = verbose, ...)
    inclusion_result <- (probabilities >= accurate.method.threshold)
    attr(inclusion_result, "NAIdx") <- naIdx
    return(inclusion_result)
  } else {
    stop("Unsupported method for fast.or.accurate.")
  }
}
