hypervolume_to_data_frame <- function(hv, remove_zeroes = TRUE) {
  if (inherits(hv, "Hypervolume")) {
    res <- data.frame(hv@RandomPoints, ValueAtRandomPoints = hv@ValueAtRandomPoints)
    if (remove_zeroes) {
      res <- res[res[, "ValueAtRandomPoints"] != 0, ]
    }
  }

  if (inherits(hv, "HypervolumeList")) {
    hyper_methods <- sapply(hv@HVList, function(x) x@Method)
    hyper_rp <- lapply(hv@HVList, function(x) x@RandomPoints)
    hyper_vrp <- lapply(hv@HVList, function(x) x@ValueAtRandomPoints)
    hyper_lab <- lapply(hv@HVList, function(x) x@Name)
    hyper_lab_length <- lapply(hyper_vrp, length)
    hyper_lab <- mapply(function(x, y) rep(x, y), hyper_lab, hyper_lab_length, SIMPLIFY = FALSE)
    hyper_rp <- do.call(rbind, hyper_rp)
    hyper_vrp <- do.call(c, hyper_vrp)
    hyper_lab <- do.call(c, hyper_lab)
    res <- data.frame(
      Name = hyper_lab,
      hyper_rp, ValueAtRandomPoints = hyper_vrp
    )
    if(any(hyper_methods %in% c("n_occupancy", "n_occupancy_test", "occupancy_to_union", "occupancy_to_intersection", "occupancy_to_unshared"))) {
      if (remove_zeroes) {
        res <- res[res[, "ValueAtRandomPoints"] != 0, ]
      }
    }
  }

  rownames(res) <- NULL
  res
}
