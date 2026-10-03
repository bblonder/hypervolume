# Scenario definitions shared by the golden-snapshot generator
# (tools/make_golden.R) and the golden-snapshot test (test-10-golden.R).
#
# Each scenario is a zero-argument function returning a named list of numeric
# values. Keep them cheap, deterministic, and additive: append new scenarios
# rather than editing existing ones, so an old fixture stays comparable.

golden_scenarios <- function() {
  list(

    # --- proposal #2: the density evaluation itself --------------------------
    calculate_density = function() {
      means <- toy_data(15L, 3L, seed = HV_SEED)
      s     <- toy_data(9L, 3L, seed = HV_SEED + 4L)
      bw    <- c(0.5, 0.7, 0.9)
      w     <- rep(1 / nrow(means), nrow(means))
      list(density = hypervolume:::calculate_density(
        s = s, means = means, kernel_sd = bw, weight = w, verbose = FALSE
      ))
    },

    # --- proposal #1/#5/#6: kd-tree query results ---------------------------
    ball_counts = function() {
      data <- toy_data(60L, 3L, seed = HV_SEED)
      pts  <- toy_data(25L, 3L, seed = HV_SEED + 7L)
      list(
        r_0_75 = as.numeric(quietly(hypervolume:::evalfspherical(
          data, 0.75, pts, verbose = FALSE
        ))),
        r_1_50 = as.numeric(quietly(hypervolume:::evalfspherical(
          data, 1.50, pts, verbose = FALSE
        )))
      )
    },

    # --- proposal #3: threshold statistics ----------------------------------
    threshold_stats = function() {
      pts <- box_points(400L, c(0, 0), c(1, 1), seed = HV_SEED)
      vals <- with_seed(HV_SEED + 1L, stats::runif(400L))
      hv <- methods::new(
        "Hypervolume",
        Name = "golden", Method = "test fixture",
        Data = pts, Dimensionality = 2,
        Volume = 1, PointDensity = 400, Parameters = list(),
        RandomPoints = pts, ValueAtRandomPoints = vals
      )
      res <- quietly(hypervolume_threshold(
        hv, num.thresholds = 25, plot = FALSE, verbose = FALSE
      ))
      list(
        volume = res$Statistics$volume,
        volume_quantile = res$Statistics$volume_quantile,
        probability_quantile = res$Statistics$probability_quantile
      )
    },

    # --- proposal #7: sampler volume accounting -----------------------------
    sampler_volume = function() {
      data <- toy_data(5L, 2L, seed = HV_SEED)
      res <- quietly(with_seed(HV_SEED, hypervolume:::sample_model_ellipsoid(
        predict_function = function(x) rep(1, nrow(x)),
        data = data, scales = c(1, 1), min.value = 0,
        samples.per.point = 100, chunk.size = 400, verbose = FALSE
      )))
      list(volume = res$volume, n_samples = nrow(res$samples))
    },

    # --- proposals #2 + #3 end to end ---------------------------------------
    gaussian_2d = function() {
      data <- toy_data(12L, 2L, seed = HV_SEED)
      hv <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
        data,
        kde.bandwidth = estimate_bandwidth(data, method = "fixed", value = 0.6),
        samples.per.point = 40, chunk.size = 300, verbose = FALSE
      )))
      list(
        volume = hv@Volume,
        n_points = nrow(hv@RandomPoints),
        point_density = hv@PointDensity,
        value_quantiles = as.numeric(stats::quantile(
          hv@ValueAtRandomPoints, probs = c(0, 0.25, 0.5, 0.75, 1)
        )),
        centroid = as.numeric(colMeans(hv@RandomPoints))
      )
    },

    gaussian_3d = function() {
      data <- toy_data(10L, 3L, seed = HV_SEED + 3L)
      hv <- with_seed(HV_SEED, quietly(hypervolume_gaussian(
        data,
        kde.bandwidth = estimate_bandwidth(data, method = "fixed", value = 0.7),
        samples.per.point = 30, chunk.size = 300, verbose = FALSE
      )))
      list(volume = hv@Volume, n_points = nrow(hv@RandomPoints))
    },

    # --- proposal #8: set operations ----------------------------------------
    set_operations = function() {
      hv1 <- make_box_hv(c(0, 0), c(1, 1), n = 4000L, seed = HV_SEED, name = "A")
      hv2 <- make_box_hv(c(0.5, 0), c(1.5, 1), n = 4000L, seed = HV_SEED + 1L, name = "B")
      res <- quietly(with_seed(HV_SEED, hypervolume_set(
        hv1, hv2, check.memory = FALSE, verbose = FALSE
      )))
      list(volumes = unname(get_volume(res)))
    },

    # --- proposal #4: probability estimation --------------------------------
    estimate_probability = function() {
      pts <- box_points(500L, c(0, 0), c(1, 1), seed = HV_SEED)
      hv <- methods::new(
        "Hypervolume",
        Name = "golden", Method = "test fixture",
        Data = pts, Dimensionality = 2,
        Volume = 1, PointDensity = 500, Parameters = list(),
        RandomPoints = pts, ValueAtRandomPoints = 0.2 + 0.6 * pts[, 1]
      )
      query <- rbind(c(0.25, 0.5), c(0.5, 0.5), c(0.75, 0.5))
      colnames(query) <- c("X1", "X2")
      # Seeded like the other scenarios: hypervolume_estimate_probability calls
      # sample() internally, which would otherwise perturb the global stream and
      # reorder the summation.
      list(prob = as.numeric(with_seed(HV_SEED, quietly(
        hypervolume_estimate_probability(hv, query, verbose = FALSE)
      ))))
    },

    # --- proposal #9: distance matrix ---------------------------------------
    fastpdist = function() {
      A <- toy_data(12L, 3L, seed = HV_SEED)
      B <- toy_data(9L, 3L, seed = HV_SEED + 1L)
      list(distances = as.numeric(hypervolume:::fastPdist2(A, B)))
    }
  )
}

golden_fixture_path <- function() {
  testthat::test_path("fixtures", "golden-baseline.rds")
}
