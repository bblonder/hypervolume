find_optimal_occupancy_thin <- function(..., verbose = TRUE, sequence = seq(0, 1, 0.1), n = 10, res_type = "raw") {

  # get random state from global environment
  old_state <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  on.exit(
    {
      # restore random state
      if (!is.null(old_state)) {
        assign(".Random.seed", old_state, envir = .GlobalEnv, inherits = FALSE)
      }
    },
    add = TRUE
  )

  store_list <- list()
  seq_expanded <- rep(sequence, each = n)
  seed_set <- rep(sample(1:(.Machine$integer.max), n), times = length(sequence))

  if (verbose) {
    pb <- progress_bar$new(total = length(seed_set))
  }

  for (i in 1:length(seq_expanded)) {
    if (seq_expanded[i] == 0) {
      temp <- hypervolume_n_occupancy(..., thin = FALSE, verbose = FALSE, print_log = TRUE, seed = seed_set[i])
    } else {
      temp <- hypervolume_n_occupancy(..., thin = TRUE, verbose = FALSE, quant.thin = seq_expanded[i], print_log = TRUE, seed = seed_set[i])
    }
    store_list[[i]] <- read.table("log_occupancy.txt", header = TRUE)

    if (verbose) {
      pb$tick()
    }
  }

  if (verbose) {
    pb$terminate()
  }

  res_fun <- sapply(store_list, function(x) sqrt(mean((x[, "input"] - x[, "re_computed"])^2)))
  res_fun <- data.frame(seed_set, seq_expanded, res_fun)
  colnames(res_fun) <- c("seed", "quant.thin", "rmse")

  if (identical(res_type, "raw")) {
    res_final <- res_fun[order(res_fun[, "rmse"]), ]
  }

  if (identical(res_type, "summary")) {
    res_final <- aggregate(rmse ~ quant.thin, res_fun, FUN = function(x) c(mn = mean(x), sd = sd(x)))
  }

  rownames(res_final) <- NULL
  res_final
}
