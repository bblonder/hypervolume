# Golden-value regression tests.
#
# These compare today's numbers against a snapshot taken from the UNMODIFIED
# package (the `main` branch). They are the safety net for the performance work:
# every proposal in hypervolume_speedup_notes.md is meant to be a pure speedup,
# so every number here should survive it.
#
# The fixture is not committed by default -- generate it on your own machine,
# on the baseline commit, with:
#
#     Rscript tools/make_golden.R
#
# If the fixture is missing these tests skip rather than fail, so a fresh
# checkout still gets a green run.

test_that("golden fixture, if present, still reproduces", {
  path <- golden_fixture_path()
  skip_if(
    !file.exists(path),
    "no golden fixture; run `Rscript tools/make_golden.R` on the baseline commit"
  )

  golden <- readRDS(path)
  scenarios <- golden_scenarios()

  # Tolerances are per-scenario. Anything driven purely by arithmetic is held
  # tight; anything that consumes the RNG stream is held tight too, because the
  # seed is fixed -- a change there means the sequence of random draws moved,
  # which for a pure-performance change is itself a finding worth seeing.
  tolerances <- c(
    calculate_density     = 1e-10,
    ball_counts           = 0,
    threshold_stats       = 1e-12,
    sampler_volume        = 1e-8,
    gaussian_2d           = 1e-8,
    gaussian_3d           = 1e-8,
    set_operations        = 1e-8,
    estimate_probability  = 1e-8,
    fastpdist             = 1e-10
  )

  for (nm in names(golden$values)) {
    if (!nm %in% names(scenarios)) {
      # Fixture is older than the code. Not a failure; just nothing to compare.
      next
    }
    tol <- if (nm %in% names(tolerances)) unname(tolerances[[nm]]) else 1e-8
    observed <- scenarios[[nm]]()

    expect_equal(
      observed,
      golden$values[[nm]],
      tolerance = tol,
      info = sprintf(
        "scenario '%s' drifted from the baseline recorded on %s (pkg %s, %s)",
        nm, golden$meta$date, golden$meta$package_version, golden$meta$r_version
      )
    )
  }
})

test_that("every scenario runs and returns finite numbers", {
  # Runs regardless of whether a fixture exists, so the scenario definitions
  # themselves cannot silently rot.
  scenarios <- golden_scenarios()

  for (nm in names(scenarios)) {
    out <- scenarios[[nm]]()
    expect_type(out, "list")
    expect_gt(length(out), 0L)
    for (el in names(out)) {
      v <- out[[el]]
      expect_true(is.numeric(v), info = paste(nm, el))
      expect_true(all(is.finite(v)), info = paste(nm, el))
    }
  }
})
