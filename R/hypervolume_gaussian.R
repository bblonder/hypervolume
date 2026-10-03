# s is a test point
# means is the centers of each data point
# kernel_sd is the kernel bandwidth (in natural units, not squared, diagonal only)
calculate_density <- function(s, means, kernel_sd, weight=NULL, chunksize=10, verbose=TRUE) {

  n_points <- nrow(means)

  # Preserves a documented trap: the original body does `weight[i] * dmvnorm(...)`,
  # and NULL[i] is NULL, so the accumulator collapses to numeric(0) on the first
  # iteration. hypervolume_gaussian always fills weight; only test-04 relies on
  # this branch. Preserve rather than fix so a rewrite of the semantics is its
  # own commit -- see TESTING.md's known-defects list.
  if (is.null(weight)) return(numeric(0))

  d <- ncol(means)
  s <- as.matrix(s)
  n_s <- nrow(s)

  # The diagonal covariance factorises: k(z, m) = prod_j dnorm(z_j; m_j, h_j).
  # Scale coordinates once by 1/kernel_sd so the exponent becomes
  #   -0.5 * sum_j ((z_j - m_j) / h_j)^2 = -0.5 * ||Z - M||^2
  # and pull the normalising constant out of the sum.
  Z <- sweep(s,     2L, kernel_sd, "/")
  M <- sweep(means, 2L, kernel_sd, "/")
  norm_const <- exp(-0.5 * d * log(2 * pi) - sum(log(kernel_sd)))

  # Block over query rows to keep the n_s x n_points distance matrix
  # bounded in memory (~ n_s_block * n_points * 8 bytes). At the flagship
  # scale (samples_retained ~ 1e4-1e5, n_points ~ 1e2-1e4) an unbounded
  # matrix runs to hundreds of MB per call. Block size chosen so a single
  # D2 chunk stays under ~50 MB. Result is bit-identical to the unblocked
  # form because each output row is a fully independent inner product.
  target_bytes <- 50 * 1024 * 1024
  block_rows <- max(1L, min(n_s, as.integer(target_bytes / max(1L, n_points) / 8L)))

  # chunksize is preserved as a signature parameter but drives only the
  # progress bar's tick count (matching the number of data-point chunks
  # the original code paced over), so callers passing non-default values
  # -- and test-04's chunksize-invariance pin -- see identical results.
  num_chunks <- max(1L, ceiling(n_points / chunksize))
  if (verbose) {
    pb <- progress_bar$new(total = num_chunks)
    pb$tick(0)
  }

  density <- numeric(n_s)
  for (start in seq.int(1L, n_s, by = block_rows)) {
    end <- min(start + block_rows - 1L, n_s)
    D2 <- fastPdist2(Z[start:end, , drop = FALSE], M)^2
    density[start:end] <- as.numeric(exp(-0.5 * D2) %*% weight)
  }

  if (verbose) {
    for (i in seq_len(num_chunks)) pb$tick()
  }

  # weight is a mixture-weight vector (sum-to-1 in the hypervolume_gaussian
  # path), so sum_i weight_i * k(x; mean_i, sd) is already the mixture
  # density and integrates to 1. Blonder 2014 Box 2 fixes the normalising
  # constant as alpha = 1 / (m * prod(h)); dmvnorm-style k already carries
  # 1/prod(h), and weight already carries 1/m via the caller's normalisation
  # -- so no further division is needed. The previous body had a trailing
  # `/ n_points` which was a second 1/m, producing values 1/n times too
  # small (see PROPOSAL2_VERIFICATION.md sec 1). Fixed here; downstream
  # consumers reading @ValueAtRandomPoints in absolute terms
  # (hypervolume_estimate_probability, ValueAtRandomPoints comparisons
  # across hvs with different n) now see actual probability density values.
  # Volume is unaffected because acceptance is deterministic (density > 0)
  # and thresholding is quantile-based.
  density * norm_const
}


hypervolume_gaussian <- function(data, name=NULL, weight=NULL, samples.per.point=ceiling((10^(3+sqrt(ncol(data))))/nrow(data)), kde.bandwidth=estimate_bandwidth(data), sd.count=3, quantile.requested=0.95, quantile.requested.type="probability", chunk.size=1e3, verbose=TRUE, ...)
{
  data = as.matrix(data)
  d = ncol(data)
  np = nrow(data)
  
  if (is.null(weight))
  {
    weight <- rep(1/nrow(data),nrow(data))
  }
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data)[[2]] <- paste("X",1:ncol(data))
  }
  
  if (sd.count<3)
  {
    warning(sprintf("Values of sd.count (%d) is low.\nRecommended minimum value is 3, with higher values giving better performance.\nBoundaries and volumes may be inaccurate.",sd.count))
  }
  
  # do error check on kde.bandwidth vector
  if(ncol(data)!=length(kde.bandwidth))
  {
    stop('data and kde.bandwidth must have same dimensionality')
  }
  
  if (any(kde.bandwidth==0))
  {
    stop('Bandwidth must be non-zero.')
  }
  
  names(kde.bandwidth) <- dimnames(data)[[2]]
  
  if (length(weight)!=nrow(data))
  {
    stop("The length of the weights must be equal to the number of observations.")
  }
  if (abs(sum(weight)-1)>10^-10)
  {
    warning("The sum of the weights must be equal to 1. Normalizing the weights.")
    weight <- weight / sum(weight)
  }
  
  kde.method = attr(kde.bandwidth,"method")
  if (is.null(kde.method))
  {
    stop('KDE bandwidth must be set using estimate_bandwidth(), not entered manually.') 
  }
  
  predict_function_gaussian <- function(x)
  {
    return(calculate_density(s=x,means=data,kernel_sd=kde.bandwidth,weight=weight))
  }
  
  # do elliptical sampling
  samples_all = sample_model_ellipsoid(
    predict_function = predict_function_gaussian,
    data = data,
    scales = kde.bandwidth * (sd.count), # distance out to sample
    min.value = 0,
    samples.per.point = samples.per.point, 
    chunk.size=chunk.size, 
    verbose=verbose)
  
  random_points = samples_all$samples[,1:d,drop=FALSE]
  values_accepted = samples_all$samples[,d+1]
  volume <- samples_all$volume
  
  point_density = nrow(random_points) / volume

  hv_gaussian <- new("Hypervolume",
                Data=data,
                Method = 'Gaussian kernel density estimate',
                RandomPoints= random_points,
                PointDensity= point_density,
                Volume= volume,
                Dimensionality=d,
                ValueAtRandomPoints=values_accepted,
                Name=ifelse(is.null(name), "untitled", toString(name)),
                Parameters = list(kde.bandwidth=kde.bandwidth,
                                  kde.method=attr(kde.bandwidth,'method'),
                                  samples.per.point=samples.per.point, 
                                  sd.count=sd.count, 
                                  quantile.requested=quantile.requested, 
                                  quantile.requested.type=quantile.requested.type))
  
  # apply desired threshold type
  hv_gaussian_thresholded = hypervolume_threshold(hv_gaussian, 
                                                  quantile.requested=quantile.requested,
                                                  quantile.requested.type=quantile.requested.type, 
                                                  uniform.density = FALSE, # keep probability values
                                                  verbose=verbose,
                                                  num.thresholds = 1000, # give higher resolution to obtain the desired quantile
                                                  plot=FALSE, ...)[["HypervolumesThresholded"]]
  hv_gaussian_thresholded@Name = hv_gaussian@Name
  
  return(hv_gaussian_thresholded)	
}