hypervolume_threshold <- function(hv, thresholds=NULL, num.thresholds=20, quantile.requested=NULL,quantile.requested.type="volume",uniform.density=TRUE, plot=TRUE, verbose=TRUE)
{
  values <- hv@ValueAtRandomPoints
  n_total <- length(values)

  if (is.null(thresholds))
  {
    thresholds <- seq(min(values), max(values), length.out=num.thresholds)
  }

  # Statistics for every threshold from a single sort. findInterval on the
  # sorted values gives count(<= t) in O(m log n) instead of m separate O(n)
  # comparisons, and the running sum of the sorted values then gives
  # sum(v > t) = sum_total - cumsum[idx] with no per-threshold rescan.
  o <- order(values)
  v_sorted <- values[o]
  csum <- cumsum(v_sorted)
  sum_total <- if (n_total > 0L) csum[n_total] else 0

  idx <- findInterval(thresholds, v_sorted)
  counts_above <- n_total - idx
  volumes <- counts_above / hv@PointDensity
  probability_sums <- sum_total - c(0, csum)[idx + 1L]

  volume_quantiles <- volumes / max(volumes)
  probability_quantiles <- probability_sums / max(probability_sums)

  stats <- data.frame(threshold=thresholds,volume=volumes,volume_quantile=volume_quantiles,probability_quantile=probability_quantiles)
  row.names(stats) <- NULL

  # Materialise one sub-hypervolume at index i. Split out so the
  # quantile.requested path can build only the one it returns.
  make_hv <- function(i) {
    t_i <- thresholds[i]
    keep <- values > t_i
    hv_new <- hv
    hv_new@RandomPoints <- hv@RandomPoints[keep,,drop=FALSE]
    hv_new@Volume <- volumes[i]
    if (uniform.density==TRUE)
    {
      # Preserves a defect in the original: `length(hv_new@ValueAtRandomPoints)`
      # was evaluated when hv_new@ValueAtRandomPoints was still the parent's
      # values, so this vector has length n_total, not sum(keep). Pinned by
      # test-05's "uniform.density controls whether values are flattened" test
      # (item 0 in TESTING.md's known-defects list). Fix separately.
      hv_new@ValueAtRandomPoints <- rep(1, n_total)
    }
    else
    {
      hv_new@ValueAtRandomPoints <- values[keep]
    }
    hv_new@Name <- sprintf("%s: t: %.3f - qv: %.3f - qp: %.3f", hv@Name, stats$threshold[i], stats$volume_quantile[i], stats$probability_quantile[i])
    hv_new
  }

  if (plot==TRUE & length(thresholds) > 1)
  {
    plot(volume_quantiles~thresholds,type='s',xlab="Threshold value",ylab="Quantile",lwd=2,yaxt='n',cex.axis=0.5,ylim=c(0,1),col='purple')
    lines(probability_quantiles~thresholds,type='s',lwd=2,yaxt='n',cex.axis=0.5,ylim=c(0,1),col='blue')
    axis(side=2,at=seq(0,1,by=0.1),cex.axis=0.5)
    axis(side=3,at=ceiling(seq(thresholds[1],thresholds[length(thresholds)],length.out=5)),cex.axis=0.5)
    mtext(side=3,line=3,"Index")
  }

  if (is.null(quantile.requested) & length(thresholds) > 1)
  {
    if (verbose==TRUE)
    {
      pb <- progress_bar$new(total = length(thresholds))
      pb$tick(0)
    }
    result <- vector(mode="list",length=length(thresholds))
    for (i in seq_along(thresholds))
    {
      if (verbose==TRUE)
      {
        pb$tick()
      }
      result[[i]] <- make_hv(i)
    }
    if (verbose==TRUE)
    {
      pb$terminate()
    }
    names(result) <- thresholds

    if (plot==TRUE)
    {
      legend('topright',c("Volume","Probability"),lty=1,col=c("purple","blue"),lwd=2,cex=0.75,bty='n')
    }

    return(list(
      HypervolumesThresholded=new("HypervolumeList",HVList=result),
      Statistics=stats
    ))
  }
  else
  {
    # choose the first hypervolume that fits the criterion
    if (quantile.requested.type=="volume")
    {
      quantile.obtained.index <- head(which(stats$volume_quantile <= quantile.requested),1)
    }
    else if (quantile.requested.type=="probability")
    {
      quantile.obtained.index <- head(which(stats$probability_quantile <= quantile.requested),1)
    }
    else
    {
      stop('Unrecognized value for argument quantile.requested.type.')
    }
    quantile.obtained = stats[quantile.obtained.index,ifelse(quantile.requested.type=="probability","probability_quantile","volume_quantile")]
    threshold.obtained = stats[quantile.obtained.index,"threshold"]

    if (plot==TRUE & length(thresholds) > 1)
    {
      abline(h=quantile.obtained, col='green',lwd=2)
      abline(v=threshold.obtained)

      legend('topright',c("Volume","Probability"),lty=1,col=c("purple","blue"),lwd=2,cex=0.75,bty='n')
    }

    if (verbose==TRUE)
    {
      cat(sprintf('Requested %s quantile %f, obtained %f - setting threshold value %f.\n For a closer match, you can increase num.thresholds in hypervolume_threshold.\n', quantile.requested.type, quantile.requested, quantile.obtained, threshold.obtained))
    }

    return(list(
      HypervolumesThresholded=make_hv(quantile.obtained.index),
      Statistics=stats
    ))
  }
}
