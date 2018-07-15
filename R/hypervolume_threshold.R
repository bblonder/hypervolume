hypervolume_threshold <- function(hv, thresholds=NULL, num.thresholds=20, quantile.requested=NULL,quantile.requested.type="volume",uniform.density=TRUE, plot=TRUE, verbose=TRUE)
{
  if (is.null(thresholds))
  {
    thresholds <- seq(min(hv@ValueAtRandomPoints), max(hv@ValueAtRandomPoints), length.out=num.thresholds)
  }
  quantiles <- rep(NA, length(thresholds))
  
  result <- vector(mode="list",length=length(thresholds))
  volume_quantiles <- rep(NA, length(thresholds))
  probability_sums <- rep(NA, length(thresholds))

  if (verbose==TRUE)
  {
    pb <- progress_bar$new(total = length(thresholds))
    pb$tick(0)
  }
  for (i in 1:length(thresholds))
  {
    if (verbose==TRUE)
    {
      pb$tick()
    }
    
    hv_new <- hv
    hv_new@RandomPoints <- hv@RandomPoints[hv@ValueAtRandomPoints > thresholds[i],,drop=F]
    hv_new@Volume <- length(which(hv@ValueAtRandomPoints > thresholds[i])) / hv@PointDensity
    if (uniform.density==TRUE)
    {
      hv_new@ValueAtRandomPoints <- rep(1, length(hv_new@ValueAtRandomPoints))
    }
    else
    {
      hv_new@ValueAtRandomPoints <- hv@ValueAtRandomPoints[hv@ValueAtRandomPoints > thresholds[i]]
    }
    result[[i]] <- hv_new
    
    # also calculate probability density enclosed
    probability_sums[i] <- sum(hv@ValueAtRandomPoints[hv@ValueAtRandomPoints > thresholds[i] ])
  }
  if (verbose==TRUE)
  {
    pb$terminate()
  }
  names(result) <- thresholds
  
  volumes <- sapply(result, function(x) {x@Volume})	
  volume_quantiles <- volumes/max(volumes)
  probability_quantiles <- probability_sums / max(probability_sums)
  
 
  
  stats <- data.frame(threshold=thresholds,volume=volumes,volume_quantile=volume_quantiles,probability_quantile=probability_quantiles)
  row.names(stats) <- NULL
  
  for (i in 1:length(result))
  {
    result[[i]]@Name <- sprintf("%s: t: %.3f - qv: %.3f - qp: %.3f", result[[i]]@Name, stats$threshold[i], stats$volume_quantile[i], stats$probability_quantile[i])
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
      HypervolumesThresholded=result[[quantile.obtained.index]],
      Statistics=stats
    ))
  }
}