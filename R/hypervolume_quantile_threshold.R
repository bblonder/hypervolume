hypervolume_quantile_threshold <- function(hv, thresholds=NULL, num.thresholds=100, plot=TRUE, verbose=TRUE, quantile.requested=NULL)
{
  if (is.null(thresholds))
  {
    thresholds <- seq(min(hv@ProbabilityDensityAtRandomUniformPoints), max(hv@ProbabilityDensityAtRandomUniformPoints), length.out=num.thresholds)
  }
  quantiles <- rep(NA, length(thresholds))

  if(all(thresholds==thresholds[1]))
  {
    warning('All probability values in hypervolume are equivalent - no threshold range available.')
  }
  
  result <- vector(mode="list",length=length(thresholds))
  for (i in 1:length(thresholds))
  {
    if (verbose==TRUE)
    {
      cat(sprintf('%.3f ',i/length(thresholds)))
      if (i%%20==0)
      {
        cat('\n')
      }
    }
    hv_new <- hv
    hv_new@RandomUniformPointsThresholded <- hv@RandomUniformPointsThresholded[hv@ProbabilityDensityAtRandomUniformPoints > thresholds[i],,drop=F]
    hv_new@Volume <- length(which(hv@ProbabilityDensityAtRandomUniformPoints > thresholds[i])) / hv@PointDensity
    hv_new@ProbabilityDensityAtRandomUniformPoints <- normalize_probability(rep(1, length(hv_new@ProbabilityDensityAtRandomUniformPoints)),hv_new@PointDensity)
    result[[i]] <- hv_new
  }
  names(result) <- thresholds
  
  volumes <- sapply(result, function(x) {x@Volume})	
  quantiles <- (volumes/max(volumes))
  
  if (plot==TRUE)
  {
    plot(quantiles~thresholds,type='s',xlab="Threshold value",ylab="Volume quantile",lwd=2,yaxt='n',cex.axis=0.5,ylim=c(0,1))
    axis(side=2,at=seq(0,1,by=0.1),cex.axis=0.5)
    axis(side=3,at=seq(thresholds[1],thresholds[length(thresholds)],length.out=10),labels=seq(1,length(thresholds),length.out=10),cex.axis=0.5)
    mtext(side=3,line=3,"Index")
    abline(h=0.95,col='red',lwd=2)
    abline(h=0.50,col='orange',lwd=2)
    if(!is.null(quantile.requested))
    {
      abline(h=quantile.requested, col='green',lwd=2)
      legend('topright',c("50% quantile","95% quantile","Requested quantile"),col=c("orange","red","green"),lty=1,lwd=2,cex=0.75)
    }
    else
    {
      legend('topright',c("50% quantile","95% quantile"),col=c("orange","red"),lty=1,lwd=2,cex=0.75)
    }

  }
  
  stats <- data.frame(threshold=thresholds,quantile=quantiles,volume=volumes)
  row.names(stats) <- NULL
  
  for (i in 1:length(result))
  {
    result[[i]]@Name <- sprintf("%s: t: %.3f - q: %.3f", result[[i]]@Name, stats$threshold[i], stats$quantile[i])
  }
  
  if (is.null(quantile.requested))
  {
    return(list(
      HypervolumesThresholded=new("HypervolumeList",HVList=result),
      Statistics=stats
    ))
  }
  else
  {
    # choose the first hypervolume that fits the criterion
    quantile.obtained.index <- head(which(stats$quantile <= quantile.requested),1)

    if (verbose==TRUE)
    {
      cat(sprintf('Requested quantile %f, obtained %f.\n For a closer match, you may be able to increase num.thresholds.\n', quantile.requested, quantiles[quantile.obtained.index]))
    }
    
    return(result[[quantile.obtained.index]])
  }
}