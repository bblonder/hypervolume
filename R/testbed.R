
estimate_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

estimate_local_density <- function(samples, np, box.widths)
{
  if(!is.list(box.widths))
  {
    stop('box.widths must be of class list')
  }
  
  d = ncol(samples)
  # don't try to do more samples than than exist
  np = min(np, nrow(samples))
  
  for (i in 1:length(box.widths))
  {
    if (length(box.widths[[i]])==1)
    {
      box.widths[[i]] = rep(box.widths[[i]], d)
    }
    if (length(box.widths[[i]])!=d)
    {
      stop('Box widths must have either dimensionality 1 or that of the data.')
    }
  }
  
  box_volumes = sapply(box.widths, function(box.width) { prod(2*box.width) }) 
  
  points_at = samples[sample(1:nrow(samples),np),]
  
  # iterate over each box.width
  counts = evalfrectangular_multiple(samples, box.widths, points_at)
  
  point_density = matrix(NA,nrow=nrow(counts),ncol=ncol(counts))
  for (i in 1:length(box.widths))
  {
    point_density[,i] = counts[,i] / box_volumes[i]
  }
  
  return(point_density)
  
  pairs(points_at, col=rainbow(15)[cut(counts,10)])
  
  
  
  modal_density = estimate_mode(point_density)
  
  plot(density(point_density))
  abline(v=modal_density,col='red')
  
  return(list(modal_density=modal_density,points_at=points_at,counts=counts))
}

importance_sample_integral <- function(model, data, N.samples)
{
  d = ncol(data)
  
  mean = apply(data, 2, mean, na.rm=T)
  sd = apply(data, 2, function(x) { quantile(x, 0.75,na.rm=T) - quantile(x,0.25,na.rm=T) })
  
  ys = rmvnorm(n = N.samples, mean=start, sigma=sd*diag(d))
  
  yc = ys
  r = rep(NA, d)
  for (i in 1:d)
  {
    qlo = quantile(ys[,i],0.25)
    qhi = quantile(ys[,i], 0.75)
    
    yc = yc[yc[,i] > qlo & yc[,i] < qhi,]
    
    r[i] = qhi - qlo
  }
  
  vals = predict(model, newdata=yc)
  
  integrand = sum(vals)/nrow(ys)
  
  normalization = 1/prod(r)
  
  volume = integrand * normalization
  
  return(volume)
}