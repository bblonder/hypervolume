rgb2rgba <- function(colorlist, alphalist)
{
  result <- rep(NA, length(colorlist))
  for (i in 1:length(colorlist))
  {
    vals <- col2rgb(colorlist[i],alpha=FALSE)
  
    result[i] <- rgb(vals[1]/255, vals[2]/255, vals[3]/255, alphalist[i])
  }
  
  return(result)
}

rgb2rgbdark <- function(colorlist, darkfactor)
{
  if (length(colorlist) == 0)
  {
    colorlist <- 'red'
  }
  
  result <- rep(NA, length(colorlist))
  for (i in 1:length(colorlist))
  {
    vals <- col2rgb(colorlist[i],alpha=TRUE)
    
    result[i] <- rgb(vals[1]/255*darkfactor, vals[2]/255*darkfactor, vals[3]/255*darkfactor, vals[4]/255)
  }
  
  return(result)
}
