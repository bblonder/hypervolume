rgb_2_rgba <- function(color, alpha)
{
  vals <- col2rgb(color,alpha=FALSE)

  result <- rgb(vals[1]/255, vals[2]/255, vals[3]/255, alpha)
  
  return(result)
}

rgb_2_set_hsv <- function(color, h=NULL, s=NULL, v=NULL)
{
 if (is.null(color))
 {
   color <- 'red'
 }
  
  vals <- col2rgb(color)
  vals <- rgb2hsv(vals[1],vals[2],vals[3])
  
  if (!is.null(h))
  {
    vals[1] = h
  }    
  if (!is.null(s))
  {
    vals[2] = s
  }
  if (!is.null(v))
  {
    vals[3] = v
  }
  result <- hsv(vals[1],vals[2],vals[3])

  return(result)
}
