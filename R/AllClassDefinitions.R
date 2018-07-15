setClass("Hypervolume", slots=c(
    Name="character",
    Method="character",
    Data="matrix",
    Dimensionality="numeric",
    Volume="numeric",
    PointDensity="numeric",
    Parameters="list",
    RandomPoints="matrix",
    ValueAtRandomPoints="numeric"
    ))

setClass("HypervolumeList", slots=c(
    HVList="list"
  ))


summary.Hypervolume <- function(object, ...)
{
  cat(sprintf("***** Object of class %s *****\n",class(object)))
  cat(sprintf("Name: %s\n",object@Name))
  cat(sprintf("Method: %s\n",object@Method))
  cat(sprintf("Number of data points (after weighting): %d\n",ifelse(all(is.nan(object@Data)), 0, nrow(object@Data))))
  cat(sprintf("Dimensionality: %d\n",object@Dimensionality))
  cat(sprintf("Volume: %f\n",object@Volume))
  cat(sprintf("Random point density: %f\n",object@PointDensity))
  cat(sprintf("Number of random points: %d\n",nrow(object@RandomPoints)))
  
  varp <- object@ValueAtRandomPoints
  if (length(varp)==0)
  {
    varp <- NA
  }
  cat(sprintf("Random point values:\n\tmin: %.3f\n\tmean: %.3f\n\tmedian: %.3f\n\tmax:%.3f\n",
                min(varp),
                mean(varp),
                median(varp),
                max(varp)))
  cat(sprintf("Parameters:\n"))
  if (length(object@Parameters) > 0)
  {
    lapply(1:length(object@Parameters), function(x) {
      cat(sprintf("\t%s: %s\n",names(object@Parameters)[x], paste(format(object@Parameters[[x]]),collapse=" ")))
    })
  }
  else
  {
    cat('\t(none)\n')
  }
  
  invisible()
}

summary.HypervolumeList <- function(object, ...)
{
  cat(sprintf("HypervolumeList with %d elements:\n\n", length(object@HVList)))
  
  if (length(object@HVList)>0)
  {
    for (i in 1:length(object@HVList))
    {
      whichhv <- object@HVList[[i]]
      cat(sprintf("\n$`%s` (%d / %d) \n", names(object@HVList)[[i]],i, length(object@HVList)))
      summary.Hypervolume(whichhv)
    }
  }
}

print.Hypervolume <- function(x, ...) {summary.Hypervolume(x)}
print.HypervolumeList <- function(x, ...) {summary.HypervolumeList(x)}

show.Hypervolume <- function(object) {summary.Hypervolume(object)}
show.HypervolumeList <- function(object) {summary.HypervolumeList(object)}

setMethod("show","Hypervolume", function(object) {summary.Hypervolume(object)})
setMethod("show","HypervolumeList", function(object) {summary.HypervolumeList(object)})



`[[.HypervolumeList` <- function(x, i, ...) {
  if (length(i)==1)
  {
    return(x@HVList[i, ...][[1]])
  }
  else
  {
    xss <- x
    xss@HVList <- xss@HVList[i, ...]
    return(xss)
  }
}

`[[<-.HypervolumeList` <- function (x, i, ..., value) 
{
  x@HVList[c(i,...)] <- value
}
