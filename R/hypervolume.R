hypervolume <- function(data, method="gaussian", ...)
{
  data <- as.matrix(data)
  
  if (nrow(data) == 0)
  {
    stop('Hypervolume cannot be computed with empty input data.')
  }
  
  if (any(is.na(data)))
  {
    stop('Hypervolume cannot be computed with missing values. Reduce dimensionality or remove observations with missing values.')
  }  
  
  ### CHECK FOR WARNINGS
  corMatrix = cor(data)
  corMatrix[lower.tri(corMatrix,diag=TRUE)] <- NA
  correlationThreshold = 0.8
  corInd = which(abs(corMatrix) > correlationThreshold,arr.ind=TRUE)
  finalstring = ""
  if (length(corInd) >= 1)
  {
    for (i in 1:nrow(corInd))
    {
      finalstring <- c(finalstring, sprintf('\n\tDimensions %s and %s are highly correlated (r=%.2f)', 
                                            names(data)[corInd[i,1]], names(data)[corInd[i,2]], corMatrix[corInd[i,1], corInd[i,2]]))
    }
    finalstring <- c(finalstring, "\nConsider removing some axes.")
    warning(finalstring)
  }
  
  if (log(nrow(data)) <= ncol(data))
  {
    warning(sprintf('Log number of observations (%.2f) is less than or equal to the number of dimensions (%d).\nYou may not have enough data to accurately estimate a hypervolume with this dimensionality.\nConsider reducing the dimensionality of the analysis.',
                    log(nrow(data)), ncol(data)))
  }
  
  if (nrow(data) > 2)
  {
    sds = apply(data, 2, sd)
    sdsds <- sd(sds)
    if (!any(is.na(sdsds)))
    {
      if (sd(sds) > 5)
      {
        finalstring = 'Some dimensions have much more higher standard deviations than others:\n'
        for (i in 1:length(sds))
        {
          finalstring <- c(finalstring, sprintf('\t%s %.2f\n',names(sds)[i], sds[i]))
        }
        warning(c(finalstring, 'Consider rescaling axes before analysis.'))
      }
    }
  }
  else
  {
    warning('Data contain two or fewer observations - cannot perform scaling checks. ')
  }
  
  if (method=="box")
  {
    hv = (hypervolume_box(data=data, ...))
  }
  else if (method=="svm")
  {
    hv = (hypervolume_svm(data=data, ...))
  }
  else if (method=="gaussian")
  {
    hv = (hypervolume_gaussian(data=data, ...))
  }
  else
  {
    stop(sprintf("Method %s not recognized.",method))
  }
  
  npmin = ceiling(10^sqrt(hv@Dimensionality))
  if (nrow(hv@RandomPoints) < npmin)
  {
    warning(sprintf("Hypervolume is represented by a low number of random points (%d) - suggested minimum %d.\nConsider increasing samples.per.point to improve accuracy.",nrow(hv@RandomPoints),npmin))
  }
  
  return(hv)
}