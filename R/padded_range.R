padded_range <- function(data, multiply.interval.amount=0, add.amount=0)
{
  if (length(multiply.interval.amount) != ncol(data))
  {
    multiply.interval.amount <- rep(multiply.interval.amount, ncol(data))
  }
  if (length(add.amount) != ncol(data))
  {
    add.amount <- rep(add.amount, ncol(data))
  }
  
  result_final = matrix(NA, nrow=2,ncol=ncol(data))
  for (i in 1:ncol(data))
  {
    range_this = range(data[,i])
    delta = diff(range_this)
    
    result_final[1,i] = range_this[1] - delta*multiply.interval.amount[i] - add.amount[i]
    result_final[2,i] = range_this[2] + delta*multiply.interval.amount[i] + add.amount[i]
  }
  
  if (is.null(dimnames(data)[[2]]))
  {
    dimnames(data) <- list(NULL,paste("X",1:ncol(data),sep=""))
  }
  dimnames(result_final) <- list(c("low","high"),dimnames(data)[[2]])
  
  return(result_final)
}
