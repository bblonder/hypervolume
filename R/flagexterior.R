expectednumber <- function(density, rmin, n)
{
  pi^(n/2) * rmin^n / gamma(n/2+1) * density
}

flagexterior <- function(input, density, fuzzfactor)
{  
  distances <- as.matrix(dist(input))
  
  n <- ncol(input)
  
  criticaldistance <- fuzzfactor * (density)^(-1/n)
  
  finalstat <- rep(NA, nrow(input))
  
  lop <- rep(NA, nrow(input))
  
  for (i in 1:nrow(input))
  {
    # find all the neighbors in the critical ball
    nnids <- which(distances[i,] < criticaldistance & distances[i,] > 0)
    
    print(length(nnids))
    str(input[nnids,])
    
    vec_this <- (input[nnids,] - input[i,])
    
    lop[i] <- var(colMeans(vec_this^2))
    #lop[i] <- sqrt(sum(colMeans(vec_this)^2))
    
    # compare to expected number of points
    finalstat[i] <- length(nnids)
  }
  finalstat <- finalstat / max(finalstat)
  lop <- lop / max(lop)
  
  return(cbind(lop, finalstat))
}

#data("iris")
#td <- iris[,1:3]#matrix(runif(10000),nrow=500,ncol=2)

#hv1 <- hypervolume(td)
#fe <- flagexterior(td, density=10000, fuzzfactor=10); 
#plot(td, col=rainbow(30)[as.numeric(cut(fe[,1]^2+fe[,2]^2, 20))])