estimate_threshold_gaussian <- function(sd.count=1, kde.bandwidth)
{
  k = length(kde.bandwidth) # dimensionality
  # standardize distance to a given number of s.d.s away
  probs <- rep(NA, length(kde.bandwidth))
  for (i in 1:length(kde.bandwidth))
  {
    probs[i] <- 1/sqrt(2*pi*(2*kde.bandwidth[i])^2) * exp(-1/2*sd.count^2)
  }
  # total probability is product of probability along each axis (no covariance)
  prob_final <- prod(probs)
  return(prob_final)
}
