hypervolume_funnel <- function(input_path, title = NULL, func = get_volume, CI = .95, as_table = FALSE) {
  upperq = c()
  sample_mean = c()
  lowerq = c()
  n = c()
  paths = list.files(input_path)
  foreach(i = paths) %do% {
    dat <- foreach(j = list.files(file.path(input_path, i), pattern = 'resample*'), .combine = c) %do% {
      func(readRDS(file.path(input_path, i, j)))
    }
    upperq = c(upperq, quantile(dat, .5 + CI/2))
    sample_mean = c(sample_mean, mean(dat))
    lowerq = c(lowerq, quantile(dat, (1 - CI)/2))
    n = c(n, as.numeric(strsplit(i, split = " ")[[1]][3]))
  }
  dat = data.frame(upperq, sample_mean, lowerq, n)
  if (as_table) {
    return(dat)
  }
  plot = ggplot(dat, aes(x = n)) + geom_line(aes(y = lowerq)) + 
    geom_line(aes(y = sample_mean), col = 'blue') + 
    geom_line(aes(y = upperq)) + 
    ggtitle(title, paste('Confidence interval:', as.character(CI))) +
    xlab('Resample size') +
    ylab('Parameter')
  return(plot)
}