if (exists('doBootstrapSeqIris')==TRUE){

  require("ggplot2")
  if (!exists('doBootstrapDemoCores')) {
    doBootstrapDemoCores = 1
  }
  
  data(iris)
  hv = hypervolume(iris[,1:4])
  # Make sure to set working directory before running. "Objects" directory gets created in current working directory.
  resample_seq_path = hypervolume_resample("iris_hvs", hv, method = "bootstrap seq", n = 50, seq = c(25, 50, 75, 100, 125, 150), cores = doBootstrapDemoCores)
  hypervolume_funnel(resample_seq_path, title = "Volume of Hypervolumes at Different Resample Sizes") + ylab("Volume")
  
  # Compute non parametric confidence interval for mean
  hypervolume_funnel(resample_seq_path, title = "Mean of Sepal Length at Different Resample Sizes",
                     func = function(x) {get_centroid(x)["Sepal.Length"]}) + 
    ylab("Sepal Length")

  hypervolume_funnel(resample_seq_path, title = "Mean of Sepal Width at Different Resample Sizes",
                     func = function(x) {get_centroid(x)["Sepal.Width"]}) + 
    ylab("Sepal Width")

  hypervolume_funnel(resample_seq_path, title = "Mean of Petal Length at Different Resample Sizes",
                     func = function(x) {get_centroid(x)["Petal.Length"]}) + 
    ylab("Petal Length")
  
  hypervolume_funnel(resample_seq_path, title = "Mean of Sepal Width at Different Resample Sizes",
                     func = function(x) {get_centroid(x)["Petal.Width"]}) + 
    ylab("Petal Width")
  
  
} else {
  message('Demo does not run by default to meet CRAN runtime requirements.')
  message('This demo runs faster by multithreading.')
  message('will take approximately 30 minutes to run with 12 cores.')
  message('To run the demo, type')
  message('\tdoBootstrapSeqIris=TRUE')
  message('\tdemo(quercus)')
  message('at the R command line prompt.')
  message('Optionally type the following for set number of cores used for demo')
  message('\tdoBootstrapDemoCores = <number of cores>   (default cores = 1)')
}