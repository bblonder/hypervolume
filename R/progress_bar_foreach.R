progress_bar_foreach <- function(iterator, fun = NULL, clear = TRUE){
  # credits to:
  # https://gist.github.com/kvasilopoulos/d49499ea854541924a8a4cc43a77fed0
  pb <- progress_bar$new(total = iterator, show_after = 0, clear = clear)
  pb$tick(0)
  
  count <- 1
  function(...) {
    count <<- count + length(list(...)) - 1
    pb$update(count/iterator)
    fun(...)
    
  }
}
