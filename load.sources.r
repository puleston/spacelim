load.sources <- function(x)
 {
  #preliminary identification of functions and function locations
  prefix <- SS
  source.file <- paste(prefix, "/", x, sep="")
  function.files <- source(source.file)$value
  function.paths <- paste(prefix, "/",function.files, sep="")

  #loop to load all the functions
  for(i in 1:length(function.paths)) {source(function.paths[i])}
}