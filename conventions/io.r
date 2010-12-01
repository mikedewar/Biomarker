library(yaml)
structure <- yaml.load_file('../conventions/structure.yaml')
data.dir  <- paste(structure$root,structure$data,sep="")
cache.dir <- paste(structure$root,structure$cache,sep="")
cat("\nthe following variables have been defined:\n")
cat("\tstructure: list describing the directory structure\n")
cat("\tdata.dir: string containing the data directory\n")
cat("\tcache.dir: string containing the cache directory\n")
cat("\n")

data.load <- function(filename){
  filename = paste(structure$root,structure$data,filename,sep="")
  print(filename)
  load(filename, sys.frame(1))
}

cache.save <- function(cache_variable,filename){
    filename = paste(structure$root,structure$cache,filename,sep="")
    # look - every variable we save will be called "cache variable"
    save(cache_variable,file=filename)
}

cache.load <- function(filename){
  filename = paste(structure$root,structure$cache,filename,sep="")
  load(filename) # this will be called cache_variable
  return(cache_variable)
}
