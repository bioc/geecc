.onAttach <- function(lib,pkg){
		ver <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
		ver <- as.character(ver)
		packageStartupMessage(paste("geecc", ver, "loaded
"))
		}
.onUnload <- function(libpath){ library.dynam.unload("geecc", libpath) }
