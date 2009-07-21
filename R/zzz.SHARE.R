## without NAMESPACE
##.First.lib <- function(lib, pkg) {
##   library.dynam("SHARE", pkg, lib)
##}

## with NAMESPACE
.onLoad <- function(lib, pkg) {
   library.dynam("SHARE", pkg, lib)
}
