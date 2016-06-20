.onAttach <- function(libname, pkgname) {
    # Runs when attached to search() path such as by library() or require()
    if (interactive()) {
        packageStartupMessage('v', as.character(packageVersion("ganalyse")), 
        ', type vignette("ganalyse-vignette", package="ganalyse") to get started.')
    }
}
