.onLoad <- function(libname, pkgname) {
    ganalyse_verbose = "FALSE"
    # global options
    opts = c(
              "ganalyse_verbose" = ganalyse_verbose
            )
    for (i in setdiff(names(opts), names(options())) ) {
        text = paste('options(', i, '=', opts[i], ')', sep="")
        eval(parse(text=text))
    }
    invisible()
}
