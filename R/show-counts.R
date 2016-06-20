#' @title Display all counts in a rectangular format
#' 
#' @description This function is for interactive purposes to inspect the 
#' read counts. The first column corresponds to feature id and all other 
#' columns correspond to read counts for each sample.
#' 
#' @param x An object of class \code{fpkm_counts} or \code{raw_counts}. See 
#' \code{?gather_counts}.
#' @return An object of class \code{data.table} containing the feature and 
#' associated read counts for each sample in a tabular format.
#' @seealso \code{\link{rnaseq}} \code{\link{limma_dge}} 
#' \code{\link{edger_dge}} \code{\link{as.eset}} \code{\link{gather_counts}}
#' \code{\link{construct_design}} \code{\link{construct_contrasts}}
#' \code{\link{volcano_plot}} \code{\link{write_dge}} 
#' \code{\link{density_plot}}
#' @export
#' @examples 
#' path = system.file("tests", package="ganalyse")
#' 
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
#' fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
#' head(show_counts(fpkm_counts))
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' head(show_counts(raw_counts))
show_counts <- function(x) {
    UseMethod("show_counts")
}

#' @rdname show_counts
#' @export
show_counts.default <- function(x) {
    stop("'show_counts' expects input to be xects of class 'fpkm_counts ", 
            "or 'raw_counts'. See ?gather_counts for more.")
}

#' @rdname show_counts
#' @export
show_counts.fpkm_counts <- function(x) {

    check_counts(x, "data.table")
    ans = bind_counts(x)
    ans = dcast(ans, id ~ sample, value.var="fpkm", fill=0)
    setcolorder(ans, c("id", x[["sample"]]))[]
}

#' @rdname show_counts
#' @export
show_counts.raw_counts <- function(x) {

    check_counts(x, "data.table")
    ans = bind_counts(x)
    ans = dcast(ans, id ~ sample, value.var="raw", fill=0)
    setcolorder(ans, c("id", x[["sample"]]))[]
}

# ----- internal functions ----- #

check_counts <- function(x, class) {
    if (!"counts" %in% names(x)) {
        stop("Couldn't find 'counts' col. Please run gather_counts() first.")
    }
    first_class <- function(x) class(x)[1L]
    classes = lapply(x[["counts"]], first_class)
    if (length(invalid <- which(classes != class))) {
        stop("Rows corresponding to samples [", 
            paste(x[["sample"]][invalid], collapse=","), 
            "] do not have have a ", class, " object for 'counts'.")
    }
}

bind_counts <- function(x) {
    ll = x[["counts"]]
    names(ll) = x[["sample"]]
    rbindlist(ll, idcol="sample")
}
