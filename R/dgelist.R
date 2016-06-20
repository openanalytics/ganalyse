#' @title Convert to a \code{DGEList} object
#' 
#' @description At the moment, these objects can be converted to a 
#' \code{DGEList} object directly by doing \code{as.dgelist(x)}.
#' \itemize{
#' \item \code{fpkm_counts} -- object created by calling 
#' \code{ganalyse::gather_counts} on a \code{fpkm} object.
#' 
#' \item \code{raw_counts}  -- object created by calling 
#' \code{ganalyse::gather_counts} on a \code{raw} object. 
#' }
#' 
#' @details The \code{as.dgelist} method is a simple wrapper for converting 
#' \code{raw_counts} and \code{fpkm_counts} objects to \code{DGEList} object. 
#' No additional transformations or filtering is done, i.e., it assumes that 
#' the counts are already filtered and (log) transformed as required. Use the 
#' arguments in \code{gather_counts()} method if necessary to prepare the 
#' counts as necessary.
#' 
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.eset}}
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' 
#' @param x An object of class \code{fpkm_counts} or \code{raw_counts}.
#' @param group column name in 'x' corresponding to the group each sample 
#' belongs to.
#' @param ... Additional arguments specific to method.
#' @return An object of class \code{DGEList}.
#' @examples
#' path = system.file("tests", package="ganalyse")
#'
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
#' fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
#' (as.dgelist(fpkm_counts))
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' (as.dgelist(raw_counts))
#' 
#' @export
as.dgelist <- function(x, ...) {
    UseMethod("as.dgelist")
}

#' @rdname as.dgelist
#' @export
as.dgelist.default <- function(x, ...) {
    stop("as.dgelist method is only implemented for ", 
        "objects of class 'fpkm_counts' and 'raw_counts'")
}

#' @rdname as.dgelist
#' @export
as.dgelist.fpkm_counts <- function(x, group=NULL, ...) {

    counts = show_counts(x)
    stopifnot("id" %in% names(counts))
    counts_mat = as.matrix(counts[, !"id", with=FALSE])
    rownames(counts_mat) = counts[["id"]]
    edgeR::DGEList(counts = counts_mat, 
            genes  = data.frame(id = counts[["id"]], 
                        row.names = counts[["id"]], 
                        stringsAsFactors=FALSE), 
            group  = if (is.null(group)) NULL else x[[group]])
}

#' @rdname as.dgelist
#' @export
as.dgelist.raw_counts <- function(x, group=NULL, ...) {
    
    counts = show_counts(x)
    stopifnot("id" %in% names(counts))
    counts_mat = as.matrix(counts[, !"id", with=FALSE])
    rownames(counts_mat) = counts[["id"]]
    edgeR::DGEList(counts = counts_mat, 
            genes  = data.frame(id = counts[["id"]], 
                        row.names = counts[["id"]], 
                        stringsAsFactors=FALSE), 
            group  = if (is.null(group)) NULL else x[[group]])
}

# ----- Internal functions for dgelist ----- #

setcolreorder <- function(x, names) {
    rest = setdiff(names(x), names)
    setcolorder(x, c(names, rest))
}

is.dgelist <- function(x) inherits(x, 'DGEList')
