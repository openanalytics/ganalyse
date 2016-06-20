#' @title Convert objects to an \code{ExpressionSet} object
#' 
#' @description At the moment, these objects can be converted to an 
#' \code{ExpressionSet} object directly by doing \code{as.eset(x)}.
#' \itemize{
#' \item \code{fpkm_counts} -- object created by calling 
#' \code{ganalyse::gather_counts} on a \code{fpkm_format} object.
#' 
#' \item \code{raw_counts}  -- object created by calling 
#' \code{ganalyse::gather_counts} on a \code{raw_format} object. 
#'
#' \item \code{DGEList} -- object created by calling \code{edgeR::DGEList} 
#' function.
#' }
#' 
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.dgelist}}
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' 
#' @param x An object of class \code{fpkm_format}, \code{raw_format} or 
#' \code{DGEList}.
#' @param groups columns in 'x' corresponding to the group each sample 
#' belongs to.
#' @param ... Additional arguments specific to methods.
#' @return An object of class \code{ExpressionSet}.
#' @examples 
#' path = system.file("tests", package="ganalyse")
#'
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
#' fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
#' (as.eset(fpkm_counts))
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' (as.eset(raw_counts))
#' 
#' @export
as.eset <- function(x, ...) {
    UseMethod("as.eset")
}

#' @rdname as.eset
#' @export
as.eset.default <- function(x, ...) {
    stop("No known method for converting input object ", 
          "to an Expression Set object")
}

#' @rdname as.eset
#' @export
as.eset.DGEList <- function(x, groups=NULL, ...) {
    # from RNASeq v0.10, Author: Heather Turner
    if (is.null(x$genes)) {
        return(ExpressionSet(assayData = x$counts, 
                phenoData = AnnotatedDataFrame(x$samples)))
    }
    ExpressionSet(assayData = x$counts, 
                    phenoData = AnnotatedDataFrame(x$samples), 
                    featureData = AnnotatedDataFrame(x$genes))
}

#' @rdname as.eset
#' @export
as.eset.fpkm_counts <- function(x, groups=NULL, ...) {

    eset = as.eset(as.dgelist(x))
    pdata = pData(phenoData(eset))
    if (!is.null(groups)) {
        cols = setdiff(names(pdata), "group")
        pdata = cbind(pdata[, cols, drop=FALSE], 
                        data.table::setDF(x[, groups, with=FALSE]))
        setcolreorder(pdata, groups)
    } else groups = "group" # created from as.dgelist
    pdata["sample"] = sampleNames(eset)
    paste2 <- function(...) paste(..., sep=",")
    subgroups = do.call("paste2", pdata[, groups, drop=FALSE])
    pdata["sub_group"] = as.integer(factor(subgroups))
    pdata["sample_colour"] = rainbow(max(pdata$sub_group))[pdata$sub_group]
    # pdata[attributes(x)[["group_name"]]] = paste(pdata[["sub_group"]], 
    #                                         ":", subgroups, sep=" ")
    pData(phenoData(eset)) = pdata
    eset
}

#' @rdname as.eset
#' @export
as.eset.raw_counts <- function(x, groups=NULL, ...) {

    eset = as.eset(as.dgelist(x))
    pdata = pData(phenoData(eset))
    if (!is.null(groups)) {
        cols = setdiff(names(pdata), "group")
        pdata = cbind(pdata[, cols, drop=FALSE], 
                        data.table::setDF(x[, groups, with=FALSE]))
        setcolreorder(pdata, groups)
    } else groups = "group" # created from as.dgelist
    pdata["sample"] = sampleNames(eset)
    paste2 <- function(...) paste(..., sep=",")
    subgroups = do.call("paste2", pdata[, groups, drop=FALSE])
    pdata["sub_group"] = as.integer(factor(subgroups))
    pdata["sample_colour"] = rainbow(max(pdata$sub_group))[pdata$sub_group]
    # pdata[attributes(x)[["group_name"]]] = paste(pdata[["sub_group"]], 
    #                                         ":", pdata[["group"]], sep=" ")
    pData(phenoData(eset)) = pdata
    eset
}

# ----- internal functions ----- #

is.eset <- function(x) inherits(x, 'ExpressionSet')
