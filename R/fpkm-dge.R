#' @title Using \code{limma_dge} for \code{raw} and \code{fpkm} counts
#' 
#' @description \code{limma_dge} uses the \code{limma} package for DGE 
#' analysis on \code{fpkm_counts} and \code{raw_counts} object, and optionally 
#' also directly on an \code{ExpressionSet} object.
#' 
#' @section Note: The term differential gene expression or \code{DGE} is not 
#' used in a restrictive manner and applies to genomic features in general, 
#' i.e., genes, transcripts, exons etc.
#' 
#' @param x An object of class \code{fpkm_counts}, \code{raw_counts} or 
#' \code{ExpressionSet}.
#' @param ... Additional arguments.
#' @param design A design matrix. See \code{\link{construct_design}}.
#' @param contrast Contrast to compute DGE for. See 
#' \code{\link{construct_contrasts}}.
#' @param voom Logical. If \code{TRUE}, applies \code{voom} transformation. 
#' Default is \code{FALSE}. Only relevant for \code{limma_dge} method.
#' 
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' 
#' @return An object of class \code{dge} consisting of top tables from 
#' running the corresponding \code{limma_dge} method.
#' @export 
#' @aliases limma_dge fpkm-dge
#' @examples
#' path = system.file("tests", package="ganalyse")
#'
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
#' fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
#' fpkm_design = construct_design(fpkm_counts, ~ 0 + condition)
#' fpkm_contrasts = construct_contrasts(
#'                      design = fpkm_design, 
#'                      treatA.vs.control = conditiontreatA-conditioncontrol, 
#'                      treatB.vs.control = conditiontreatB-conditioncontrol
#'                  )
#' # DE genes between treatA and control
#' limma_dge(fpkm_counts, design=fpkm_design, 
#'              contrast=fpkm_contrasts[, "treatA.vs.control"])
#' # DE genes between treatB vs control
#' limma_dge(fpkm_counts, design=fpkm_design, 
#'              contrast=fpkm_contrasts[, "treatB.vs.control"])
#' # DE genes between any of the treatments
#' limma_dge(fpkm_counts, design=fpkm_design, contrast=fpkm_contrasts)
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' raw_design = construct_design(raw_counts, ~ 0 + condition)
#' raw_contrasts = construct_contrasts(
#'                      design = raw_design, 
#'                      treatA.vs.control = conditiontreatA-conditioncontrol, 
#'                      treatB.vs.control = conditiontreatB-conditioncontrol
#'                  )
#' # DE genes between treatA and control
#' limma_dge(raw_counts, design=raw_design, 
#'              contrast=raw_contrasts[, "treatA.vs.control"])
#' # DE genes between treatB and control
#' limma_dge(raw_counts, design=raw_design, 
#'              contrast=raw_contrasts[, "treatB.vs.control"])
#' # DE genes between any of the treatments
#' limma_dge(raw_counts, design=raw_design, contrast=raw_contrasts)
#' 
limma_dge <- function(x, ...) {
    UseMethod("limma_dge")
}

#' @rdname limma_dge
#' @export
limma_dge.default <- function(x, ...) {
    stop("No methods available for objects of class, ", class(x)[1L])
}

#' @rdname limma_dge
#' @export
limma_dge.ExpressionSet <- function(x, design, contrast, voom=FALSE, ...) {

    voom = suppressWarnings(as.logical(voom[1L]))
    if (is.na(voom))
        stop("'voom' argument must be logical TRUE/FALSE.")
    if (missing(design)) 
        stop("design matrix is missing. Please see ?construct_design.")
    if (!is.matrix(design)) stop("design argument must be a matrix.")
    if (missing(contrast))
        stop("contrast argument is missing. Please see construct_contrasts.")
    if (!is.matrix(contrast)) contrast = as.matrix(contrast)

    if (voom) {
        # fpkm doesn't need voom (which log2 transforms)
        eset = limma::vooma(x, design)
    } else eset = x
    fit = lmFit(eset, design)
    fit = contrasts.fit(fit, contrast)
    fit = eBayes(fit, trend=TRUE)
    fit$genes = Biobase::fData(x)[["id"]]
    if (is.null(fit$genes)) {
        warning("Biobase::fData(x)$id returned NULL. Column 'id' ", 
                    "is non-existent in the input eset object.")
    }
    format = if (voom) "limma_voom_fpkm" else "limma_fpkm"
    limma_toptables_dge(fit, contrast, format)
}


#' @rdname limma_dge
#' @export
limma_dge.fpkm_counts <- function(x, design, contrast, voom=FALSE, ...) {

    x = suppressWarnings(as.eset(x))
    limma_dge(x, design, contrast, voom, ...)
}


#' @rdname limma_dge
#' @export
limma_dge.raw_counts <- function(x, design, contrast, voom=TRUE, ...) {
    
    voom = suppressWarnings(as.logical(voom[1L]))
    if (is.na(voom))
        stop("'voom' argument must be logical TRUE/FALSE.")
    if (missing(design)) 
        stop("design matrix is missing. Please see ?construct_design.")
    if (!is.matrix(design)) stop("design argument must be a matrix.")
    if (missing(contrast))
        stop("contrast argument is missing. Please see construct_contrasts.")
    if (!is.matrix(contrast)) contrast = as.matrix(contrast)

    # can do this directly with as.eset() as well, but this is fine
    dlist = as.dgelist(x)
    dlist = edgeR::calcNormFactors(dlist)

    eset = voom(dlist, design, plot=FALSE)
    fit = lmFit(eset, design)
    fit = contrasts.fit(fit, contrast)
    fit = eBayes(fit)
    limma_toptables_dge(fit, contrast, 
            if (voom) "limma_voom_raw" else "limma_raw")
}

#' @title Write top table results to file.
#' 
#' @description This is a convenience function for writing multiple top 
#' tables to files. It takes a list of top tables and corresponding file 
#' names and writes each table to that file.
#' 
#' @param x A list of top tables
#' @param filename A character vector of same length as \code{x} containing 
#' filenames corresponding to each element in \code{x}.
#' @seealso \code{\link{rnaseq}} \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}}
#' \code{\link{as.dgelist}} \code{\link{as.eset}} 
#' \code{\link{volcano_plot}}
#' 
#' @return It returns a \code{NULL} value if writing file to disk 
#' completes successfully.
#' @examples
#' path = system.file("tests", package="ganalyse")
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' raw_design = construct_design(raw_counts, ~ 0 + condition)
#' raw_contrasts = construct_contrasts(
#'                      design = raw_design, 
#'                      treatA.vs.control = conditiontreatA-conditioncontrol, 
#'                      treatB.vs.control = conditiontreatB-conditioncontrol
#'                  )
#' # DE genes between treatA and control
#' ans = edger_dge(raw_counts, design=raw_design, 
#'              contrast=raw_contrasts[, "treatA.vs.control"])
#' \dontrun{
#' write_dge(ans, "treatA_vs_control.tsv")
#' }
#' @export
write_dge <- function(x, filename) {
    if (is.data.frame(x)) x = list(x)
    if (is.atomic(filename)) filename = as.list(filename)
    stopifnot(length(x) == length(filename))
    # TODO write.table -> fwrite after v1.9.8 hits CRAN
    Map(function(dt, pth) write.table(dt, pth, quote=FALSE, sep="\t", 
            row.names=FALSE, col.names=TRUE), x, filename)
    NULL
}

# ----- internal function ----- #

limma_toptables_dge <- function(fit, contrast, format) {
    rows = nrow(fit$coefficients)
    stopifnot(is.matrix(contrast))
    ans = setDT(topTable(fit, number=rows))
    old = c("ProbeID", "ID", "P.Value", "adj.P.Val", "logFC", "AveExpr")
    new = c("id", "id", "pval", "padj", "log2fc", "ave_expr")
    idx = which(old %chin% names(ans))
    if (length(idx)) setnames(ans, old[idx], new[idx])
    setnames(ans, gsub("\\.", "_", names(ans)))
    setnames(ans, gsub("logFC", "log2fc", names(ans)))
    setattr(ans, 'class', unique(c("dge", class(ans))))
    ans[, "method" := format][]
}
