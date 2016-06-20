#' @title Using \code{edger_dge} for \code{raw} counts
#' 
#' @description \code{edger_dge} uses the \code{edgeR} package to perform DGE 
#' analysis on \code{raw_counts} object.
#' 
#' @section Note: The term differential gene expression or \code{DGE} is not 
#' used in a restrictive manner and applies to genomic features in general, 
#' i.e., genes, transcripts, exons etc.
#' 
#' @param x An object of class \code{raw_counts}.
#' @param design A design matrix. See \code{\link{construct_design}}.
#' @param contrast Contrast to compute DGE for. See 
#' \code{\link{construct_contrasts}}.
#' @param ... Additional arguments.
#' 
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' 
#' @return An object of class \code{dge} consisting of top tables from 
#' running the corresponding \code{edger_dge} method.
#' 
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
#' edger_dge(raw_counts, design=raw_design, 
#'              contrast=raw_contrasts[, "treatA.vs.control"])
#' # DE genes between treatB and control
#' edger_dge(raw_counts, design=raw_design, 
#'              contrast=raw_contrasts[, "treatB.vs.control"])
#' # DE genes between any of the treatments
#' edger_dge(raw_counts, design=raw_design, contrast=raw_contrasts)
#' 
#' @export 
edger_dge <- function(x, ...) {
    UseMethod("edger_dge")
}

#' @rdname edger_dge
#' @export
edger_dge.default <- function(x, ...) {
    stop("No default method available for object of class 'edger_dge'")
}

#' @rdname edger_dge
#' @export
edger_dge.raw_counts <- function(x, design, contrast, ...) {

    # TODO add support for by='exon' in gather_counts. edgeR has a 
    # method for exons in its user manual
    if (missing(design)) 
        stop("design matrix is missing. Please see ?construct_design.")
    if (!is.matrix(design)) stop("design argument must be a matrix.")
    if (missing(contrast))
        stop("contrast argument is missing. Please see construct_contrasts.")
    if (!is.matrix(contrast)) contrast = as.matrix(contrast)

    dlist = as.dgelist(x)
    dlist = calcNormFactors(dlist)
    dlist = estimateGLMCommonDisp(dlist, design)
    dlist = estimateGLMTrendedDisp(dlist, design)
    dlist = estimateGLMTagwiseDisp(dlist, design)

    fit = glmFit(dlist, design)
    edger_toptables_dge(fit, contrast=contrast, "edger_raw")
}

# ----- internal functions ----- #

edger_toptables_dge <- function(fit, contrast, format) {
    rows = nrow(fit$coefficients)
    stopifnot(is.matrix(contrast))
    ans = glmLRT(fit, contrast=contrast)
    ans = setDT(topTags(ans, rows)$table)
    old = c("PValue", "FDR", "logFC", "logCPM", "LR")
    new = c("pval", "padj", "log2fc", "logcpm", "lr")
    idx = which(old %chin% names(ans))
    if (length(idx)) setnames(ans, old[idx], new[idx])
    setnames(ans, gsub("\\.", "_", names(ans)))
    setnames(ans, gsub("logFC", "log2fc", names(ans)))
    setattr(ans, 'class', unique(c("dge", class(ans))))
    ans[, "method" := format][]
}
