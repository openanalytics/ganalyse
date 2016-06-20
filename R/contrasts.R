#' @title Construct \code{contrasts} matrix for comparison
#' 
#' @description \code{construct_contrasts} is a wrapper function to 
#' \code{limma::makeContrasts}. It allows to specify which comparisons from 
#' the linear model fit are to be extracted.
#' 
#' @param design A design matrix. See \code{\link{construct_design}}.
#' @param ... A set of contrasts that should be constructed, preferably in 
#' the format \code{name1 = expression1, name2 = expression2, ...} as it makes 
#' it clearer to understand the comparisons when looking at the code later 
#' on. See \code{Examples} section for more. 
#' 
#' See \code{limma::makeContrasts} for other possible ways to specify 
#' contrasts.
#' 
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{write_dge}} 
#' \code{\link{as.dgelist}} \code{\link{as.eset}} \code{\link{density_plot}}
#' 
#' @return An object of class matrix with columns corresponding to 
#' contrasts
#' @export
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
#' 
construct_contrasts <- function(design, ...) {
    UseMethod("construct_contrasts")
}

#' @rdname construct_contrasts
#' @export
construct_contrasts.default <- function(design, ...) {
    stopifnot(is.matrix(design))
    limma::makeContrasts(..., levels = design)
}
