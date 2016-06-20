#' @title Construct design matrix
#' 
#' @description \code{construct_design} is a helper/wrapper function for easy 
#' construction of design matrices directly from an object of class 
#' raw/raw_counts or fpkm/fpkm_counts. 
#'
#' Please refer to the \code{Details} section for more.
#' 
#' Note that you can construct the design matrix directly (usually using 
#' \code{model.matrix}) without using \code{construct_design}.
#' 
#' @details To describe treatment conditions, a design matrix is required. It 
#' is constructed by providing a \code{formula} using the columns in the 
#' input raw/raw_counts or fpkm/fpkm_counts object that correspond to the 
#' groups and/or treatment each sample belongs to.
#' 
#' If the columns do not already exist in the input object, it should be added 
#' first before calling \code{construct_design}.
#' 
#' Please ensure that columns that are categorical are of type \code{factor} 
#' with appropriate reference level (see \code{\link{relevel}}). See 
#' \code{Examples} for more.
#' 
#' @param x A \code{fpkm_counts} or \code{raw_counts} object.
#' @param formula Formula to construct the \code{model matrix} with.
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_contrasts}} \code{\link{write_dge}} 
#' \code{\link{as.dgelist}} \code{\link{as.eset}} 
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' @return An object of class \code{matrix} is returned.
#' @examples
#' path = system.file("tests", package="ganalyse")
#'
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
#' fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
#' (fpkm_design = construct_design(fpkm_counts, ~ 0 + condition))
#'
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' (raw_design = construct_design(raw_counts, ~ 0 + condition))
#' @export
construct_design <- function(x, formula) {
    if (!inherits(x, "raw") && !inherits(x, "fpkm"))
        stop("x must be a raw/raw_counts or fpkm/fpkm_counts object")
    if (is.character(formula)) 
        formula = as.formula(formula)
    if (!is.formula(formula)) 
        stop("Argument formula is not of type formula or character")
    mm = model.matrix(formula, x)
    if (length(idx <- which(colnames(mm) == "(Intercept)"))) {
        colnames(mm)[idx] = "Intercept"
    }
    colnames(mm) = make.names(colnames(mm)) # avoid limma::makeContrasts warn
    rownames(mm) = x[["sample"]]
    # this_samples = x[, .(sample, rep = paste(do.call(paste, c(.BY, sep=".")), 
    #                         seq_len(.N), sep="-")), by=c(all.vars(formula))]
    # rownames(mm) = this_samples[x, rep, on="sample"]
    mm
}

# ----- internal functions ----- #

is.formula <- function(x) inherits(x, 'formula')
is.design <- function(x) inherits(x, 'design')
