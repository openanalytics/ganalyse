#' @title Create an RNASeq experiment object
#' 
#' @description Loads all sample info related to an experiment on to an object 
#' of class \code{fpkm} or \code{raw}.
#' 
#' 
#' @param sample_info Complete path to file that contains \code{sample} names, 
#' columns that identify the \code{groups} each \code{sample} belongs to 
#' and path to the file that contains \emph{read counts} for each sample. 
#' See \code{details} section for more.
#' @param experiment Name of the experiment. Default is "example".
#' @param format Is the input data \code{fpkm} or \code{raw} counts?
#' @param verbose Logical. Default is \code{FALSE}. If \code{TRUE}, sends 
#' useful status messages to the console. 
#' 
#' @details The file provided to \code{sample_info} (say 
#' \code{annotation.txt}) \bold{must} have columns named \code{sample} 
#' (corresponding to sample names) and \code{path} (containing the file path 
#' to \code{raw} or \code{fpkm} counts). The path for each sample can either 
#' be a complete path or just the filename (with extension). If it is just a 
#' file name, it is assumed to be in the same path as \code{annotation.txt}.
#'
#' The other columns usually identify which groups and/or treatments each 
#' sample belongs to. This allows easy construction of design matrix to 
#' describe treatment conditions later on, with the help of 
#' \code{\link{construct_design}}. These columns are optional at this point 
#' and can be added manually later on as well.
#' 
#' The easiest way would be to place \code{annotation.txt} and the counts in 
#' the folder.
#'  
#' @return An object of class \code{fpkm} or \code{raw}. It 
#' inherits from \code{data.table}.
#' @aliases rnaseq
#' @seealso \code{\link{gather_counts}} \code{\link{show_counts}}
#' \code{\link{limma_dge}} \code{\link{edger_dge}} 
#' \code{\link{construct_design}} \code{\link{construct_contrasts}}
#' \code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
#' \code{\link{volcano_plot}} \code{\link{density_plot}}
#' @export 
#' @examples
#' path = system.file("tests", package="ganalyse")
#'
#' # ----- fpkm ----- # 
#' fpkm_path = file.path(path, "fpkm", "annotation.txt")
#' (fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample"))
#' class(fpkm_obj)
#'
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' (raw_obj = rnaseq(raw_path, format="raw", experiment="sample"))
#' class(raw_obj)
rnaseq <- function(sample_info, experiment="example", format=c("fpkm", "raw"), 
                        verbose=FALSE) {

    experiment = experiment[1L]
    if (!is.character(experiment) || 
        length(experiment) != 1L || 
        is.na(experiment))
        stop("experiment must be a character vector of length=1.")
    token = tokenise(sample_info, experiment, match.arg(format), 
                verbose)
    ans = rnaseq_obj(token)
    ans[]
}

# ----- Helper/Internal functions for rnaseq -----

is.fpkm <- function(x) inherits(x, 'fpkm')
is.raw <- function(x) inherits(x, 'raw')

tokenise <- function(sample_info, experiment, format, verbose) {

    tokens = structure(list(sample_info=sample_info, experiment=experiment, 
                verbose=verbose), class=paste0(format, ""))
    tokens
}

rnaseq_obj <- function(token) {
    UseMethod("rnaseq_obj")
}

rnaseq_obj.default <- function(token) {
    stop("No method available for format ", 
        gsub("(.*)$", "", class(token)))
}

rnaseq_obj.fpkm <- function(token) {

    obj = load_experiment(token[["sample_info"]], token[["experiment"]])
    data.table::setattr(obj, 'class', c("fpkm", class(obj)))
}

rnaseq_obj.raw <- function(token) {

    obj = load_experiment(token[["sample_info"]], token[["experiment"]])
    data.table::setattr(obj, 'class', c("raw", class(obj)))
}

load_experiment <- function(sample_info, experiment) {

    obj = fread(sample_info)
    invalid = setdiff(c("sample", "path"), names(obj))
    if (length(invalid)) stop("Column(s) [", paste(invalid, collapse=","), 
        "] not found in the file provided to argument sample_info.")
    only_fn = which(basename(obj[["path"]]) == obj[["path"]])
    path=NULL
    if (length(only_fn)) {
        obj[(only_fn), "path" := file.path(dirname(sample_info), path)]
    }
    invalid = which(!file.exists(obj[["path"]]))
    if (length(invalid)) stop("Files corresponding to samples [", 
            paste(obj[["sample"]][invalid], collapse=","), 
            "] do not exist. Please recheck.")
    if (anyDuplicated(obj[["sample"]]))
        stop("Sample names must be unique.")
    # expand file paths by default
    obj[, "path" := normalizePath(path)]
    obj[, "experiment" := experiment]
    setcolorder(obj, c("experiment", setdiff(names(obj), "experiment")))
}
