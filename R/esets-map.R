#' @title Detailed \code{ExpressionSet} object
#' 
#' @description \code{as.eset_map()} converts an object of class 
#' \code{fpkm_counts} or \code{raw_counts} to an \code{ExpressionSet} object, 
#' with a lot of other details that are normally required for DGE analysis, 
#' but might be necessary at times. 
#' 
#' This function is still in development, hence not exported.
#' 
#' @param x An object of class \code{fpkm_counts}, \code{raw_counts} 
#' or \code{DGEList}.
#' @param entrez_map Full path to file that contains mapping of gene names to 
#' ENTREZ ids along with corresponding gene ids . It must be a three column 
#' \code{data.table}, with column names \code{gene_name}, \code{entrez_id} 
#' and \code{gene_id} respectively. 
#' 
#' If the argument is \code{NULL} (default), the column \code{ENTREZID} in 
#' object \code{pData(featureData(<eset_obj>))} will be set to \code{NA}.
#' @param genome_annotation An object of class \code{gtf} or \code{gff}.
#' @param gene_name_col A character vector of length=1 indicating the name of 
#' the column corresponding to gene names in \code{genome_annotation}
#' @param reference_group A character vector of length=1 indicating the 
#' reference sample.
#' @param ... Additional arguments specific to methods.
#' @return An object of class \code{ExpressionSet} tuned for MAP.
#' @examples 
#' \dontrun{
#'   # my_gtf is obtained using gread::read_format.
#'   eset_obj = as.eset_map(x, genome_annotation = my_gtf, 
#'                    gene_name_col = "gene_name")
#' }
as.eset_map <- function(x, entrez_map=NULL, reference_group = "", ...) {
    UseMethod("as.eset_map")
}

#' @rdname as.eset_map
as.eset_map.default <- function(x, entrez_map=NULL, ...) {
    stop("No known method for converting input object to an Expression ", 
            "Set object")
}

#' @rdname as.eset_map
as.eset_map.fpkm_counts <- function(x, entrez_map=NULL, genome_annotation, 
                            gene_name_col, reference_group = "", ...) {
    stopifnot(is.data.table(genome_annotation))
    stopifnot(gene_name_col %in% names(genome_annotation))
    eset = as.eset(x)
  
    # pheno data
    pdata = setDT(pData(phenoData(eset)))
    setnames(pdata, c("sub_group", "sample_colour"), 
                    c("subGroup", "sampleColor"))
    pdata[, c("group", "lib.size", "norm.factors") := NULL]
    pData(phenoData(eset)) = setDF(pdata, rownames = pdata[["sample"]])

    # experiment data
    create_reference <- function(pdata, condition) {
        ExpSubGroupID = NULL
        ans = data.table(VariableID = NA_real_, # not used in our case
               var_type = "categorical", # hard-coded for now
               name = condition,
               ExpSubGroupID = unique(pdata[["subGroup"]]), 
               ReferenceFor = -1)
        ref_group = unique(pdata[["subGroup"]][grepl(reference_group, 
                            pdata[["sample"]])])
        ans[ExpSubGroupID %in% ref_group, "ReferenceFor" := 0]
        list(setDF(ans))
    }
    experiment_samples = as.list(tail(attributes(x)$names, -1L))
    setattr(experiment_samples, 'names', sampleNames(eset))
    refs = create_reference(pdata, attributes(x)$group_name)

    experiment_data = experimentData(eset)
    experiment_data@name = "<NA>" # name Investigator
    experiment_data@lab = "<NA>"
    experiment_data@contact = "<NA>"
    experiment_data@title = "<NA>" # name/title Experiment
    experiment_data@abstract = "<NA>"
    experiment_data@url = "TBD"
    experiment_data@samples = experiment_samples
    experiment_data@preprocessing = list("<NA>")
    experiment_data@other = list(references=refs)
    experimentData(eset) = experiment_data

    # feature data
    pdata = pData(featureData(eset))
    names(pdata)[names(pdata) == "agg_id"] = "gene_id"
    # if ensembl map is given, fill ENTREZID
    entrez = FALSE
    entrez_id = NULL
    if (!is.null(entrez_map) && is.data.table(entrez_map) && 
            c("entrez_id", "ensembl_id") %in% names(entrez_map)) {
        entrez = TRUE
        pdata["ENTREZID"] = entrez_map[pdata, entrez_id, 
                                on="gene_id", mult="first"]
    }

    pdata["SYMBOL"] = genome_annotation[pdata, get(gene_name_col), 
                        on="gene_id", mult="first"]
    setnames(pdata, "gene_id", "ENSEMBL")

    na_cols = c(if (!entrez) "ENTREZID", "GENEFAMILY", 
                    "INICALLS", "PROBENUM", "GENENAME")
    pdata[na_cols] = NA_character_
    pData(featureData(eset)) = pdata
    eset
}

#' @rdname as.eset_map
as.eset_map.raw_counts <- function(x, entrez_map=NULL, genome_annotation, 
                            gene_name_col, reference_group = "", ...) {
    stopifnot(is.data.table(genome_annotation))
    stopifnot("gene_id" %in% names(genome_annotation))
    stopifnot(gene_name_col %in% names(genome_annotation))
    eset = as.eset(x)

  # pheno data
  pdata = setDT(pData(phenoData(eset)))
  setnames(pdata, c("sub_group", "sample_colour"), 
                    c("subGroup", "sampleColor"))
  pdata[, c("group", "lib.size", "norm.factors") := NULL]
  pData(phenoData(eset)) = setDF(pdata, rownames = pdata$sample)

    # experiment data
    create_reference <- function(pdata, condition) {
        ExpSubGroupID = NULL
        ans = data.table(VariableID = NA_real_, # not used in our case
               var_type = "categorical", # hard-coded for now
               name = condition,
               ExpSubGroupID = unique(pdata[["subGroup"]]), 
               ReferenceFor = -1)
        ref_group = unique(pdata[["subGroup"]][grepl(reference_group, 
                                pdata[["sample"]])])
        ans[ExpSubGroupID %in% ref_group, "ReferenceFor" := 0]
        list(setDF(ans))
    }
    experiment_samples = as.list(tail(attributes(x)$names, -1L))
    setattr(experiment_samples, 'names', sampleNames(eset))
    refs = create_reference(pdata, attributes(x)$group_name)

    experiment_data = experimentData(eset)
    experiment_data@name = "<NA>" # name Investigator
    experiment_data@lab = "<NA>"
    experiment_data@contact = "<NA>"
    experiment_data@title = "<NA>" # name/title Experiment
    experiment_data@abstract = "<NA>"
    experiment_data@url = "TBD"
    experiment_data@samples = experiment_samples
    experiment_data@preprocessing = list("<NA>")
    experiment_data@other = list(references=refs)
    experimentData(eset) = experiment_data

    # feature data
    pdata = pData(featureData(eset))
    names(pdata)[names(pdata) == "agg_id"] = "gene_id"
    # if ensembl map is given, fill ENTREZID
    entrez = FALSE
    entrez_id = NULL
    if (!is.null(entrez_map) && is.data.table(entrez_map) && 
            c("entrez_id", "ensembl_id") %in% names(entrez_map)) {
        entrez = TRUE
        pdata["ENTREZID"] = entrez_map[pdata, entrez_id, 
                                on="gene_id", mult="first"]
    }

    pdata["SYMBOL"] = genome_annotation[pdata, get(gene_name_col), 
                            on="gene_id", mult="first"]
    setnames(pdata, "gene_id", "ENSEMBL")

    na_cols = c(if (!entrez) "ENTREZID", "GENEFAMILY", "INICALLS", 
                    "PROBENUM", "GENENAME")
    pdata[na_cols] = NA
    pData(featureData(eset)) = pdata
    eset
}

# ----- internal functions ----- #
