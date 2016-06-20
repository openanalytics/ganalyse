#' ganalyse: Easy Analysis of RNASeq DE
#' 
#' \code{ganalyse} is an R package that aims to make analysis on genomics 
#' data easy. It is implemented only for RNA-Seq data at the moment.
#' 
#' The current organisation is as follows:
#' 
#' \itemize{
#' \item Provide \code{experiment}, \code{sample_info} and \code{format} to 
#' the function \code{\link{rnaseq}} to generate an object of class 
#' \code{fpkm} or \code{raw}.
#' 
#' \item Collect counts using \code{\link{gather_counts}} to load the 
#' corresponding count data for the samples. It returns an object of class 
#' \code{fpkm_counts} or \code{raw_counts} respectively.
#' 
#' \item Construct \code{design matrix} and \code{contrasts}. See 
#' \code{construct_design} and \code{construct_contrasts}.
#' 
#' \item Perform different types on analyses (only Differential Gene 
#' Expression (DGE) analysis is supported currently) on the 
#' \code{fpkm_counts} or \code{raw_counts} object.
#' 
#' For \code{fpkm_counts} object, use \code{limma_dge} to perform DGE 
#' analysis. Similarly, for \code{raw_counts} object, use \code{limma_dge} or 
#' \code{edger_dge} methods.
#' }
#' 
#' @section Note: The term differential gene expression or \code{DGE} is not 
#' used in a restrictive manner and applies to genomic features in general, 
#' i.e., genes, transcripts, exons etc.
#' 
#' @seealso \code{\link{rnaseq}} \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{construct_design}} 
#' \code{\link{construct_contrasts}} \code{\link{limma_dge}} 
#' \code{\link{edger_dge}} \code{\link{as.eset}}
#' @docType package
#' @name ganalyse
#' @import tools
#' @importFrom grDevices rainbow
#' @importFrom stats model.matrix as.formula
#' @importFrom utils packageVersion head tail capture.output write.table
#' @import methods
#' @import data.table
#' @import limma
#' @importFrom tools file_ext
#' @importFrom edgeR DGEList calcNormFactors estimateGLMCommonDisp glmLRT
#' @importFrom edgeR estimateGLMTrendedDisp estimateGLMTagwiseDisp glmFit
#' @importFrom edgeR topTags
#' @importFrom Biobase ExpressionSet AnnotatedDataFrame fData pData 
#' @importFrom Biobase featureData phenoData pData<- phenoData<- 
#' @importFrom Biobase featureData<- sampleNames experimentData 
#' @importFrom Biobase experimentData<-
NULL
