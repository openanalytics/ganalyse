% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raw-dge.R
\name{edger_dge}
\alias{edger_dge}
\alias{edger_dge.default}
\alias{edger_dge.raw_counts}
\title{Using \code{edger_dge} for \code{raw} counts}
\usage{
edger_dge(x, ...)

\method{edger_dge}{default}(x, ...)

\method{edger_dge}{raw_counts}(x, design, contrast, ...)
}
\arguments{
\item{x}{An object of class \code{raw_counts}.}

\item{...}{Additional arguments.}

\item{design}{A design matrix. See \code{\link{construct_design}}.}

\item{contrast}{Contrast to compute DGE for. See 
\code{\link{construct_contrasts}}.}
}
\value{
An object of class \code{dge} consisting of top tables from 
running the corresponding \code{edger_dge} method.
}
\description{
\code{edger_dge} uses the \code{edgeR} package to perform DGE 
analysis on \code{raw_counts} object.
}
\section{Note}{
 The term differential gene expression or \code{DGE} is not 
used in a restrictive manner and applies to genomic features in general, 
i.e., genes, transcripts, exons etc.
}
\examples{
path = system.file("tests", package="ganalyse")

# ----- raw ----- # 
raw_path = file.path(path, "raw", "annotation.txt")
raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
raw_design = construct_design(raw_counts, ~ 0 + condition)
raw_contrasts = construct_contrasts(
                     design = raw_design, 
                     treatA.vs.control = conditiontreatA-conditioncontrol, 
                     treatB.vs.control = conditiontreatB-conditioncontrol
                 )
# DE genes between treatA and control
edger_dge(raw_counts, design=raw_design, 
             contrast=raw_contrasts[, "treatA.vs.control"])
# DE genes between treatB and control
edger_dge(raw_counts, design=raw_design, 
             contrast=raw_contrasts[, "treatB.vs.control"])
# DE genes between any of the treatments
edger_dge(raw_counts, design=raw_design, contrast=raw_contrasts)

}
\seealso{
\code{\link{rnaseq}}, \code{\link{gather_counts}} 
\code{\link{show_counts}} \code{\link{limma_dge}} 
\code{\link{construct_design}} \code{\link{construct_contrasts}} 
\code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
\code{\link{volcano_plot}} \code{\link{density_plot}}
}

