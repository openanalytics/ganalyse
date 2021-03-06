% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contrasts.R
\name{construct_contrasts}
\alias{construct_contrasts}
\alias{construct_contrasts.default}
\title{Construct \code{contrasts} matrix for comparison}
\usage{
construct_contrasts(design, ...)

\method{construct_contrasts}{default}(design, ...)
}
\arguments{
\item{design}{A design matrix. See \code{\link{construct_design}}.}

\item{...}{A set of contrasts that should be constructed, preferably in 
the format \code{name1 = expression1, name2 = expression2, ...} as it makes 
it clearer to understand the comparisons when looking at the code later 
on. See \code{Examples} section for more. 

See \code{limma::makeContrasts} for other possible ways to specify 
contrasts.}
}
\value{
An object of class matrix with columns corresponding to 
contrasts
}
\description{
\code{construct_contrasts} is a wrapper function to 
\code{limma::makeContrasts}. It allows to specify which comparisons from 
the linear model fit are to be extracted.
}
\examples{
path = system.file("tests", package="ganalyse")

# ----- fpkm ----- # 
fpkm_path = file.path(path, "fpkm", "annotation.txt")
fpkm_obj = rnaseq(fpkm_path, format="fpkm", experiment="sample")
fpkm_counts = gather_counts(fpkm_obj, by="gene-id", log_base=2L)
fpkm_design = construct_design(fpkm_counts, ~ 0 + condition)
fpkm_contrasts = construct_contrasts(
                     design = fpkm_design, 
                     treatA.vs.control = conditiontreatA-conditioncontrol, 
                     treatB.vs.control = conditiontreatB-conditioncontrol
                 )

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

}
\seealso{
\code{\link{rnaseq}}, \code{\link{gather_counts}} 
\code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}} 
\code{\link{construct_design}} \code{\link{write_dge}} 
\code{\link{as.dgelist}} \code{\link{as.eset}} \code{\link{density_plot}}
}

