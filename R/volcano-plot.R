#' @title Volcano plot for DGE analysis
#' 
#' @description \code{volcano_plot} takes an object of class \code{dge} and 
#' returns a volcano plot. If \code{filename} is provided, the plot is also 
#' saved to the file. 
#' 
#' @details Both \code{\link{limma_dge}} or \code{\link{edger_dge}} 
#' methods return \code{dge} object which can be directly passed to 
#' \code{volcano_plot}. It inherits from \code{data.table}. 
#'
#' @param x An object of class \code{dge} or a list of objects of class 
#' \code{dge}.
#' @param interactive Default is \code{FALSE}. If \code{TRUE}, uses 
#' \code{ggplotly::ggplotly()} to return an interactive plot.
#' @param title Plot title.
#' @param labels Integer vector of length-1 indicating the number of 
#' top genes that should be labeled. Default is 6.
#' @param filename Default is to plot to screen. If a file name 
#' is provided, the plot is saved to file and the plot object will be 
#' returned. The type of graphic device is auto-detected from file extension 
#' (using \code{ggsave}).
#' @param height Height of the plot, default is 8 inches.
#' @param width Width of the plot, default is 12 inches.
#' @param ... Optional arguments to pass to \code{ggsave} call.
#' 
#' @return A \code{ggplot2} or \code{plotly} object containing the 
#' volcano plot of input \code{dge} object.
#' 
#' @export
#' @aliases volcano_plot volcano-plot
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}}
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
#' \code{\link{density_plot}}
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
#' volcano_plot(ans, interactive=FALSE) # ggplot2 plot
#' volcano_plot(ans, interactive=FALSE, title="A.vs.control")
#' volcano_plot(ans, interactive=TRUE, title="A.vs.control") # ggplotly plot
#' \dontrun{
#' volcano_plot(ans, filename="tmp.png") # write to file and return plot obj
#' }
volcano_plot <- function(x, interactive=FALSE, title, labels=6L, 
                            filename, height=8L, width=12L, ...) {
  UseMethod("volcano_plot")
}

#' @rdname volcano_plot
#' @export
volcano_plot.default <- function(x, interactive=FALSE, title, labels=6L, 
                            filename, height=8L, width=12L, ...) {
    stop("No default method available for object of class 'volcano_plot'")
}

#' @rdname volcano_plot
#' @export
volcano_plot.dge <- function(x, interactive=FALSE, title, labels=6L, 
                        filename, height=8L, width=12L, ...) {
    interactive = suppressWarnings(as.logical(interactive[1L]))
    if (!interactive %in% c(FALSE, TRUE))
        stop("Argument interactive must be logical TRUE/FALSE.")
    if (!all(c("padj", "log2fc") %in% names(x)))
        stop("Columns padj and log2fc must be present in input object.")
    padj=log2fc=log10padj=id=NULL
    x = copy(x)
    x[, "log10padj" := -log10(padj)]
    idx = head(which(x$padj < 0.05), max(0L, labels, na.rm=TRUE))
    # disabling coloring for now
    # [, "colour" := paste0(v1*1L, v2*1L), 
    #   by=.(v1=(2^log2FC >= 2 | 2^log2FC <= 0.5), v2=padj <= 0.05)][]
    if (!requireNamespace("ggplot2"))
        stop("Package ggplot2 is not available.")
    if (interactive && !requireNamespace("plotly"))
        stop("Package plotly is not available.")
    # if (interactive && !requireNamespace("htmltools"))
    #     stop("Package htmltools is not available.")
    if (interactive || !labels) {
        repel = FALSE
    } else {
        # if labels=0L, then there's no need for 'ggrepel'
        if (!requireNamespace("ggrepel")) {
            warning("Package ggrepel is not available. ", 
                "Top genes will not be marked.")
        }
        repel = TRUE
    }
    ggp = ggplot2::ggplot(x, ggplot2::aes(x=log2fc, y=log10padj, 
                        text=paste("id:", id))) + 
            ggplot2::geom_point(color="#B8A6B1", alpha=0.5) + 
            (if (sum(x[["padj"]] < 0.05))
                ggplot2::geom_point(fill="#3F3D45", color="#3F3D45", 
                    shape=21L, alpha=0.8, data=x[padj<0.05])) + 
            (if (length(idx)) ggplot2::geom_point(shape=21L, colour="#5E586E", 
                fill="#3F3D45", alpha=0.8, size=4L, data=x[(idx)])) + 
            ggplot2::xlab("log2 fold change") + 
            ggplot2::ylab("log10 p-value") + 
            ggplot2::theme_minimal(base_size = 15L) +
            ggplot2::theme(legend.text=ggplot2::element_text(size = 12), 
                    legend.background=ggplot2::element_blank(),
                    legend.key=ggplot2::element_blank(),
                    panel.grid.major=ggplot2::element_line(colour="#333333", 
                                                            size=0.12), 
                    panel.grid.minor=ggplot2::element_line(colour="#999999", 
                                                            size=0.06)
                    # color="#C986AF", 
                    # color="#A75486", 
                    # axis.title=element_text(face="bold"), 
                    # axis.title=element_text(size=16,face="bold"), 
                    # legend.title = element_text(size=14), 
                    # axis.text=element_text(size=14)
                    # strip.text=element_text(size = 14), 
            # scale_colour_manual(name = "FC >= 2 | <= 0.5,\npval <= 0.05", 
            #           values = c(brewer.pal(8, "Set1"))) + 
            # geom_hline(aes(yintercept=-log10(0.05)), linetype=3, 
                        # colour="#641144", alpha=0.5) + 
            # facet_grid( method ~ comparison, scales="free") + 
                )
    if (!missing(title)) ggp = ggp + ggplot2::ggtitle(title)
    if (repel & length(idx)) {
        ggp = ggp + ggrepel::geom_text_repel(data=x[(idx)], 
                ggplot2::aes(label=id), size=3L, segment.color = "#4F444B")
    }
    if (!missing(filename)) 
        ggplot2::ggsave(filename, plot=ggp, height=height, width=width, ...)
    if (interactive) 
        ggp = plotly::ggplotly(ggp, tooltop = c("id", "log2fc", "padj"))
    ggp
}

# TODO: add method for list od 'dge' objects
# #' @rdname volcano_plot
# #' @export
# volcano_plot.list <- function(x, interactive=FALSE, filename, 
#                         height=8L, width=12L, ...) {
#     if (!all(vapply(x, is.dge, TRUE)))
#         stop("All elements of input list must be objects of class dge.")
#     x = rbindlist(x, idcol=".title")
#     volcano_plot(x, ...)
# }

# ----- Internals ----- #

is.dge <- function(x) inherits(x, "dge")
