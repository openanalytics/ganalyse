#' @title Density plots of count data facetted by groups
#' 
#' @description \code{density_plot} function takes an object of class 
#' \code{fpkm_counts} or \code{raw_counts} and returns density plots of all 
#' samples facetted by the groups specified.
#' 
#' @param x An object of class \code{dge} or a list of objects of class 
#' \code{dge}.
#' @param interactive Default is \code{FALSE}. If \code{TRUE}, uses 
#' \code{ggplotly::ggplotly()} to return an interactive plot.
#' @param title Plot title.
#' @param groups Columns indicating those samples that needs to be 
#' selected together while facetting.
#' @param facet_cols Number of columns to provide to \code{facet_wrap}.
#' @param log Should values be log-transformed? Default is \code{FALSE}.
#' @param filename Default is to plot to screen. If a file name 
#' is provided, the plot is saved to file and the plot object will be 
#' returned. The type of graphic device is auto-detected from file extension 
#' (using \code{ggsave}).
#' @param height Height of the plot, default is 8 inches.
#' @param width Width of the plot, default is 12 inches.
#' @param ... Optional arguments to pass to \code{ggsave} call.
#' 
#' @return A \code{ggplot2} or \code{plotly} object containing the 
#' density plot of input \code{fpkm_counts} or \code{raw_counts} object.
#' 
#' @export
#' @aliases density_plot density-plot
#' @seealso \code{\link{rnaseq}}, \code{\link{gather_counts}} 
#' \code{\link{show_counts}} \code{\link{limma_dge}} \code{\link{edger_dge}}
#' \code{\link{construct_design}} \code{\link{construct_contrasts}} 
#' \code{\link{write_dge}} \code{\link{as.dgelist}} \code{\link{as.eset}}
#' @examples
#' path = system.file("tests", package="ganalyse")
#' 
#' # ----- raw ----- # 
#' raw_path = file.path(path, "raw", "annotation.txt")
#' raw_obj = rnaseq(raw_path, format="raw", experiment="sample")
#' raw_counts = gather_counts(raw_obj, by="gene-id", threshold=1L)
#' 
#' # DE genes between treatA and control
#' density_plot(raw_counts, interactive=FALSE, groups="condition") # ggplot2
#' density_plot(raw_counts, interactive=FALSE, groups="condition", log=TRUE)
#' density_plot(raw_counts, interactive=FALSE, 
#'         groups="condition", title="Raw count density")
#' density_plot(raw_counts, interactive=TRUE, 
#'         groups="condition", log=TRUE) # ggplotly
#' \dontrun{
#' # write to file and return plot obj
#' density_plot(raw_counts, filename="tmp.png")
#' }
density_plot <- function(x, interactive=FALSE, title, groups, facet_cols=1L, 
                        log=FALSE, filename, height=8L, width=12L, ...) {
  UseMethod("density_plot")
}

#' @rdname density_plot
#' @export
density_plot.default <- function(x, interactive=FALSE, title, groups, 
                            facet_cols=1L, log=FALSE, filename, height=8L, 
                            width=12L, ...) {
    stop("No default method available for object of class 'density_plot'")
}

#' @rdname density_plot
#' @export
density_plot.fpkm_counts <- function(x, interactive=FALSE, title, groups, 
    facet_cols=1L, log=FALSE, filename, height=8L, width=12L, ...) {

    interactive = suppressWarnings(as.logical(interactive[1L]))
    if (!interactive %in% c(FALSE, TRUE))
        stop("Argument interactive must be logical TRUE/FALSE.")
    if (!is.character(groups) || !all(groups %in% names(x))) 
        stop("groups must be a character vector of column names ", 
            "in input object 'x', corresponding to the groups each ", 
            "sample belongs to")
    log = suppressWarnings(as.logical(log[1L]))
    if (is.na(log) || !length(log))
        stop("log must be logical TRUE/FALSE")

    facet_cols = suppressWarnings(as.integer(facet_cols[1L]))
    if (is.na(facet_cols) || !length(facet_cols) || facet_cols < 1L)
        stop("facet_cols must be a positive integer value.")

    if (!requireNamespace("ggplot2"))
        stop("Package ggplot2 is not available.")
    if (interactive && !requireNamespace("plotly"))
        stop("Package plotly is not available.")
    value=.GROUPS=NULL
    mx = melt(show_counts(x), id="id", variable.name="sample")
    mx[x, c(groups) := mget(paste0("i.", groups)), on="sample"]
    mx[, ".GROUPS" := do.call(paste, c(.SD, sep="_")), .SDcols=groups]
    facet_formula = as.formula(paste("~", paste(groups, collapse="+")))
    ggp = ggplot2::ggplot(mx, ggplot2::aes(value)) + 
            ggplot2::geom_density(ggplot2::aes(group=sample, colour=.GROUPS)) + 
            ggplot2::facet_wrap(facet_formula, ncol=facet_cols) + 
            ggplot2::xlab("Counts") + 
            ggplot2::ylab("Density") + 
            ggplot2::theme_minimal(base_size = 15L) +
            ggplot2::theme(legend.text=ggplot2::element_text(size = 12), 
                    legend.background=ggplot2::element_blank(),
                    legend.key=ggplot2::element_blank(),
                    panel.grid.major=ggplot2::element_line(colour="#333333", 
                                                            size=0.08), 
                    panel.grid.minor=ggplot2::element_line(colour="#999999", 
                                                            size=0.04)) + 
            ggplot2::guides(colour=FALSE) + 
            (if (log) ggplot2::scale_x_log10()) + ggplot2::xlab("log(Counts)")
    if (!missing(title)) ggp = ggp + ggplot2::ggtitle(title)
    if (!missing(filename)) 
        ggplot2::ggsave(filename, plot=ggp, height=height, width=width, ...)
    if (interactive) 
        ggp = plotly::ggplotly(ggp)
    ggp
}

#' @rdname density_plot
#' @export
density_plot.raw_counts <- density_plot.fpkm_counts

# TODO: add method for list od 'dge' objects
# #' @rdname density_plot
# #' @export
# density_plot.list <- function(x, interactive=FALSE, filename, 
#                         height=8L, width=12L, ...) {
#     if (!all(vapply(x, is.dge, TRUE)))
#         stop("All elements of input list must be objects of class dge.")
#     x = rbindlist(x, idcol=".title")
#     density_plot(x, ...)
# }

# ----- Internals ----- #

is.dge <- function(x) inherits(x, "dge")
