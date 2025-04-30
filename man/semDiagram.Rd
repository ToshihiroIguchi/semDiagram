\name{semDiagram}
\alias{semDiagram}
\title{Visualize SEM Results as a Path Diagram}
\description{
  The \code{semDiagram} function visualizes structural equation modeling (SEM) results as a path diagram using the DiagrammeR package.
}
\usage{
semDiagram(fitted_model, digits = 3, standardized = TRUE, alpha = 0.05,
  min_width = 1, max_width = 5, pos_color = "blue", neg_color = "red",
  fontname = "Helvetica", node_fontsize = 11, edge_fontsize = 9,
  show_residuals = FALSE, show_intercepts = FALSE, show_fit = TRUE,
  layout = "LR", curvature = 0.3)
}
\arguments{
  \item{fitted_model}{An object of class \code{lavaan}, representing the SEM model fit.}
  \item{digits}{An integer specifying the number of digits for parameter estimates in the diagram.}
  \item{standardized}{Logical. Whether to use standardized estimates. Defaults to \code{TRUE}.}
  \item{alpha}{Significance level for coloring edges based on p-values. Defaults to 0.05.}
  \item{min_width}{Minimum width of edges in the diagram. Defaults to 1.}
  \item{max_width}{Maximum width of edges in the diagram. Defaults to 5.}
  \item{pos_color}{Color for positive parameter estimates. Defaults to \code{"blue"}.}
  \item{neg_color}{Color for negative parameter estimates. Defaults to \code{"red"}.}
  \item{fontname}{Font family to use in the diagram. Defaults to \code{"Helvetica"}.}
  \item{node_fontsize}{Font size for node labels. Defaults to 11.}
  \item{edge_fontsize}{Font size for edge labels. Defaults to 9.}
  \item{show_residuals}{Logical. Whether to display residual variances as separate nodes. Defaults to \code{FALSE}.}
  \item{show_intercepts}{Logical. Whether to display intercepts as separate nodes. Defaults to \code{FALSE}.}
  \item{show_fit}{Logical. Whether to display model fit indices. Defaults to \code{TRUE}.}
  \item{layout}{Layout direction of the diagram. Options include \code{"LR"} (left-to-right) and \code{"TB"} (top-to-bottom). Defaults to \code{"LR"}.}
  \item{curvature}{Curvature of bidirectional edges (for covariances). Defaults to 0.3.}
}
\value{
  A graph object rendered by \code{DiagrammeR::grViz}, showing the SEM diagram.
}
\details{
  This function checks for required packages and installs any missing ones. It extracts parameter estimates and fit measures from a \code{lavaan} model and constructs a diagram using GraphViz. Edge width represents effect size, transparency indicates significance, and colors show sign.
}
\seealso{
  \code{\link[lavaan]{lavaan}}, \code{\link[DiagrammeR]{grViz}}
}
\examples{
\dontrun{
library(lavaan)
model <- ' visual  =~ x1 + x2 + x3
           textual =~ x4 + x5 + x6
           speed   =~ x7 + x8 + x9 '
fit <- cfa(model, data = HolzingerSwineford1939)
semDiagram(fit)
}
}
\author{Toshihiro Iguchi}
