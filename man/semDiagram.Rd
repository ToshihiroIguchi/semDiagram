\name{semDiagram}
\alias{semDiagram}
\title{Visualize Structural Equation Models as Diagrams}
\description{
  This function visualizes a structural equation model (SEM) fitted by the \code{lavaan} package using the \code{DiagrammeR} package. It includes latent and observed variables, standardized coefficients, significance-based edge transparency, and optional residuals, intercepts, and fit indices.
}
\usage{
semDiagram(fitted_model,
           digits = 3,
           standardized = TRUE,
           alpha = 0.05,
           min_width = 1, max_width = 5,
           pos_color = "blue", neg_color = "red",
           fontname = "Helvetica",
           node_fontsize = 11,
           edge_fontsize = 9,
           show_residuals = FALSE,
           show_intercepts = TRUE,
           show_fit = TRUE,
           layout = "LR",
           curvature = 0.3)
}
\arguments{
  \item{fitted_model}{A fitted \code{lavaan} model object.}
  \item{digits}{Number of digits to round edge labels (default is 3).}
  \item{standardized}{Logical. If \code{TRUE}, standardized estimates are used.}
  \item{alpha}{Significance threshold for edge transparency.}
  \item{min_width}{Minimum edge width for smallest effect.}
  \item{max_width}{Maximum edge width for largest effect.}
  \item{pos_color}{Color for positive coefficients.}
  \item{neg_color}{Color for negative coefficients.}
  \item{fontname}{Font name used in the diagram.}
  \item{node_fontsize}{Font size for node labels.}
  \item{edge_fontsize}{Font size for edge labels.}
  \item{show_residuals}{Logical. If \code{TRUE}, residual variances are shown.}
  \item{show_intercepts}{Logical. If \code{TRUE}, intercepts are shown.}
  \item{show_fit}{Logical. If \code{TRUE}, model fit indices are displayed.}
  \item{layout}{Graph layout direction. One of \code{"LR"}, \code{"TB"}, etc.}
  \item{curvature}{Not used currently. Reserved for future curvature options.}
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
