\name{semDiagram}
\alias{semDiagram}
\title{Draw SEM Diagrams with Fit Indices Using lavaan and DiagrammeR}
\description{
  This function visualizes structural equation models (SEM) estimated using \code{lavaan}.
  It uses \code{DiagrammeR} to draw the path diagram and \code{scales} to manage edge color transparency.
  Fit indices are also displayed at the top of the graph with color-coded evaluation.
}
\usage{
semDiagram(fitted_model,
           digits = 3,
           standardized = TRUE,
           alpha = 0.05,
           min_width = 1,
           max_width = 5,
           pos_color = "#1f78b4",
           neg_color = "#e41a1c",
           fontname = "Helvetica",
           node_fontsize = 11)
}
\arguments{
  \item{fitted_model}{A model object fitted using \code{lavaan::sem()} or related functions.}
  \item{digits}{Number of digits to round edge labels. Default is 3.}
  \item{standardized}{Logical. If \code{TRUE}, standardized estimates are used.}
  \item{alpha}{Significance level used to determine edge transparency. Default is 0.05.}
  \item{min_width}{Minimum edge width.}
  \item{max_width}{Maximum edge width.}
  \item{pos_color}{Edge color for positive estimates.}
  \item{neg_color}{Edge color for negative estimates.}
  \item{fontname}{Font name used in nodes and edges.}
  \item{node_fontsize}{Font size used in the graph.}
}
\details{
  The function uses Graphviz via \code{DiagrammeR::grViz} to draw a structural equation model (SEM) diagram.
  Latent variables are displayed as ellipses, and observed variables as boxes.
  Edge color and transparency reflect the sign and statistical significance of path coefficients.
  Fit indices such as RMSEA, SRMR, CFI, GFI, etc. are shown at the top of the diagram with threshold-based color coding.
}
\value{
  A \code{grViz} object that renders the diagram in an interactive R environment.
}
\examples{
# Load example data
library(lavaan)
model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
  visual ~ textual + speed
'
fit <- sem(model, data = HolzingerSwineford1939)

# Draw standardized SEM diagram
semDiagram(fit, standardized = TRUE)
}
\seealso{
  \code{\link[lavaan]{sem}}, \code{\link[DiagrammeR]{grViz}}, \code{\link[scales]{alpha}}
}
\author{Toshihiro Iguchi}
