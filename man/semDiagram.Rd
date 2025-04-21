\name{semDiagram}
\alias{semDiagram}
\title{Visualize SEM Model with Fit Indices}
\description{
  Generates a path diagram of a fitted \code{lavaan} SEM model using \code{DiagrammeR}, with edge widths and colors representing standardized estimates, and fit indices displayed with color-coded thresholds.
}
\usage{
semDiagram(fitted_model, digits = 3, min_width = 1, max_width = 5,
           pos_color = "#1f78b4", neg_color = "#e41a1c")
}
\arguments{
  \item{fitted_model}{A fitted \code{lavaan} model object.}
  \item{digits}{Number of decimal places to display for parameter estimates.}
  \item{min_width}{Minimum edge width in the diagram.}
  \item{max_width}{Maximum edge width in the diagram.}
  \item{pos_color}{Color for positive parameter estimates.}
  \item{neg_color}{Color for negative parameter estimates.}
}
\value{
  An object of class \code{grViz} representing the SEM diagram.
}
\seealso{
  \code{\link[lavaan]{parameterEstimates}}, \code{\link[DiagrammeR]{grViz}}
}
\examples{
library(lavaan)
model <- ' visual  =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed   =~ x7 + x8 + x9 '
fit <- cfa(model, data = HolzingerSwineford1939)
semDiagram(fit)
}
