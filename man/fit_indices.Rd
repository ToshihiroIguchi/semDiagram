\name{fit_indices}
\alias{fit_indices}
\title{Extract Fit Indices from a lavaan Model}
\description{
  Extracts key fit indices from a fitted \code{lavaan} model object and returns them as a named list.
}
\usage{
fit_indices(model)
}
\arguments{
  \item{model}{A fitted \code{lavaan} model object.}
}
\value{
  A named list containing the following fit indices:
  \itemize{
    \item \code{pvalue} — Model chi-square test p-value
    \item \code{SRMR} — Standardized Root Mean Square Residual
    \item \code{RMSEA} — Root Mean Square Error of Approximation
    \item \code{AIC} — Akaike Information Criterion
    \item \code{GFI} — Goodness of Fit Index
    \item \code{AGFI} — Adjusted Goodness of Fit Index
    \item \code{NFI} — Normed Fit Index
    \item \code{CFI} — Comparative Fit Index
  }
}
\seealso{
  \code{\link[lavaan]{fitMeasures}}
}
\examples{
library(lavaan)
model <- ' visual  =~ x1 + x2 + x3
            textual =~ x4 + x5 + x6
            speed   =~ x7 + x8 + x9 '
fit <- cfa(model, data = HolzingerSwineford1939)
fit_indices(fit)
}
\author{Toshihiro Iguchi}
