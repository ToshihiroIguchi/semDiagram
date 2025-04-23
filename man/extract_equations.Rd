\name{extract_equations}
\alias{extract_equations}
\title{Generate Simple Equations from a lavaan SEM Model}
\description{
  Extracts regression equations from a fitted \code{lavaan} SEM model and formats them as human-readable strings.
}
\usage{
extract_equations(fit, standardized = TRUE, digits = 3)
}
\arguments{
  \item{fit}{A fitted \code{lavaan} model object, typically obtained using \code{\link[lavaan]{sem}} or \code{\link[lavaan]{cfa}}.}
  \item{standardized}{Logical; if \code{TRUE}, standardized estimates are used. Defaults to \code{TRUE}.}
  \item{digits}{Integer; number of decimal places to round the coefficients. Defaults to \code{3}.}
}
\value{
  A named character vector where each element represents a regression equation in the format \code{"DV = coef1 * IV1 + coef2 * IV2 + ..."}. Each name corresponds to a dependent variable in the model.
}
\examples{
library(lavaan)

# Define a simple SEM model
model <- '
  # Measurement model
  Latent1 =~ x1 + x2 + x3
  Latent2 =~ y1 + y2 + y3
  Latent3 =~ z1 + z2 + z3

  # Structural model
  Latent2 ~ 0.6 * Latent1
  Latent3 ~ -0.4 * Latent1 + 0.7 * Latent2
'

# Simulate data
set.seed(123)
n <- 500
Latent1 <- rnorm(n)
Latent2 <- 0.6 * Latent1 + rnorm(n, sd = 0.5)
Latent3 <- -0.4 * Latent1 + 0.7 * Latent2 + rnorm(n, sd = 0.5)
x1 <- 0.8 * Latent1 + rnorm(n, sd = 0.3)
x2 <- 0.7 * Latent1 + rnorm(n, sd = 0.3)
x3 <- 0.6 * Latent1 + rnorm(n, sd = 0.3)
y1 <- 0.9 * Latent2 + rnorm(n, sd = 0.3)
y2 <- 0.8 * Latent2 + rnorm(n, sd = 0.3)
y3 <- 0.7 * Latent2 + rnorm(n, sd = 0.3)
z1 <- 0.85 * Latent3 + rnorm(n, sd = 0.3)
z2 <- 0.75 * Latent3 + rnorm(n, sd = 0.3)
z3 <- 0.65 * Latent3 + rnorm(n, sd = 0.3)
data <- data.frame(x1, x2, x3, y1, y2, y3, z1, z2, z3)

# Fit the model
fit <- sem(model, data = data)

# Extract equations
equations <- extract_equations(fit)
print(equations)
}
\seealso{
  \code{\link[lavaan]{parameterEstimates}}, \code{\link[lavaan]{standardizedSolution}}, \code{\link[lavaan]{sem}}, \code{\link[lavaan]{cfa}}
}
\author{
  Toshihiro Iguchi
}
\keyword{sem}
\keyword{lavaan}
\keyword{equation}
\keyword{regression}


