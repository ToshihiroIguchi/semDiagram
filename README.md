# semDiagram

R functions for visualizing structural equation models (SEM) with fit indices and path diagrams.

## Features

- Extracts model fit indices from lavaan models
- Draws SEM diagrams using DiagrammeR
- Automatically color-codes coefficients and fit indices based on thresholds

## Installation

```r
# Install required packages
install.packages(c("lavaan", "DiagrammeR", "scales"))

# Load functions
source("semDiagram.R")
