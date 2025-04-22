# Define the function to extract equations
extract_equations <- function(fit, standardized = TRUE, digits = 3) {
  # Extract parameter estimates
  params <- parameterEstimates(fit, standardized = standardized)

  # Filter for regression paths
  regressions <- subset(params, op == "~")

  # Initialize a list to store equations
  equations <- list()

  # Loop through each unique dependent variable
  for (dv in unique(regressions$lhs)) {
    # Subset for current dependent variable
    subset_reg <- subset(regressions, lhs == dv)

    # Construct the right-hand side of the equation
    rhs_terms <- paste0(
      format(round(subset_reg$est, digits), nsmall = digits),
      " * ",
      subset_reg$rhs
    )

    # Combine into a single equation
    equation <- paste(dv, "=", paste(rhs_terms, collapse = " + "))

    # Append to the list
    equations[[dv]] <- equation
  }

  # Return the list of equations
  return(equations)
}
