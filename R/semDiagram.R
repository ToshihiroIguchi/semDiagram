
#' Visualize SEM Results as a Path Diagram (optimized: no library()/install, fixed vapply)
#'
#' @param fitted_model An object of class lavaan
#' @param digits Integer: number of digits for estimates
#' @param standardized Logical: use standardized estimates?
#' @param alpha Numeric: significance threshold for edge transparency
#' @param min_width Numeric: minimum edge width
#' @param max_width Numeric: maximum edge width
#' @param pos_color,neg_color Character: colors for positive/negative estimates
#' @param fontname Character: font family for labels
#' @param node_fontsize,edge_fontsize Numeric: font sizes
#' @param show_residuals,show_intercepts,show_fit Logical: whether to show residuals/intercepts/fit indices
#' @param layout Character: "LR" or "TB"
#' @param ratio Character or numeric: GraphViz 'ratio' attribute (e.g., "fill", numeric aspect ratio)
#' @param curvature Numeric: curvature for bidirectional edges
#' @export
semDiagram <- function(fitted_model,
                       digits = 3,
                       standardized = TRUE,
                       alpha = 0.05,
                       min_width = 1,
                       max_width = 5,
                       pos_color = "blue",
                       neg_color = "red",
                       fontname = "Helvetica",
                       node_fontsize = 11,
                       edge_fontsize = 9,
                       show_residuals = FALSE,
                       show_intercepts = FALSE,
                       show_fit = TRUE,
                       layout = "LR",
                       ratio = "fill",
                       curvature = 0.3) {

  # 1. Check that fitted_model is lavaan object
  if (!inherits(fitted_model, "lavaan")) {
    stop("`fitted_model` must be a lavaan model object.")
  }

  # 2. Dependency check: no install or library calls at runtime
  pkgs <- c("scales", "DiagrammeR", "lavaan")
  invisible(lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Please install it via install.packages('%s').", pkg, pkg))
    }
  }))

  # 3. Retrieve parameter estimates
  params <- lavaan::parameterEstimates(fitted_model, standardized = standardized)

  # 4. Retrieve fit measures
  fit <- lavaan::fitMeasures(fitted_model, c("pvalue", "srmr", "rmsea", "gfi", "agfi", "nfi", "cfi", "aic", "bic"))

  # 5. Colorizing function for fit indices
  colorize <- function(value, threshold, inverse = FALSE) {
    if (is.na(value)) return("gray50")
    if (inverse) {
      if (value > threshold) "red" else "gray20"
    } else {
      if (value <= threshold) "red" else "gray20"
    }
  }

  # 6. Build fit label
  if (show_fit) {
    p_col    <- colorize(fit["pvalue"], 0.05)
    s_col    <- colorize(fit["srmr"],   0.08, inverse = TRUE)
    r_col    <- colorize(fit["rmsea"],  0.08, inverse = TRUE)
    gfi_col  <- colorize(fit["gfi"],    0.90)
    agfi_col <- colorize(fit["agfi"],   0.90)
    nfi_col  <- colorize(fit["nfi"],    0.90)
    cfi_col  <- colorize(fit["cfi"],    0.90)

    fit_label <- paste0(
      "<<font point-size='12'>",
      sprintf("<font color='%s'>p = %.3f</font>",   p_col,    fit["pvalue"]), " | ",
      sprintf("<font color='%s'>SRMR = %.3f</font>", s_col,    fit["srmr"]),   " | ",
      sprintf("<font color='%s'>RMSEA = %.3f</font>",r_col,    fit["rmsea"]),  " | ",
      sprintf("AIC = %.1f", fit["aic"]),             " | ",
      sprintf("BIC = %.1f", fit["bic"]),             "<br/>",
      sprintf("<font color='%s'>GFI = %.3f</font>",  gfi_col,  fit["gfi"]),    " | ",
      sprintf("<font color='%s'>AGFI = %.3f</font>", agfi_col, fit["agfi"]),   " | ",
      sprintf("<font color='%s'>NFI = %.3f</font>",  nfi_col,  fit["nfi"]),    " | ",
      sprintf("<font color='%s'>CFI = %.3f</font>",  cfi_col,  fit["cfi"]),
      "</font>>"
    )
  } else {
    fit_label <- ""
  }

  # 7. Define nodes
  latent_vars   <- unique(as.character(params$lhs[params$op == "=~"]))
  observed_vars <- unique(as.character(c(params$lhs, params$rhs)))
  observed_vars <- setdiff(observed_vars, c(latent_vars, "1", ""))

  nodes <- list()
  for (lv in latent_vars) {
    nodes[[lv]] <- list(
      shape     = "ellipse",
      label     = lv,
      style     = "filled",
      fillcolor = "#F0F0F0",
      fontname  = fontname,
      fontsize  = node_fontsize
    )
  }
  for (ov in observed_vars) {
    nodes[[ov]] <- list(
      shape    = "box",
      label    = ov,
      fontname = fontname,
      fontsize = node_fontsize
    )
  }

  # 8. Define edges
  edges <- list()
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    if (p$op %in% c("=~", "~", "~~") && p$lhs != p$rhs) {
      pen        <- abs(p$est) * (max_width - min_width) + min_width
      alpha_edge <- ifelse(p$pvalue < alpha, 1, 0.3)
      col        <- scales::alpha(ifelse(p$est >= 0, pos_color, neg_color), alpha_edge)

      edge_def <- switch(
        p$op,
        "=~" = list(from = p$lhs, to = p$rhs, arrowhead = "vee"),
        "~"  = list(from = p$rhs, to = p$lhs, arrowhead = "vee"),
        "~~" = list(from = p$lhs, to = p$rhs,
                    arrowhead = "vee", arrowtail = "vee",
                    dir = "both", style = "dashed")
      )
      # Use standardized estimate if requested
      value <- if (standardized) p$std.all else p$est
      edge_def$label    <- sprintf("%.*f", digits, value)
      edge_def$penwidth <- pen
      edge_def$color    <- col
      edge_def$fontsize <- edge_fontsize
      edge_def$fontname <- fontname
      if (p$op == "~~") {
        edge_def$constraint <- FALSE
        edge_def$dir        <- "both"
      }
      edges[[length(edges) + 1]] <- edge_def
    }
  }

  # 9. Add residual nodes
  if (show_residuals) {
    resids <- params[params$op == "~~" & params$lhs == params$rhs, ]
    for (i in seq_len(nrow(resids))) {
      r       <- resids[i, ]
      node_id <- paste0("res_", r$lhs)
      nodes[[node_id]] <- list(
        shape    = "circle",
        label    = paste0("ε_", r$lhs),
        fontname = fontname,
        fontsize = node_fontsize - 1,
        width    = 0.5,
        height   = 0.5
      )
      edges[[length(edges) + 1]] <- list(
        from     = node_id,
        to       = r$lhs,
        label    = sprintf("%.*f", digits, r$est),
        arrowhead= "vee",
        style    = "dotted",
        penwidth = 1,
        color    = "gray50",
        fontsize = edge_fontsize
      )
    }
  }

  # 10. Add intercept nodes
  if (show_intercepts) {
    ints <- params[params$op == "~1", ]
    if (nrow(ints) > 0) {
      for (i in seq_len(nrow(ints))) {
        p       <- ints[i, ]
        node_id <- paste0("int_", p$lhs)
        nodes[[node_id]] <- list(
          shape     = "triangle",
          label     = "1",
          style     = "filled",
          fillcolor = "#DDDDDD",
          fontname  = fontname,
          fontsize  = node_fontsize * 0.5,
          width     = 0,
          height    = 0
        )
        edges[[length(edges) + 1]] <- list(
          from     = node_id,
          to       = p$lhs,
          label    = sprintf("%.*f", digits, p$est),
          arrowhead= "vee",
          style    = "solid",
          penwidth = 1,
          color    = "gray50",
          fontsize = edge_fontsize,
          fontname = fontname
        )
      }
    }
  }

  # 11. Generate GraphViz code
  graph_code <- paste0(
    "digraph {\n",
    "  rankdir = ", layout, ";\n",
    "  graph [overlap=false, labelloc='t', labeljust='c', label=", fit_label, ", ratio=", ratio, "];\n",
    "  node [fontname='", fontname, "', margin=0.05];\n",
    "  edge [fontname='", fontname, "', fontcolor='#333333'];\n\n",
    paste(sapply(names(nodes), function(n) {
      attrs <- paste(
        names(nodes[[n]]),
        vapply(nodes[[n]], function(x) if (is.numeric(x)) as.character(x) else paste0("\"", x, "\""), character(1)),
        sep = "=",
        collapse = ", "
      )
      sprintf("  \"%s\" [%s];", n, attrs)
    }), collapse = "\n"), "\n\n",
    paste(sapply(edges, function(e) {
      attrs <- paste(
        names(e)[-1:-2],
        vapply(e[-1:-2], function(x) if (is.numeric(x)) as.character(x) else paste0("\"", x, "\""), character(1)),
        sep = "=",
        collapse = ", "
      )
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }), collapse = "\n"), "\n",
    "}\n"
  )

  # 12. Render graph
  DiagrammeR::grViz(graph_code)
}
