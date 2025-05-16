#' Visualize SEM Results as a Path Diagram (no external alpha dependency; robust against NA)
#'
#' @description
#' Draw a publication-ready path diagram from a fitted \code{lavaan} SEM object
#' without relying on additional packages and with sensible defaults that remain
#' robust when statistics are \code{NA}. Supports both standardized and
#' unstandardized coefficients, automatic edge transparency by p-value, fit
#' indices in the title, and several Graphviz layout engines.
#'
#' @param fitted_model    A \code{lavaan} model object.
#' @param digits          Number of digits for printed estimates.
#' @param standardized    If \code{TRUE}, use fully standardized coefficients (\code{std.all})
#'                        for labels and widths; otherwise raw estimates.
#' @param alpha           p-value threshold for edge opacity.
#' @param low_alpha       Transparency for non-significant edges and \code{NA} p-values.
#' @param min_width       Minimum edge width.
#' @param max_width       Maximum edge width.
#' @param pos_color       Color for positive edges.
#' @param neg_color       Color for negative edges.
#' @param fontname        Font family for all labels and node text.
#' @param node_fontsize   Base font size for nodes.
#' @param edge_fontsize   Base font size for edge labels.
#' @param show_residuals  Whether to draw residual-variance edges.
#' @param show_intercepts Whether to draw intercept edges.
#' @param show_fit        Whether to include fit indices label at the top.
#' @param layout          Graph orientation: \code{"LR"} (left→right) or \code{"TB"} (top→bottom).
#' @param ratio           Graphviz \code{ratio} attribute (\code{"fill"} or numeric).
#' @param curvature       Curvature for bidirectional edges.
#' @param engine          Graphviz layout engine (\code{"dot"}, \code{"neato"}, \code{"fdp"},
#'                        \code{"circo"}, \code{"twopi"}).
#' @param twopi_compact   If \code{TRUE} (default) compacts \code{twopi} by shrinking
#'                        \code{ranksep}, \code{nodesep}, and \code{sep}.
#'
#' @return A DiagrammeR \code{htmlwidget} displaying the path diagram.
#' @export
semDiagram <- function(
    fitted_model,
    digits          = 3,
    standardized    = TRUE,
    alpha           = 0.05,
    low_alpha       = 0.2,
    min_width       = 1,
    max_width       = 5,
    pos_color       = "blue",
    neg_color       = "red",
    fontname        = "Helvetica",
    node_fontsize   = 11,
    edge_fontsize   = 9,
    show_residuals  = FALSE,
    show_intercepts = FALSE,
    show_fit        = TRUE,
    layout          = "LR",
    ratio           = "fill",
    curvature       = 0.3,
    engine          = "dot",
    twopi_compact   = TRUE) {

  ## ---------- 1. Validate & helpers ----------
  engine <- match.arg(engine, c("dot","neato","fdp","circo","twopi"))

  alpha_color <- function(col, alpha_val) {
    if (is.na(alpha_val) || alpha_val < 0 || alpha_val > 1) alpha_val <- 1
    rgb_mat <- grDevices::col2rgb(col) / 255
    grDevices::rgb(rgb_mat[1,], rgb_mat[2,], rgb_mat[3,], alpha = alpha_val)
  }

  if (!inherits(fitted_model, "lavaan"))
    stop("`fitted_model` must be a lavaan object.")

  invisible(lapply(c("DiagrammeR","lavaan"), function(p)
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", p))))

  ## ---------- 2. lavaan statistics ----------
  params    <- lavaan::parameterEstimates(fitted_model, standardized = TRUE)
  scale_col <- if (standardized) "std.all" else "est"

  fit   <- lavaan::fitMeasures(
    fitted_model,
    c("pvalue","srmr","rmsea","gfi","agfi","nfi","cfi","aic","bic"))
  n_obs <- lavaan::lavInspect(fitted_model, "nobs")

  ## ---------- 3. Edge-width scaling ----------
  edge_rows <- params[params$op %in% c("=~","~","~~") & params$lhs!=params$rhs, ]
  vals <- abs(edge_rows[[scale_col]]); vals <- vals[!is.na(vals)]
  # When standardized, fix max_abs to 1 so scaling respects absolute bounds
  if (standardized) {
    max_abs <- 1
  } else {
    max_abs <- if (length(vals)==0 || !is.finite(max(vals))) 1 else max(vals)
  }

  ## ---------- 4. Color helper ----------
  colorize_fit <- function(v, thr, inv = FALSE) {
    if (is.na(v)) "gray50"
    else if (inv && v > thr) "red"
    else if (!inv && v <= thr) "red"
    else "gray20"
  }

  ## ---------- 5. Top-label construction ----------
  # 5-1) Fit indices (inner only)
  if (show_fit) {
    fit_inner <- paste0(
      sprintf("N = %d | ", n_obs),
      sprintf("<font color='%s'>p = %.3f</font> | ",
              colorize_fit(fit["pvalue"], 0.05), fit["pvalue"]),
      sprintf("<font color='%s'>SRMR = %.3f</font> | ",
              colorize_fit(fit["srmr"], 0.08, TRUE), fit["srmr"]),
      sprintf("<font color='%s'>RMSEA = %.3f</font> | ",
              colorize_fit(fit["rmsea"], 0.08, TRUE), fit["rmsea"]),
      sprintf("AIC = %.1f | BIC = %.1f<BR/>", fit["aic"], fit["bic"]),
      sprintf("<font color='%s'>GFI = %.3f</font> | ",
              colorize_fit(fit["gfi"], 0.90), fit["gfi"]),
      sprintf("<font color='%s'>AGFI = %.3f</font> | ",
              colorize_fit(fit["agfi"], 0.90), fit["agfi"]),
      sprintf("<font color='%s'>NFI = %.3f</font> | ",
              colorize_fit(fit["nfi"], 0.90), fit["nfi"]),
      sprintf("<font color='%s'>CFI = %.3f</font>",
              colorize_fit(fit["cfi"], 0.90), fit["cfi"])
    )
  }

  # 5-2) Annotation row
  coeff_text <- if (standardized)
    "Coefficients: <b>Standardized</b>"
  else
    "Coefficients: <b>Unstandardized</b>"
  inter_text <- if (show_intercepts)
    "Intercepts: <b>Shown</b>"
  else
    " | Intercepts: <b>Hidden</b>"
  annot_inner <- sprintf("<font point-size='12'>%s | %s</font>",
                         coeff_text, inter_text)

  # 5-3) Combine and wrap with single < … >
  top_label <- if (show_fit)
    sprintf("<%s<BR/>%s>", annot_inner, fit_inner)
  else
    sprintf("<%s>", annot_inner)

  ## ---------- 6. Nodes ----------
  latents   <- unique(params$lhs[params$op == "=~"])
  observeds <- setdiff(unique(c(params$lhs, params$rhs)),
                       c(latents, "1", ""))
  nodes <- list()
  for (lv in latents)
    nodes[[lv]] <- list(shape = "ellipse", label = lv, style = "filled",
                        fillcolor = "#F0F0F0",
                        fontname = fontname, fontsize = node_fontsize)
  for (ov in observeds)
    nodes[[ov]] <- list(shape = "box", label = ov,
                        fontname = fontname, fontsize = node_fontsize)

  ## ---------- 7. Edges ----------
  edges <- list()
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    if (p$op %in% c("=~","~","~~") && p$lhs != p$rhs) {
      value <- if (standardized) p$std.all else p$est
      if (is.na(value)) value <- p$est
      pen <- (abs(value) / max_abs) * (max_width - min_width) + min_width
      if (!is.finite(pen)) pen <- min_width

      alpha_edge <- if (is.na(p$pvalue)) low_alpha
      else if (p$pvalue < alpha) 1
      else low_alpha
      col <- alpha_color(if (value >= 0) pos_color else neg_color,
                         alpha_edge)

      e_base <- switch(p$op,
                       "=~" = list(from = p$lhs, to = p$rhs, arrowhead = "vee"),
                       "~"  = list(from = p$rhs, to = p$lhs, arrowhead = "vee"),
                       "~~" = list(from = p$lhs, to = p$rhs, arrowhead = "vee",
                                   arrowtail = "vee",
                                   dir = "both", style = "dashed")
      )
      e_base$label    <- sprintf("%.*f", digits, value)
      e_base$penwidth <- pen
      e_base$color    <- col
      e_base$fontsize <- edge_fontsize
      e_base$fontname <- fontname
      if (p$op == "~~") {
        e_base$constraint <- FALSE
        e_base$dir        <- "both"
      }
      edges[[length(edges) + 1]] <- e_base
    }
  }

  ## ---------- 8. DOT-code construction ----------
  node_defs <- paste(
    vapply(names(nodes), function(n) {
      attrs <- paste(names(nodes[[n]]),
                     vapply(nodes[[n]], function(x)
                       if (is.numeric(x)) as.character(x)
                       else sprintf("\"%s\"", x),
                       character(1)),
                     sep = "=", collapse = ", ")
      sprintf("  \"%s\" [%s];", n, attrs)
    }, character(1)), collapse = "\n")

  edge_defs <- paste(
    vapply(edges, function(e) {
      attrs <- paste(names(e)[-1:-2],
                     vapply(e[-1:-2], function(x)
                       if (is.numeric(x)) as.character(x)
                       else sprintf("\"%s\"", x),
                       character(1)),
                     sep = "=", collapse = ", ")
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }, character(1)), collapse = "\n")

  radial_opts <- if (engine == "circo") {
    c("splines=true", "nodesep=0.4", "sep=\"+4\"", "mindist=1")
  } else if (engine == "twopi" && twopi_compact) {
    c("splines=true", "nodesep=0.2", "sep=\"+2\"",
      "ranksep=0.5", "normalize=true")
  } else character(0)
  radial_opts <- paste(radial_opts, collapse = ", ")

  graph_code <- sprintf(
    "digraph {\n  rankdir=%s;\n  graph [layout=%s%s%s, overlap=false,\n         labelloc=\"t\", labeljust=\"c\", label=%s, ratio=%s];\n  node  [fontname=\"%s\", margin=0.05];\n  edge  [fontname=\"%s\", fontcolor=\"#333333\"];\n\n%s\n\n%s\n}",
    layout, engine,
    if (nchar(radial_opts)) ", " else "", radial_opts,
    top_label, ratio,
    fontname, fontname,
    node_defs, edge_defs)

  ## ---------- 9. Render ----------
  DiagrammeR::grViz(graph_code, engine = engine)
}
