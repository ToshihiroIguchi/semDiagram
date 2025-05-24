#' Visualize SEM Results as a Path Diagram (no external alpha dependency; robust against NA)
#'
#' @description
#' Draw a publication-ready path diagram from a fitted \code{lavaan} SEM object.
#' Added functionality: automatically computes and (optionally) shows
#'   * Condition Number (max / min eigenvalue ratio)
#'   * Max Condition Index (sqrt(max eigenvalue / eigenvalue_i))
#' Both indices help diagnose model-level multicollinearity.
#'
#' @param fitted_model     A \code{lavaan} model object.
#' @param digits           Number of digits for printed estimates.
#' @param standardized     If \code{TRUE}, use fully standardized coefficients (\code{std.all})
#'                         for labels and widths; otherwise raw estimates.
#' @param alpha            p-value threshold for edge opacity.
#' @param low_alpha        Transparency for non-significant edges and \code{NA} p-values.
#' @param min_width        Minimum edge width.
#' @param max_width        Maximum edge width.
#' @param pos_color        Color for positive edges.
#' @param neg_color        Color for negative edges.
#' @param fontname         Font family for all labels and node text.
#' @param node_fontsize    Base font size for nodes.
#' @param edge_fontsize    Base font size for edge labels.
#' @param show_residuals   Whether to draw residual-variance edges.
#' @param show_intercepts  Whether to draw intercept edges.
#' @param show_fit         Whether to include fit indices label at the top.
#' @param show_collinearity Whether to show Condition Number & Max Condition Index
#'                         (default = \code{TRUE}).
#' @param layout           Graph orientation: \code{"LR"} (left→right) or \code{"TB"} (top→bottom).
#' @param ratio            Graphviz \code{ratio} attribute (\code{"fill"} or numeric).
#' @param curvature        Curvature for bidirectional edges.
#' @param engine           Graphviz layout engine (\code{"dot"}, \code{"neato"}, \code{"fdp"},
#'                         \code{"circo"}, \code{"twopi"}).
#' @param twopi_compact    If \code{TRUE} compacts \code{twopi} by shrinking
#'                         \code{ranksep}, \code{nodesep}, and \code{sep}.
#'
#' @return A DiagrammeR \code{htmlwidget} displaying the path diagram.
#' @export
semDiagram <- function(
    fitted_model,
    digits            = 3,
    standardized      = TRUE,
    alpha             = 0.05,
    low_alpha         = 0.2,
    min_width         = 1,
    max_width         = 5,
    pos_color         = "blue",
    neg_color         = "red",
    fontname          = "Helvetica",
    node_fontsize     = 11,
    edge_fontsize     = 9,
    show_residuals    = FALSE,
    show_intercepts   = FALSE,
    show_fit          = TRUE,
    show_collinearity = TRUE,      # ← ★ 新規追加：多重共線性指標の表示スイッチ
    layout            = "LR",
    ratio             = "fill",
    curvature         = 0.3,
    engine            = "dot",
    twopi_compact     = TRUE) {

  ## ---------- 1. Validate arguments & helper ----------

  # Ensure Graphviz engine is valid
  engine <- match.arg(engine, c("dot","neato","fdp","circo","twopi"))

  # Simple RGBA helper that adds alpha to base color
  alpha_color <- function(col, alpha_val) {
    if (is.na(alpha_val) || alpha_val < 0 || alpha_val > 1) alpha_val <- 1
    rgb_mat <- grDevices::col2rgb(col) / 255
    grDevices::rgb(rgb_mat[1,], rgb_mat[2,], rgb_mat[3,], alpha = alpha_val)
  }

  # Check for lavaan object
  if (!inherits(fitted_model, "lavaan"))
    stop("`fitted_model` must be a lavaan object.")

  # Require DiagrammeR & lavaan
  invisible(lapply(c("DiagrammeR","lavaan"), function(p)
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", p))))

  ## ---------- 2. lavaan statistics & multicollinearity diagnostics ----------

  # 2-1) Basic parameter table and fit indices
  params <- lavaan::parameterEstimates(fitted_model, standardized = TRUE)
  scale_col <- if (standardized) "std.all" else "est"

  fit_measures <- lavaan::fitMeasures(
    fitted_model,
    c("pvalue","srmr","rmsea","gfi","agfi","nfi","cfi","aic","bic"))
  n_obs <- lavaan::lavInspect(fitted_model, "nobs")

  # 2-2) --- NEW --- Condition Number & Max Condition Index -------------------
  # Extract sample covariance matrix via lavInspect() (symmetric positive-definite) :contentReference[oaicite:2]{index=2}
  samp_cov <- lavaan::lavInspect(fitted_model, "sampstat")$cov

  # Convert to correlation matrix to remove scale effects
  samp_cor <- stats::cov2cor(samp_cov)  # base R function cov2cor() :contentReference[oaicite:3]{index=3}

  # Eigen-decomposition (symmetric matrix ⇒ numeric stability)
  eig_vals <- eigen(samp_cor, symmetric = TRUE, only.values = TRUE)$values  # base::eigen() :contentReference[oaicite:4]{index=4}

  # Condition Number = max(eigen) / min(eigen)
  condition_number <- max(eig_vals) / min(eig_vals)

  # Condition Indices vector = sqrt(max(eigen) / eigen_i)
  condition_indices <- sqrt(max(eig_vals) / eig_vals)

  # Max Condition Index
  max_cond_index <- max(condition_indices)

  ## ---------- 3. Edge-width scaling (as in original code) ----------

  edge_rows <- params[params$op %in% c("=~","~","~~") & params$lhs != params$rhs, ]
  vals <- abs(edge_rows[[scale_col]]); vals <- vals[!is.na(vals)]
  if (standardized) {
    max_abs <- 1                                   # fixed scaling for std coeffs
  } else {
    max_abs <- if (length(vals) == 0 || !is.finite(max(vals))) 1 else max(vals)
  }

  ## ---------- 4. Helper: colorize according to threshold ----------

  colorize_thresh <- function(v, thr, invert = FALSE) {
    if (is.na(v)) "gray50"                                 # NA → neutral gray
    else if (invert) {                                     # invert: higher is worse
      if (v > thr) "red" else "gray20"
    } else {                                               # normal: lower is worse
      if (v <= thr) "red" else "gray20"
    }
  }

  ## ---------- 5. Top-label construction (fit indices + new collinearity) ----------

  # 5-1) Fit indices block
  fit_block <- if (show_fit) {
    paste0(
      sprintf("N = %d | ", n_obs),
      # p-value
      sprintf("<font color='%s'>p = %.3f</font> | ",
              colorize_thresh(fit_measures["pvalue"], 0.05), fit_measures["pvalue"]),
      # SRMR
      sprintf("<font color='%s'>SRMR = %.3f</font> | ",
              colorize_thresh(fit_measures["srmr"], 0.08, invert = TRUE), fit_measures["srmr"]),
      # RMSEA
      sprintf("<font color='%s'>RMSEA = %.3f</font> | ",
              colorize_thresh(fit_measures["rmsea"], 0.08, invert = TRUE), fit_measures["rmsea"]),
      # AIC / BIC
      sprintf("AIC = %.1f | BIC = %.1f<BR/>", fit_measures["aic"], fit_measures["bic"]),
      # GFI
      sprintf("<font color='%s'>GFI = %.3f</font> | ",
              colorize_thresh(fit_measures["gfi"], 0.90), fit_measures["gfi"]),
      # AGFI
      sprintf("<font color='%s'>AGFI = %.3f</font> | ",
              colorize_thresh(fit_measures["agfi"], 0.90), fit_measures["agfi"]),
      # NFI
      sprintf("<font color='%s'>NFI = %.3f</font> | ",
              colorize_thresh(fit_measures["nfi"], 0.90), fit_measures["nfi"]),
      # CFI
      sprintf("<font color='%s'>CFI = %.3f</font>",
              colorize_thresh(fit_measures["cfi"], 0.90), fit_measures["cfi"])
    )
  } else ""

  # 5-2) Collinearity block (new; optional)
  coll_block <- if (show_collinearity) {
    paste0(
      sprintf("<font color='%s'>Condition Number = %.1f</font>  | ",
              colorize_thresh(condition_number, 30, invert = TRUE), condition_number),
      sprintf("<font color='%s'>Max Condition Index = %.1f</font>",
              colorize_thresh(max_cond_index, 30, invert = TRUE), max_cond_index)
    )
  } else ""

  # 5-3) Coefficient / intercept annotation
  coeff_text <- if (standardized)
    "Coefficients: <b>Standardized</b>"
  else
    "Coefficients: <b>Unstandardized</b>"
  intercept_text <- if (show_intercepts)
    "Intercepts: <b>Shown</b>"
  else
    "Intercepts: <b>Hidden</b>"

  annot_block <- sprintf("%s   | %s", coeff_text, intercept_text)

  # 5-4) Assemble <HTML-like> label
  label_parts <- c(annot_block, if (show_fit) fit_block else NULL,
                   if (show_collinearity) coll_block else NULL)
  # Join with <BR/> where appropriate
  top_label <- sprintf("<%s>", paste(label_parts, collapse = "<BR/>"))

  ## ---------- 6. Nodes ----------

  latents   <- unique(params$lhs[params$op == "=~"])
  observeds <- setdiff(unique(c(params$lhs, params$rhs)),
                       c(latents, "1", ""))
  nodes <- list()
  # Latent variable nodes (ellipses)
  for (lv in latents) nodes[[lv]] <- list(
    shape = "ellipse",
    label = lv,
    style = "filled",
    fillcolor = "#F0F0F0",
    fontname = fontname,
    fontsize = node_fontsize)
  # Observed variable nodes (boxes)
  for (ov in observeds) nodes[[ov]] <- list(
    shape = "box",
    label = ov,
    fontname = fontname,
    fontsize = node_fontsize)

  ## ---------- 7. Edges ----------

  edges <- list()
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    if (p$op %in% c("=~","~","~~") && p$lhs != p$rhs) {

      # 7-1) Determine value to plot
      value <- if (standardized) p$std.all else p$est
      if (is.na(value)) value <- p$est                     # fallback to raw

      # 7-2) Pen width scaling
      pen <- (abs(value) / max_abs) * (max_width - min_width) + min_width
      if (!is.finite(pen)) pen <- min_width

      # 7-3) Edge transparency based on p-value
      alpha_edge <- if (is.na(p$pvalue)) low_alpha
      else if (p$pvalue < alpha) 1
      else low_alpha

      # 7-4) Edge color (positive/negative) with alpha
      col <- alpha_color(if (value >= 0) pos_color else neg_color, alpha_edge)

      # 7-5) Build base edge definition
      e_base <- switch(p$op,
                       "=~" = list(from = p$lhs, to = p$rhs, arrowhead = "vee"),
                       "~"  = list(from = p$rhs, to = p$lhs, arrowhead = "vee"),
                       "~~" = list(from = p$lhs, to = p$rhs,
                                   arrowhead = "vee", arrowtail = "vee",
                                   dir = "both", style = "dashed"))
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

  ## ---------- 8. DOT graph construction ----------

  # 8-1) Node definitions
  node_defs <- paste(
    vapply(names(nodes), function(n) {
      attrs <- paste(
        names(nodes[[n]]),
        vapply(nodes[[n]], function(x)
          if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x),
          character(1)),
        sep = "=", collapse = ", ")
      sprintf("  \"%s\" [%s];", n, attrs)
    }, character(1)), collapse = "\n")

  # 8-2) Edge definitions
  edge_defs <- paste(
    vapply(edges, function(e) {
      attrs <- paste(
        names(e)[-1:-2],
        vapply(e[-1:-2], function(x)
          if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x),
          character(1)),
        sep = "=", collapse = ", ")
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }, character(1)), collapse = "\n")

  # 8-3) Graph-level options (Radial tweaks for circo / twopi)
  radial_opts <- if (engine == "circo") {
    c("splines=true", "nodesep=0.4", "sep=\"+4\"", "mindist=1")
  } else if (engine == "twopi" && twopi_compact) {
    c("splines=true", "nodesep=0.2", "sep=\"+2\"",
      "ranksep=0.5", "normalize=true")
  } else character(0)
  radial_opts <- paste(radial_opts, collapse = ", ")

  # 8-4) Assemble DOT source
  graph_code <- sprintf(
    "digraph {\n  rankdir=%s;\n  graph [layout=%s%s%s, overlap=false,\n         labelloc=\"t\", labeljust=\"c\", label=%s, ratio=%s];\n  node  [fontname=\"%s\", margin=0.05];\n  edge  [fontname=\"%s\", fontcolor=\"#333333\"];\n\n%s\n\n%s\n}",
    layout, engine,
    if (nchar(radial_opts)) ", " else "", radial_opts,
    top_label, ratio,
    fontname, fontname,
    node_defs, edge_defs)

  ## ---------- 9. Render diagram ----------

  DiagrammeR::grViz(graph_code, engine = engine)
}
