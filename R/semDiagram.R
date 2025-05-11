#' Visualize SEM Results as a Path Diagram (no external alpha dependency; robust against NA)
#'
#' @param fitted_model    A lavaan model object
#' @param digits          Number of digits for printed estimates
#' @param standardized    If TRUE, use fully standardized coefficients (std.all) for labels & widths; otherwise raw estimates
#' @param alpha           p-value threshold for edge opacity
#' @param low_alpha       Transparency for non-significant edges and NA p-values
#' @param min_width       Minimum edge width
#' @param max_width       Maximum edge width
#' @param pos_color       Color for positive edges
#' @param neg_color       Color for negative edges
#' @param fontname        Font family for all labels
#' @param node_fontsize   Base font size for nodes
#' @param edge_fontsize   Base font size for edge labels
#' @param show_residuals  Whether to draw residual‐variance edges
#' @param show_intercepts Whether to draw intercept edges
#' @param show_fit        Whether to include fit indices label at the top
#' @param layout          "LR" (left→right) or "TB" (top→bottom)
#' @param ratio           GraphViz 'ratio' attribute ("fill" or numeric)
#' @param curvature       Curvature for bidirectional edges
#' @param engine          Graphviz layout engine ("dot","neato","fdp","circo","twopi")
#' @param twopi_compact   If TRUE (default) compacts twopi by shrinking ranksep/nodesep/sep
#' @export
semDiagram <- function(fitted_model,
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

  # --- validate layout engine ------------------------------------
  engine <- match.arg(engine,
                      choices = c("dot","neato","fdp","circo","twopi"))

  # Helper: RGBA color string
  alpha_color <- function(col, alpha_val) {
    if (is.na(alpha_val) || alpha_val < 0 || alpha_val > 1) alpha_val <- 1
    rgb_mat <- grDevices::col2rgb(col) / 255
    grDevices::rgb(
      red   = rgb_mat[1, ],
      green = rgb_mat[2, ],
      blue  = rgb_mat[3, ],
      alpha = alpha_val
    )
  }

  # 1) Validate model
  if (!inherits(fitted_model, "lavaan")) {
    stop("`fitted_model` must be a lavaan object.")
  }

  # 2) Required packages
  req_pkgs <- c("DiagrammeR", "lavaan")
  invisible(lapply(req_pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed.", pkg))
    }
  }))

  # 3) Parameter estimates
  params    <- lavaan::parameterEstimates(fitted_model, standardized = TRUE)
  scale_col <- if (standardized) "std.all" else "est"

  # 4) Fit indices & N
  fit   <- lavaan::fitMeasures(fitted_model,
                               c("pvalue","srmr","rmsea","gfi","agfi","nfi","cfi","aic","bic"))
  n_obs <- lavaan::lavInspect(fitted_model, "nobs")

  # 5) Width scaling
  edge_rows <- params[params$op %in% c("=~","~","~~") & params$lhs != params$rhs, ]
  vals      <- abs(edge_rows[[scale_col]])
  vals      <- vals[!is.na(vals)]
  max_abs   <- if (length(vals)==0 || !is.finite(max(vals))) 1 else max(vals)

  # 6) Fit index color helper
  colorize_fit <- function(v, thr, inv=FALSE) {
    if (is.na(v)) return("gray50")
    if (inv) { if (v > thr) "red" else "gray20" }
    else     { if (v <= thr) "red" else "gray20" }
  }

  # 7) Fit label
  if (show_fit) {
    p_col  <- colorize_fit(fit["pvalue"], 0.05)
    s_col  <- colorize_fit(fit["srmr"],  0.08, TRUE)
    r_col  <- colorize_fit(fit["rmsea"], 0.08, TRUE)
    g_col  <- colorize_fit(fit["gfi"],   0.90)
    ag_col <- colorize_fit(fit["agfi"],  0.90)
    n_col  <- colorize_fit(fit["nfi"],   0.90)
    c_col  <- colorize_fit(fit["cfi"],   0.90)

    fit_label <- paste0(
      "<<font point-size='12'>",
      sprintf("N = %d", n_obs), " | ",
      sprintf("<font color='%s'>p = %.3f</font>",   p_col,  fit["pvalue"]), " | ",
      sprintf("<font color='%s'>SRMR = %.3f</font>",s_col,  fit["srmr"]),   " | ",
      sprintf("<font color='%s'>RMSEA = %.3f</font>",r_col, fit["rmsea"]),  " | ",
      sprintf("AIC = %.1f", fit["aic"]),             " | ",
      sprintf("BIC = %.1f", fit["bic"]),             "<br/>",
      sprintf("<font color='%s'>GFI = %.3f</font>",  g_col,  fit["gfi"]),   " | ",
      sprintf("<font color='%s'>AGFI = %.3f</font>", ag_col, fit["agfi"]),  " | ",
      sprintf("<font color='%s'>NFI = %.3f</font>",  n_col,  fit["nfi"]),   " | ",
      sprintf("<font color='%s'>CFI = %.3f</font>",  c_col,  fit["cfi"]),
      "</font>>"
    )
  } else fit_label <- ""

  # 8) Nodes
  latents   <- unique(params$lhs[params$op == "=~"])
  observeds <- setdiff(unique(c(params$lhs, params$rhs)), c(latents, "1", ""))
  nodes     <- list()
  for (lv in latents) {
    nodes[[lv]] <- list(shape="ellipse", label=lv, style="filled",
                        fillcolor="#F0F0F0", fontname=fontname, fontsize=node_fontsize)
  }
  for (ov in observeds) {
    nodes[[ov]] <- list(shape="box", label=ov,
                        fontname=fontname, fontsize=node_fontsize)
  }

  # 9) Edges
  edges <- list()
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    if (p$op %in% c("=~","~","~~") && p$lhs != p$rhs) {
      value <- if (standardized) p$std.all else p$est
      if (is.na(value)) value <- p$est
      pen <- (abs(value)/max_abs)*(max_width-min_width)+min_width
      if (!is.finite(pen)) pen <- min_width

      alpha_edge <- if (is.na(p$pvalue)) low_alpha else
        if (p$pvalue < alpha) 1 else low_alpha

      base_col <- if (value >= 0) pos_color else neg_color
      col      <- alpha_color(base_col, alpha_edge)

      e_base <- switch(p$op,
                       "=~" = list(from=p$lhs, to=p$rhs, arrowhead="vee"),
                       "~"  = list(from=p$rhs, to=p$lhs, arrowhead="vee"),
                       "~~" = list(from=p$lhs, to=p$rhs, arrowhead="vee", arrowtail="vee",
                                   dir="both", style="dashed")
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
      edges[[length(edges)+1]] <- e_base
    }
  }

  # 10) Residuals / intercepts (未実装部分；必要なら追加)

  # 11) Assemble DOT
  node_defs <- paste(
    vapply(names(nodes), function(n) {
      attrs <- paste(
        names(nodes[[n]]),
        vapply(nodes[[n]], function(x)
          if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x), character(1)),
        sep="=", collapse=", "
      )
      sprintf("  \"%s\" [%s];", n, attrs)
    }, character(1)), collapse="\n"
  )

  edge_defs <- paste(
    vapply(edges, function(e) {
      attrs <- paste(
        names(e)[-1:-2],
        vapply(e[-1:-2], function(x)
          if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x), character(1)),
        sep="=", collapse=", "
      )
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }, character(1)), collapse="\n"
  )

  # --- radial engine specific options ----------------------------
  radial_opts <- if (engine == "circo") {
    c("splines=true",
      "nodesep=0.4",
      "sep=\"+4\"",
      "mindist=1")
  } else if (engine == "twopi" && twopi_compact) {
    c("splines=true",
      "nodesep=0.2",
      "sep=\"+2\"",
      "ranksep=0.5",
      "normalize=true")
  } else character(0)
  radial_opts <- paste(radial_opts, collapse=", ")

  graph_code <- sprintf(
    "digraph {\n  rankdir=%s;\n  graph [layout=%s%s%s, overlap=false, labelloc=\"t\", labeljust=\"c\", label=%s, ratio=%s];\n  node [fontname=\"%s\", margin=0.05];\n  edge [fontname=\"%s\", fontcolor=\"#333333\"];\n\n%s\n\n%s\n}",
    layout,
    engine,
    if (nchar(radial_opts)) ", " else "",
    radial_opts,
    fit_label,
    ratio,
    fontname,
    fontname,
    node_defs,
    edge_defs
  )

  # 12) Render
  DiagrammeR::grViz(graph_code, engine=engine)
}
