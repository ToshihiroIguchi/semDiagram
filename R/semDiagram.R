# lavaanの結果から適合度指標を抜き出してlistにする関数
fit_indices <- function(model) {
  fit <- fitMeasures(model)
  list(
    pvalue = fit["pvalue"],
    SRMR   = fit["srmr"],
    RMSEA  = fit["rmsea"],
    AIC    = fit["aic"],
    GFI    = fit["gfi"],
    AGFI   = fit["agfi"],
    NFI    = fit["nfi"],
    CFI    = fit["cfi"]
  )
}

# lavaan と DiagrammeR, scales パッケージを利用した SEM 図と適合度指標の表示関数
semDiagram <- function(fitted_model, digits = 3,
                       min_width = 1, max_width = 5,
                       pos_color = "#1f78b4", neg_color = "#e41a1c") {

  if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
  if (!requireNamespace("DiagrammeR", quietly = TRUE)) install.packages("DiagrammeR")
  library(scales)
  library(DiagrammeR)

  params <- parameterEstimates(fitted_model)
  fit <- fitMeasures(fitted_model)

  # 色を決める関数（条件を満たしたら赤）
  colorize <- function(value, threshold, comp = "<=", inverse = FALSE) {
    if (inverse) {
      if (value > threshold) "red" else "gray20"
    } else {
      if (value <= threshold) "red" else "gray20"
    }
  }

  # 指標ごとの色付け（判定基準は一般的なSEMの基準）
  p_col    <- colorize(fit["pvalue"], 0.05)
  srmr_col <- colorize(fit["srmr"], 0.08, inverse = TRUE)
  rmsea_col<- colorize(fit["rmsea"], 0.08, inverse = TRUE)
  cfi_col  <- colorize(fit["cfi"], 0.90)
  gfi_col  <- colorize(fit["gfi"], 0.90)
  agfi_col <- colorize(fit["agfi"], 0.90)
  nfi_col  <- colorize(fit["nfi"], 0.90)

  # AICは情報量なので色付けしない
  aic_str <- sprintf("AIC = %.1f", fit["aic"])

  # HTMLラベル（Graphviz形式）で横並び＆色付き
  fit_label <- paste0(
    "<",
    sprintf("<FONT COLOR='%s'>p = %.3f</FONT>", p_col, fit["pvalue"]),
    " | ", sprintf("<FONT COLOR='%s'>SRMR = %.3f</FONT>", srmr_col, fit["srmr"]),
    " | ", sprintf("<FONT COLOR='%s'>RMSEA = %.3f</FONT>", rmsea_col, fit["rmsea"]),
    " | ", aic_str,
    " | ", sprintf("<FONT COLOR='%s'>GFI = %.3f</FONT>", gfi_col, fit["gfi"]),
    " | ", sprintf("<FONT COLOR='%s'>AGFI = %.3f</FONT>", agfi_col, fit["agfi"]),
    " | ", sprintf("<FONT COLOR='%s'>NFI = %.3f</FONT>", nfi_col, fit["nfi"]),
    " | ", sprintf("<FONT COLOR='%s'>CFI = %.3f</FONT>", cfi_col, fit["cfi"]),
    ">"
  )

  latent_vars <- unique(params$lhs[params$op == "=~"])
  all_vars <- unique(c(params$lhs, params$rhs))

  # ノード定義
  nodes <- setNames(lapply(all_vars, function(var) {
    list(shape = ifelse(var %in% latent_vars, "box", "ellipse"), label = var)
  }), all_vars)

  # エッジ定義
  edges <- list()
  for (i in seq_len(nrow(params))) {
    p <- params[i, ]
    if (p$op %in% c("=~", "~", "~~") && p$lhs != p$rhs) {
      penwidth <- abs(p$est) * (max_width - min_width) + min_width
      color    <- ifelse(p$est >= 0, pos_color, neg_color)
      edge_def <- switch(p$op,
                         "=~" = list(from = p$lhs, to = p$rhs, arrowhead = "odot"),
                         "~"  = list(from = p$rhs, to = p$lhs, arrowhead = "normal"),
                         "~~" = list(from = p$lhs, to = p$rhs, arrowhead = "none", dir = "both", style = "dashed"))
      edge_def$label    <- sprintf("%.*f", digits, p$est)
      edge_def$penwidth <- penwidth
      edge_def$color    <- color
      edges[[length(edges) + 1]] <- edge_def
    }
  }

  # Graphvizコード生成（HTMLラベル対応）
  graph_code <- paste0(
    "digraph {\n",
    "  rankdir = LR;\n",
    "  graph [overlap=false, fontsize=11, labelloc=\"t\", labeljust=\"c\", label=", fit_label, "];\n",
    "  node [fontname=Helvetica, width=1.5, height=0.8];\n",
    "  edge [fontname=Helvetica, fontcolor='#333333'];\n",
    paste(sapply(names(nodes), function(n) {
      sprintf("  %s [shape=%s, label=\"%s\"];", n, nodes[[n]]$shape, nodes[[n]]$label)
    }), collapse = "\n"), "\n",
    paste(sapply(edges, function(e) {
      attrs <- paste(names(e)[-1:-2], sapply(e[-1:-2], function(x) {
        if (is.numeric(x)) x else paste0("\"", x, "\"")
      }), sep = "=", collapse = ", ")
      sprintf("  %s -> %s [%s];", e$from, e$to, attrs)
    }), collapse = "\n"), "\n",
    "}\n"
  )

  DiagrammeR::grViz(graph_code)
}
