semDiagram <- function(fitted_model,
                       digits = 3,
                       standardized = TRUE,
                       alpha = 0.05,
                       min_width = 1, max_width = 5,
                       pos_color = "blue", neg_color = "red",
                       fontname = "Helvetica",
                       node_fontsize = 11,
                       edge_fontsize = 9,
                       show_residuals = FALSE,
                       show_intercepts = FALSE,
                       show_fit = TRUE,
                       layout = "LR",
                       curvature = 0.3) {

  # パッケージチェックとインストール
  if (!inherits(fitted_model, "lavaan")) {
    stop("The argument 'fitted_model' must be a lavaan model object.")
  }

  required_pkgs <- c("scales", "DiagrammeR", "lavaan")
  missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(missing_pkgs) > 0) {
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs)
  }
  invisible(lapply(required_pkgs, library, character.only = TRUE))

  # パラメータ推定値の取得
  params <- parameterEstimates(fitted_model, standardized = standardized)

  # 適合度指標の取得
  fit <- fitMeasures(fitted_model,
                     c("pvalue", "srmr", "rmsea", "gfi", "agfi", "nfi", "cfi", "aic", "bic"))

  # 適合度指標の色付け関数
  colorize <- function(value, threshold, inverse = FALSE) {
    if (is.na(value)) return("gray50")
    if (inverse) {
      if (value > threshold) "red" else "gray20"
    } else {
      if (value <= threshold) "red" else "gray20"
    }
  }

  # 適合度ラベルの作成
  if (show_fit) {
    p_col    <- colorize(fit["pvalue"], 0.05)
    s_col    <- colorize(fit["srmr"], 0.08, inverse = TRUE)
    r_col    <- colorize(fit["rmsea"], 0.08, inverse = TRUE)
    gfi_col  <- colorize(fit["gfi"], 0.90)
    agfi_col <- colorize(fit["agfi"], 0.90)
    nfi_col  <- colorize(fit["nfi"], 0.90)
    cfi_col  <- colorize(fit["cfi"], 0.90)

    fit_label <- paste0(
      "<<font point-size='12'>",
      sprintf("<font color='%s'>p = %.3f</font>", p_col, fit["pvalue"]), " | ",
      sprintf("<font color='%s'>SRMR = %.3f</font>", s_col, fit["srmr"]), " | ",
      sprintf("<font color='%s'>RMSEA = %.3f</font>", r_col, fit["rmsea"]), " | ",
      sprintf("AIC = %.1f", fit["aic"]), " | ",
      sprintf("BIC = %.1f", fit["bic"]), "<br/>",
      sprintf("<font color='%s'>GFI = %.3f</font>", gfi_col, fit["gfi"]), " | ",
      sprintf("<font color='%s'>AGFI = %.3f</font>", agfi_col, fit["agfi"]), " | ",
      sprintf("<font color='%s'>NFI = %.3f</font>", nfi_col, fit["nfi"]), " | ",
      sprintf("<font color='%s'>CFI = %.3f</font>", cfi_col, fit["cfi"]),
      "</font>>"
    )
  } else {
    fit_label <- ""
  }

  # ノードの定義
  # まず、各変数名を文字列として扱う
  latent_vars <- unique(as.character(params$lhs[params$op == "=~"]))
  observed_vars <- unique(as.character(c(params$lhs, params$rhs)))

  # "1"（およびもしあれば空文字列）や latent_vars を除外する
  observed_vars <- setdiff(observed_vars, c(latent_vars, "1", ""))

  nodes <- list()

  # 潜在変数
  for (lv in latent_vars) {
    nodes[[lv]] <- list(
      shape = "ellipse",
      label = lv,
      style = "filled",
      fillcolor = "#F0F0F0",
      fontname = fontname,
      fontsize = node_fontsize
    )
  }

  # 観測変数
  for (ov in observed_vars) {
    nodes[[ov]] <- list(
      shape = "box",
      label = ov,
      fontname = fontname,
      fontsize = node_fontsize
    )
  }

  # エッジの定義
  edges <- list()

  for (i in seq_len(nrow(params))) {
    p <- params[i, ]

    if (p$op %in% c("=~", "~", "~~") && p$lhs != p$rhs) {
      # エッジの太さ (標準化係数の絶対値に基づく)
      pen <- abs(p$est) * (max_width - min_width) + min_width

      # 透明度 (有意性に基づく)
      alpha_edge <- ifelse(p$pvalue < alpha, 1, 0.3)

      # 色 (係数の符号に基づく)
      col <- scales::alpha(ifelse(p$est >= 0, pos_color, neg_color), alpha_edge)

      # エッジタイプごとの設定
      edge_def <- switch(
        p$op,
        "=~" = list(from = p$lhs, to = p$rhs, arrowhead = "vee"),
        "~"  = list(from = p$rhs, to = p$lhs, arrowhead = "vee"),
        "~~" = list(from = p$lhs, to = p$rhs,
                    arrowhead = "vee", arrowtail = "vee",
                    dir = "both", style = "dashed")
      )

      # 共通設定
      edge_def$label <- sprintf("%.*f", digits, p$est)
      edge_def$penwidth <- pen
      edge_def$color <- col
      edge_def$fontsize <- edge_fontsize
      edge_def$fontname <- fontname

      # 曲率設定 (双方向関係の場合)
      if (p$op == "~~") {
        edge_def$constraint <- FALSE
        edge_def$dir <- "both"
      }

      edges[[length(edges) + 1]] <- edge_def
    }
  }

  # 残差分散の表示
  if (show_residuals) {
    resids <- params[params$op == "~~" & params$lhs == params$rhs, ]
    for (i in seq_len(nrow(resids))) {
      r <- resids[i, ]
      node_id <- paste0("res_", r$lhs)

      # 残差ノードの追加
      nodes[[node_id]] <- list(
        shape = "circle",
        label = paste0("ε_", r$lhs),
        fontname = fontname,
        fontsize = node_fontsize - 1,
        width = 0.5,
        height = 0.5
      )

      # 残差エッジの追加
      edges[[length(edges) + 1]] <- list(
        from = node_id,
        to = r$lhs,
        label = sprintf("%.*f", digits, r$est),
        arrowhead = "vee",
        style = "dotted",
        penwidth = 1,
        color = "gray50",
        fontsize = edge_fontsize
      )
    }
  }

  # 切片の表示
  # 切片の表示（個別表示版: ノードは定型ラベル「1」、エッジに推定値を表示）
  if (show_intercepts) {
    ints <- params[params$op == "~1", ]
    if (nrow(ints) > 0) {
      for (i in seq_len(nrow(ints))) {
        p <- ints[i, ]
        # 各切片ノードのIDはユニークになるように "int_" + 変数名
        node_id <- paste0("int_", p$lhs)

        # 個別の切片ノードを作成（ラベルは定型ラベル "1"）
        nodes[[node_id]] <- list(
          shape = "triangle",
          label = "1",
          style = "filled",
          fillcolor = "#DDDDDD",
          fontname = fontname,
          fontsize = node_fontsize * 0.5,
          width = 0,
          height = 0
        )

        # 切片ノードから対象変数ノードへのエッジを追加（エッジラベルに切片の数値を表示）
        edges[[length(edges) + 1]] <- list(
          from = node_id,
          to = p$lhs,
          label = sprintf("%.*f", digits, p$est),
          arrowhead = "vee",
          style = "solid",
          penwidth = 1,
          color = "gray50",
          fontsize = edge_fontsize,
          fontname = fontname
        )
      }
    }
  }


  # GraphVizコードの生成
  graph_code <- paste0(
    "digraph {\n",
    "  rankdir = ", layout, ";\n",
    "  graph [overlap = false, fontsize = ", node_fontsize,
    ", labelloc = 't', labeljust = 'c', label = ", fit_label, "];\n",
    "  node [fontname = '", fontname, "', margin = 0.05];\n",
    "  edge [fontname = '", fontname, "', fontcolor = '#333333'];\n\n",

    paste(sapply(names(nodes), function(n) {
      attrs <- paste(
        names(nodes[[n]]),
        sapply(nodes[[n]], function(x) {
          if (is.numeric(x)) x else paste0("\"", x, "\"")
        }),
        sep = "=",
        collapse = ", "
      )
      sprintf("  \"%s\" [%s];", n, attrs)
    }), collapse = "\n"), "\n\n",

    paste(sapply(edges, function(e) {
      attrs <- paste(names(e)[-1:-2],
                     sapply(e[-1:-2], function(x) {
                       if (is.numeric(x)) x else paste0("\"", x, "\"")
                     }),
                     sep = "=", collapse = ", ")
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }), collapse = "\n"), "\n",
    "}\n"
  )

  # グラフの描画
  DiagrammeR::grViz(graph_code)
}
