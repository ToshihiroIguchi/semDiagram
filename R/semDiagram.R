#' Visualize SEM Results as a Path Diagram (robust penwidth; no blank graph)
#'
#' @param fitted_model An object of class lavaan
#' @param digits Integer: number of digits for estimates
#' @param standardized Logical: If TRUE, both labels *and* edge widths use fully-standardized coefficients (std.all); otherwise raw estimates.
#' @param alpha Numeric: significance threshold for edge transparency
#' @param min_width,max_width Numeric: minimum / maximum edge width
#' @param pos_color,neg_color Character: colors for positive/negative estimates
#' @param fontname Character: font family for labels
#' @param node_fontsize,edge_fontsize Numeric: font sizes
#' @param show_residuals,show_intercepts,show_fit Logical switches
#' @param layout Character: "LR" (left-to-right) or "TB" (top-to-bottom)
#' @param ratio Character or numeric: GraphViz 'ratio' attribute
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

  # 1. sanity check
  if (!inherits(fitted_model, "lavaan"))
    stop("`fitted_model` must be a lavaan model object.")

  # 2. deps
  pkgs <- c("scales", "DiagrammeR", "lavaan")
  lapply(pkgs, function(pkg){
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' is required but not installed.", pkg))
  })

  # 3. estimates (always include standardized cols)
  params <- lavaan::parameterEstimates(fitted_model, standardized = TRUE)  # always TRUE
  scale_col <- if (standardized) "std.all" else "est"

  # 4. fit indices
  fit   <- lavaan::fitMeasures(fitted_model,
                               c("pvalue","srmr","rmsea","gfi","agfi","nfi","cfi","aic","bic"))
  n_obs <- lavaan::lavInspect(fitted_model, "nobs")

  # ---------- penwidth scaling ----------
  edge_rows <- params[params$op %in% c("=~","~","~~") & params$lhs != params$rhs, ]
  vals <- abs(edge_rows[[scale_col]])
  vals <- vals[!is.na(vals)]
  max_abs <- ifelse(length(vals) == 0 || !is.finite(max(vals)), 1, max(vals))
  # --------------------------------------

  # 5. fit label helper
  colorize <- function(v,t,inv=FALSE){
    if (is.na(v)) return("gray50")
    if (inv) if (v>t) "red" else "gray20" else if (v<=t) "red" else "gray20"
  }

  if (show_fit){
    p_col   <- colorize(fit["pvalue"], 0.05)
    s_col   <- colorize(fit["srmr"],   0.08, TRUE)
    r_col   <- colorize(fit["rmsea"],  0.08, TRUE)
    g_col   <- colorize(fit["gfi"],    0.90)
    ag_col  <- colorize(fit["agfi"],   0.90)
    n_col   <- colorize(fit["nfi"],    0.90)
    c_col   <- colorize(fit["cfi"],    0.90)
    fit_label <- paste0(
      "<<font point-size='12'>",
      sprintf("N = %d", n_obs)," | ",
      sprintf("<font color='%s'>p = %.3f</font>",   p_col, fit["pvalue"])," | ",
      sprintf("<font color='%s'>SRMR = %.3f</font>",s_col, fit["srmr"])," | ",
      sprintf("<font color='%s'>RMSEA = %.3f</font>",r_col,fit["rmsea"])," | ",
      sprintf("AIC = %.1f", fit["aic"])," | ",
      sprintf("BIC = %.1f", fit["bic"]),"<br/>",
      sprintf("<font color='%s'>GFI = %.3f</font>",  g_col, fit["gfi"])," | ",
      sprintf("<font color='%s'>AGFI = %.3f</font>", ag_col,fit["agfi"])," | ",
      sprintf("<font color='%s'>NFI = %.3f</font>",  n_col, fit["nfi"])," | ",
      sprintf("<font color='%s'>CFI = %.3f</font>",  c_col, fit["cfi"]),
      "</font>>")
  } else fit_label <- ""

  # 6. nodes
  latents   <- unique(params$lhs[params$op=="=~"])
  observeds <- setdiff(unique(c(params$lhs,params$rhs)), c(latents,"1",""))
  nodes <- lapply(latents, function(lv) list(
    shape="ellipse",label=lv,style="filled",fillcolor="#F0F0F0",
    fontname=fontname,fontsize=node_fontsize))
  names(nodes) <- latents
  for (ov in observeds) nodes[[ov]] <- list(
    shape="box",label=ov,fontname=fontname,fontsize=node_fontsize)

  # 7. edges
  edges <- list()
  for(i in seq_len(nrow(params))){
    p <- params[i,]
    if (p$op %in% c("=~","~","~~") && p$lhs!=p$rhs){
      value <- if (standardized) p$std.all else p$est
      if (is.na(value)) value <- p$est
      pen   <- (abs(value)/max_abs)*(max_width-min_width)+min_width
      if (is.na(pen) || !is.finite(pen)) pen <- min_width
      alpha <- ifelse(p$pvalue < alpha,1,0.3)
      col   <- scales::alpha(ifelse(value>=0,pos_color,neg_color),alpha)
      e <- switch(p$op,
                  "=~"=list(from=p$lhs,to=p$rhs,arrowhead="vee"),
                  "~" =list(from=p$rhs,to=p$lhs,arrowhead="vee"),
                  "~~"=list(from=p$lhs,to=p$rhs,arrowhead="vee",arrowtail="vee",
                            dir="both",style="dashed"))
      e$label    <- sprintf("%.*f", digits, value)
      e$penwidth <- pen
      e$color    <- col
      e$fontsize <- edge_fontsize
      e$fontname <- fontname
      if (p$op=="~~"){ e$constraint <- FALSE; e$dir <- "both" }
      edges[[length(edges)+1]] <- e
    }
  }

  # 8. residuals / intercepts （オリジナルと同じロジック。省略せず全文コピペする場合はここに挿入してください）
  # ------------------------------------------------------------------

  # 9. GraphViz dot
  node_str <- paste(
    vapply(names(nodes), function(n){
      a <- nodes[[n]]
      attrs <- paste(names(a),
                     vapply(a, function(x)
                       if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x),
                       character(1)),
                     sep="=",collapse=", ")
      sprintf("  \"%s\" [%s];", n, attrs)
    }, character(1)), collapse="\n")

  edge_str <- paste(
    vapply(edges, function(e){
      attrs <- paste(names(e)[-1:-2],
                     vapply(e[-1:-2], function(x)
                       if (is.numeric(x)) as.character(x) else sprintf("\"%s\"", x),
                       character(1)),
                     sep="=",collapse=", ")
      sprintf("  \"%s\" -> \"%s\" [%s];", e$from, e$to, attrs)
    }, character(1)), collapse="\n")

  graph_code <- sprintf(
    "digraph {\n  rankdir=%s;\n  graph [overlap=false, labelloc='t', labeljust='c', label=%s, ratio=%s];\n  node [fontname=\"%s\", margin=0.05];\n  edge [fontname=\"%s\", fontcolor=\"#333333\"];\n\n%s\n\n%s\n}",
    layout, fit_label, ratio, fontname, fontname, node_str, edge_str)

  # 10. draw
  DiagrammeR::grViz(graph_code)
}
