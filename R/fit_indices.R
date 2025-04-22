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
