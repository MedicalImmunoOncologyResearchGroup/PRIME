testingKMForDifferentValues <- function(x, min = 0.25, max = 0.75){
  # making a matrix for the results
  quanPval <- data.frame('quan' =  seq(min, max, by = 0.01), 'Pval' = 1)
  # testing and saving the p-values
  for(p in quanPval$quan){
    x <- transform(x, 'exp' = ifelse(x$Exp >= quantile(x$Exp, p), 'high', 'low'))
    try(km.fit <- survfit(Surv(time = as.numeric(OS.time), event = OS) ~ exp, data = x), silent = TRUE)
    quanPval[quanPval$quan == p , 'Pval' ] <- surv_pvalue(km.fit, data = x)[1, 2]
  }
  return(quanPval)
}