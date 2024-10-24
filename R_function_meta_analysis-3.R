
library(metafor)

## res.studies: list of results of differential expression analysis
## freq.studies: only genes available in freq.studies of the studies will be used
## res.file: text file for saving results of fixed and random effect meta analysis
#Rückgabe des rma.uni des genes mit kleinstem P-Wert
run_meta_analysis_3 <- function(res.studies,
                              freq.studies = 0.5,
                              res.file) {
  
  ## set rownames
  res.studies = lapply(res.studies, function(x) {
    rownames(x) = x$gene
    return(x)
  })
  
  ## extract genes to be used
  info.genes.name = unlist(lapply(res.studies, rownames))
  tab.name = table(info.genes.name)
  genes = sort(names(tab.name)[tab.name >= freq.studies * length(res.studies)])
  
  ## extract genes to be used
  #Alte funktion
#  info.genes = unlist(lapply(res.studies, rownames))
#  tab = table(info.genes)
#  genes = sort(names(tab)[tab >= freq.studies * length(res.studies)])
  
  ## gene level meta analysis
  res.all = NULL
  for (g in genes) {
    dat = na.omit(data.frame(
      estimate = sapply(res.studies, function(x) {x[g, "logFC"]}),
      se = sapply(res.studies, function(x) {x[g, "SE"]}),
      pvalue = sapply(res.studies, function(x) {x[g, "P.Value"]})))
    
    res.fixed = rma.uni(yi = dat$estimate,
                        sei = dat$se,
                        method = "FE")
    res.random = rma.uni(yi = dat$estimate,
                         sei = dat$se,
                         method = "DL")
    res.g = c(
      number.studies = res.fixed$k,
      estimate.fixed = res.fixed$beta,
      ci.l.fixed = res.fixed$ci.lb,
      ci.u.fixed = res.fixed$ci.ub,
      z.fixed = res.fixed$zval,
      pvalue.fixed = res.fixed$pval,
      estimate.random = res.random$beta,
      ci.l.random = res.random$ci.lb,
      ci.u.random = res.random$ci.ub,
      z.random = res.random$zval,
      pvalue.random = res.random$pval,
      tau2 = res.random$tau2,
      I2 = res.random$I2,
      Q = res.random$QE,
      pvalue.Q = res.random$QEp)
    res.all = rbind(res.all, res.g)
    
    #gen mit groester I2
    if("ENSG00000153292" == g){
      res.klpw = res.random}
    
    
  }
  res.all = data.frame(gene = genes,
                       res.all,
                       stringsAsFactors = FALSE)
  
  ## adjust P-values
  res.all$pvalue.fixed.adj = p.adjust(res.all$pvalue.fixed, method = "BH")
  res.all$pvalue.random.adj = p.adjust(res.all$pvalue.random, method = "BH")
  
  export(res.all,
         file = res.file)
  
  return(res.klpw)
  
}




