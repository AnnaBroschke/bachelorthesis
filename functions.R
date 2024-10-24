 
library(SummarizedExperiment)
library(limma)
library(variancePartition)
library(rio)
library(ggplot2)
library(ggpubr)
library(rio)
library(xtable)
library(ggthemes)

# library(biomaRt)
# # library(GEOquery)
# # library(sva)
# library(tibble)
# library(tidyr)
# library(edgeR)
# # 
# # # #library(testr)
# library(nlme)   ## linear mixed model
# library(data.table) ## fwrite()
# library(metafor)
# library(rio)
# library(ggpubr)
# library(caret)
# # # library(GGally)
# # # library(viridis)
# # # library(metafor)
# # # library(forestplot)
# # # library(plyr)
# # # library(EnsDb.Hsapiens.v75)
# # # library(UpSetR)
# 

## load study data (stored as SE objects) into list
load.se <- function(se.files) {
    se.l = lapply(se.files, 
                  readRDS)
    names = sapply(se.files,
                   function(f) {
                       temp = unlist(strsplit(
                           gsub(".rds", "", basename(f)), "_"))
                       return(temp[1])})
    names(names) = NULL
    names(se.l) = names
    return(se.l)
} 


## perform differential expression analysis for a single study
run_diff_expr_analysis_dream <- function(se,
                                         var, 
                                         covar = NULL,
                                         var.id,
                                         res.file) {
  
  ## define formula for linear mixed model
  form = paste0("~", var)
  if (!is.null(covar)) {
    form = paste0(form, "+",
                  paste(covar, collapse = "+"))
  }
  form.random = paste0("(1|", var.id, ")")
  form.final = paste(form, form.random, sep = " + ")
  
  ## extract expression data
  if (metadata(se)$technology == "rnaseq") {
    assay = "expr.voom"
    assay.voom.weights = "voom.weights"
  } else {
    if ("expr.corrected" %in% names(assays(se))) {
      assay = "expr.corrected"
    } else {
      assay = "expr"
    }
    assay.voom.weights = NULL
  }
  
  if (!(assay %in% names(assays(se)))) {
    stop(paste(assay, "not available in assays!"))
  }
  expr = assays(se)[[assay]]
  
  ## extract VOOM weights (only relevant for RNA-seq studies)
  if (!is.null(assay.voom.weights)) {
    if (!(assay.voom.weights %in% names(assays(se)))) {
      stop(paste(assay, "not available in assays!"))
    }
    weights = assays(se)[[assay.voom.weights]]
    expr = new("EList", 
               list(E = expr,
                    weights = weights))
    
  }
  
  ## extract phenotype data
  pheno = as.data.frame(colData(se)[, c(var, covar, var.id)])
  
  ## differential expression analysis
  fitmm = dream(exprObj = expr,
                formula = as.formula(form.final),
                pheno)
  
  ## extract results and adjust for multiple testing
  res = topTable(fitmm,
                 coef = 2,
                 number = nrow(expr),
                 adjust.method = "BH",
                 sort.by = "none")
  
  ## calculate SE
  res$SE = res$logFC / res$t
  
  ## add mean expression per group
  info.mean = estimate.means(gr = pheno[, var],
                             expr = expr)
  
  ## store results in data.frame
  res = data.frame(gene = rownames(se),
                   res,
                   info.mean,
                   rowData(se),
                   stringsAsFactors = FALSE)
  
  ## save in file
  export(res,
         file = res.file)
  
  #  return(res)
  
}

estimate.means <- function(gr, expr) {
  if (is.factor(gr) || is.character(gr) || length(unique(gr)) == 2) {
    if (is.matrix(expr)) {
      x = expr
    } else {
      x = expr$E
    }
    expr.l = split(data.frame(t(x)), gr)
    info.mean = sapply(expr.l, colMeans)
    colnames(info.mean) = paste0("mean.", colnames(info.mean))
    sum = apply(info.mean, 2, function(x) {sum(!(is.nan(x)))})
    info.mean = info.mean[, sum > 0]
  } else {
    info.mean = NULL
  }
  return(info.mean)
}

load.csv <- function(csv.files) {
  csv.l = lapply(csv.files, 
                import)
  names = sapply(csv.files,
                 function(f) {
                   temp = unlist(strsplit(
                     gsub(".csv", "", basename(f)), "_"))
                   return(temp[4])})
  names(names) = NULL
  names(csv.l) = names
  return(csv.l)
} 




volcano.plot <- function(data_ges, log.p = TRUE, log.fc = FALSE, main, genes.sel = NULL) {
  data = data_ges[,c("gene","P.Value","logFC")]
  colnames(data) <- c("id","p","fc")
  if (log.p) data$p = -log10(data$p)
  if (log.fc) data$fc = log2(data$fc)
  
 plot(data$fc, data$p, pch = 1, main = main,
       xlab = "logFC",
       ylab = "-log10(P-Wert) unadjustiert")
 if (!is.null(genes.sel)) {
    points(subset(data, id %in% genes.sel, c("fc", "p")),
           pch = 19, col = "red")
  }
}   

log.his.p <- function(data_ges, main){
          data_ges$m.log.p.val <- ifelse(data_ges$adj.P.Val != 0 , - log(data_ges$adj.P.Val) , NA)
          hist( data_ges$m.log.p.val, main = main,
                xlab = "-log(P-Wert)",
                ylab = "Häufigkeit")
}

his.plot.p <- function(data_ges, main){
            hist( data_ges$P.Value, main = main,
                  xlab = "unadjustierter P-Wert",
                  ylab = "Häufigkeit")
}

fallzahl.anteil.plot <- function( data_ges,main,anzahl=1){
            plot(x=sum(data_ges$logFC <= 0)/length(data_ges$gene),y=anzahl,
            xlab="Anteil logFC < 0",
            ylab="Anzahl der Studienteilnehmer")
}



info.studien <- function(se.list){
  tabelle.info <- data.frame(Jahr = 1:length(se.list),
                             Technologie = 1:length(se.list),
                             anz.gen =1:length(se.list),
                             anz.les = 1:length(se.list),
                             anz.n.les = 1:length(se.list),
                             anz.prob = 1:length(se.list))
  for (i in 1:length(se.list)) {
    tabelle.info[i,1] <- metadata(se.l[[i]])$technology
    tabelle.info[i,1] <- ifelse(tabelle.info[i,1]=="array.affy","Microarrays","RNAseq")
    if(is.null(metadata(se.l[[i]])$date)){
      tabelle.info[i,2] <- metadata(se.l[[i]])$submission_date
    }
    else{
    tabelle.info[i,2] <- metadata(se.l[[i]])$date[1]
    }
    tabelle.info[i,3] <- dim(se.l[[i]])[1]
    tabelle.info[i,4] <- sum(colData(se.l[[i]])$biomap_lesional == "lesional")
    tabelle.info[i,5] <- sum(colData(se.l[[i]])$biomap_lesional == "nonlesional")
    tabelle.info[i,6] <- max(tabelle.info[i,4],tabelle.info[i,5])
  }
  row.names(tabelle.info) <- names(se.list)
  #tabelle.info$Technologie <- unlist(lapply(names(se.list),function(x){ metadata(se.l$x)$technology}))
  #tabelle.info$Jahr <- unlist(lapply(names(se.list),function(x){ metadata(se.l$x)$date}))
  #tabelle.info$anz.les <- unlist(lapply(se.list,function(x){ sum(colData(se.l$x)$biomap_lesional == "lesional")}))
  t#abelle.info$anz.n.les <- unlist(lapply(se.list,function(x){ sum(colData(se.l$x)$biomap_lesional == "nonlesional")}))
  
 
  colnames(tabelle.info)<- c("Technologie","Datum","Anzahl Gene", "Anzahl läsional Proben","Anzahl nicht läsional Proben","anz.teil")
  return(tabelle.info)
}



#Tabelle nach diff.es analyse
#Anzahl in  ###Studie1 #Studie2
#Anzahl gene
#Anzahl p.val <0.5
#Anteil p.val < 0.5
#Anteil fc negativ
summary.tabelle <- function(data_liste,gene){
  tabelle.zus <- data.frame(Anzahl_gene = 1:length(data_liste),
                            p.val_kleiner_0.5 = 1:length(data_liste),
                            Anteil_p.val = 1:length(data_liste),
                            Anteil_fc = 1:length(data_liste))
  tabelle.zus$Anzahl_gene <- unlist(lapply(data_liste, FUN = function(x){
    anzahl_gene <- length(x$gene)}))
  tabelle.zus$p.val_kleiner_0.5 <-unlist(lapply(data_liste, FUN = function(x){
    anzahl_p.val. <- sum(x$adj.P.Val<= 0.05)}) )
  tabelle.zus$Anteil_p.val <- tabelle.zus$p.val_kleiner_0.5 / tabelle.zus$Anzahl_gene
  tabelle.zus$Anteil_fc <- unlist(lapply(data_liste, FUN = function(x){
    anzahl_p.val. <- sum(x$logFC <= 0)}) )/tabelle.zus$Anzahl_gene
  row.names(tabelle.zus) <- names(data_liste)
  colnames(tabelle.zus)<- c("Anzahl Gene","Anzahl P-Wert < 0.05", "Anteil P-Wert < 0.05","Anteil logFC < 0")
  return(tabelle.zus)
}


#Tabelle für Studien info der Studien
#funktioniert nicht, wegen ungleichmäßiger besetzung von metdata
#summary.tabelle <- function(data_liste){
#  tabelle.zus <- data.frame(Institut = 1:length(data_liste),
#                            #Datum = 1:length(data_liste),
#                            Technologie = 1:length(data_liste))
#  for (i in 1:length(data_liste)){
#  tabelle.zus[1,i] <- ifelse(is.na(metadata(data_liste[[i]]$contact_institute) ),metadata(data_liste[[i]]$contact_institute),NA) 
#  #tabelle.zus[2,i] <- ifelse(is.na(metadata(data_liste[[i]]$date)),metadata(data_liste[[i]]$date) , NA)
#  tabelle.zus[3,i] <- ifelse( is.na( metadata(data_liste[[i]]$technology)),metadata(data_liste[[i]]$technology),NA)
#  
#  }
#  knitr::kable(tabelle.zus)
#  }
