library(getmstatistic)
library(base)
library(caret)
library(corrplot)
library(getspres) 
library(ggpubr)
library(forcats)

#siginfikanten Gene aus Studien
#eingabe der Ergebnisse der diff. ex. analyse, meta-Analyse

sig_raus_stud <- function(diff.ex.l=diff.ex.l,meta=meta, rel_gen){
  
  
  
  rel_in_studie <- list(NA)
  for(i in 1:length(diff.ex.l)){
    rel_in_studie[[i]] <- diff.ex.l[[i]][diff.ex.l[[i]]$gene%in%rel_gen,]
    #rel_in_studie <- cbind(rel_in_studie, diff.ex.l[[i]]$gene[diff.ex.l[[i]]$gene%in%rel_gen])
  }
  #rel_in_studie <- rel_in_studie[,-1]
  #colnames(rel_in_studie) <- names(diff.ex.l)
  #names(rel_in_studie) <- names(diff.ex.l)
  ges_rel_studien <- NA
  for(i in 1:length(diff.ex.l)){
    studie.name = rep(names(diff.ex.l)[i],times=length(rel_in_studie[[i]][1]))
    studie <- cbind(rel_in_studie[[i]]$gene,
                    rel_in_studie[[i]]$logFC,
                    rel_in_studie[[i]]$SE,
                    rel_in_studie[[i]]$adj.P.Val,
                    studie.name)
    ges_rel_studien <- rbind(ges_rel_studien,studie)
  }
  ges_rel_studien<-ges_rel_studien[-1,]
  colnames(ges_rel_studien)<-c("gene","logFC","SE","adj.P.Val","studie.name")
  return(ges_rel_studien)
  
}



korrelation <- function(korrelation, cutoff=0.8, meta , diff.ex.l,anzahl_gene, file){

  #Gene mit Cut of bei Korrelation 0.8
  cut<-findCorrelation( korrelation, cutoff , exact = FALSE,names = TRUE)
  
  
  gene <- setdiff(meta$gene, cut) 

  

  #Gene aus allen Studien mit |logFC|>1 und Korrelation <0.8
  sig_gene <- as.data.frame(sig_raus_stud(diff.ex.l,meta,gene))

  #Berechnen der M-Statistik mit den korrelierten Genen
  mstatistic <- getmstatistic(beta_in = as.numeric(sig_gene$logFC),
                               lambda_se_in = as.numeric(sig_gene$SE),
                               study_names_in = sig_gene$gene,
                               variant_names_in = sig_gene$studie.name
                                )

  #Exportieren der unkorrelierten M-Statistik
  m_data <- mstatistic$M_dataset
  m_data$korrelation <- rep(cutoff, times= length(m_data$study_names_in))
  export(m_data,file = file)
  

}




