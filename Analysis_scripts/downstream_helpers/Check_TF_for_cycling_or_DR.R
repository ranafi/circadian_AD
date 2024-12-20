#####
# Tues Sept 5 2023
# Script for taking a df containing TF_names (from pscan but also enrichR results),
# and checking if they are found in DE, DE, DM, cycling files. 
library(ggrepel)
library(tidyverse)
augment_tf_file = function(TF_filename, Diff_expr_filename, isCyclingBHQCutoff_str, abs_path_to_cyclops_ordering, use_PRTNS = T){
  current_dir = getwd()
  TF_file = read.csv(TF_filename)
  if (colnames(TF_file)[1]== "Rank"){ #if file comes from enrichR, make it look like pscan
    TF_file$TF_NAME = TF_file$Term.name
    TF_file$FDR = TF_file$Adjusted.p.value
    TF_file = TF_file[,-1]
    TF_file$TF_NAME = str_replace(TF_file$TF_NAME , " \\(human\\)", "")
    TF_file$TF_NAME = str_replace(TF_file$TF_NAME , " \\(mouse\\)", "")
    
  }
  edgeR_toptags = read.csv(paste(abs_path_to_cyclops_ordering, Diff_expr_filename, sep = '/'))
  DR_AR1_mthd2 = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"))
  DR_AR1 = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"))
  cycling_CTL = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/cosinor_results_CTL.csv"))
  cycling_AD = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/cosinor_results_AD.csv"))
  mesor_file = list.files(path = paste0(abs_path_to_cyclops_ordering, "/downstream_output") ,pattern = "\\.*mesor_all_genes\\.csv$")
  Diff_mesor = read.csv(paste0(abs_path_to_cyclops_ordering, "/downstream_output/", mesor_file))
  #Try to read in Protein cycling/DR results:
  use_PRTN_data = F
  if(dir.exists(paste0(abs_path_to_cyclops_ordering, "/proteomics")) & use_PRTNS){
    use_PRTN_data = T
    
    cycling_PRTN_ctl_filename = list.files(path = paste0(abs_path_to_cyclops_ordering, "/proteomics/"),
                                           pattern = "cycling_in_CTL_.*csv")
    cycling_PRTN_ctl = read.csv(paste0(abs_path_to_cyclops_ordering, "/proteomics/",cycling_PRTN_ctl_filename))
    cycling_PRTN_ad_filename = list.files(path = paste0(abs_path_to_cyclops_ordering, "/proteomics/"),
                                           pattern = "cycling_in_AD_.*csv")
    cycling_PRTN_ad = read.csv(paste0(abs_path_to_cyclops_ordering, "/proteomics/",cycling_PRTN_ad_filename))
    DM_PRTN_filename = list.files(path = paste0(abs_path_to_cyclops_ordering, "/proteomics/"),
                                           pattern = "diff_mesor_.*csv")
    DM_PRTN = read.csv(paste0(abs_path_to_cyclops_ordering, "/proteomics/",DM_PRTN_filename))
    DR_PRTN_filename = list.files(path = paste0(abs_path_to_cyclops_ordering, "/proteomics/"),
                                  pattern = "diff_rhythms_Cycling.*csv")
    DR_PRTN = read.csv(paste0(abs_path_to_cyclops_ordering, "/proteomics/", DR_PRTN_filename))
  }
    
  #Look for the TF name in the files above
  TF_file$cycling_in_CTL_BHQ = cycling_CTL$BHQ[match(toupper(TF_file$TF_NAME), cycling_CTL$Gene_Symbols)]
  TF_file$cycling_in_AD_BHQ = cycling_AD$BHQ[match(toupper(TF_file$TF_NAME), cycling_AD$Gene_Symbols)]
  TF_file$DR_AR1_BHQ = DR_AR1$BHQ[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)]
  TF_file$DR_AR1_mthd2_BHQ = DR_AR1_mthd2$BHQ[match(toupper(TF_file$TF_NAME), DR_AR1_mthd2$Gene_Symbols)]
  TF_file$DR_logAmpRatio = DR_AR1$Log_AD_CTL_ampRatio[match(toupper(TF_file$TF_NAME), DR_AR1$Gene_Symbols)] 
  TF_file$diff_mesor_BHQ = Diff_mesor$BHQ[match(toupper(TF_file$TF_NAME), Diff_mesor$Gene_Symbols)]
  TF_file$edgeR_DE_BHQ = edgeR_toptags$FDR[match(toupper(TF_file$TF_NAME), edgeR_toptags$X)] 
  
  #Do the same for the PROTEIN FILES:
  if(use_PRTN_data){
    TF_file$cycling_PRTN_CTL_BHQ = cycling_PRTN_ctl$BHQ[match(toupper(TF_file$TF_NAME), cycling_PRTN_ctl$Symbol)]
    TF_file$cycling_PRTN_AD_BHQ = cycling_PRTN_ad$BHQ[match(toupper(TF_file$TF_NAME), cycling_PRTN_ad$Symbol)]
    TF_file$DR_PRTN_BHQ = DR_PRTN$BHQ[match(toupper(TF_file$TF_NAME), DR_PRTN$Symbol)]
    TF_file$DM_PRTN_BHQ = DM_PRTN$BHQ[match(toupper(TF_file$TF_NAME), DM_PRTN$Symbol)]
  }
  #write out file TF_file
  write.table(TF_file, TF_filename, row.names = F, col.names = T, sep = ',')
  
  #Create dataframes for plotting:
  DR_tfs = dplyr::filter(TF_file, FDR < 0.1 & DR_AR1_BHQ < 0.3) 
  cycling_CTL_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_in_CTL_BHQ < 0.1) 
  cycling_AD_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_in_AD_BHQ < 0.1)
  DE_tfs = dplyr::filter(TF_file, FDR < 0.1 & edgeR_DE_BHQ < 0.1) 
  DM_tfs = dplyr::filter(TF_file, FDR < 0.1 & diff_mesor_BHQ < 0.1) 
  
  if(use_PRTN_data){
    cycling_PRTN_CTL_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_PRTN_CTL_BHQ < 0.2)
    cycling_PRTN_AD_tfs = dplyr::filter(TF_file, FDR < 0.1 & cycling_PRTN_AD_BHQ < 0.2) 
    DR_PRTN_tfs = dplyr::filter(TF_file, FDR < 0.1 & DR_PRTN_BHQ < 0.3)
    DM_PRTN_tfs = dplyr::filter(TF_file, FDR < 0.1 & DM_PRTN_BHQ < 0.3)
  }
  
  if(!(dir.exists(paste0(tools::file_path_sans_ext(TF_filename), "_plots")))){
    dir.create(paste0(tools::file_path_sans_ext(TF_filename), "_plots"))
  }
  setwd(paste0(tools::file_path_sans_ext(TF_filename), "_plots"))
  #Check diff rhythms for TFs:
  if (nrow(DR_tfs) > 0){
      p1 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(DR_AR1_BHQ)))+
        geom_label_repel(data = DR_tfs,
                         aes(x = -log10(FDR), y = -log10(DR_AR1_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.35, point.padding = 0.5)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(DR BH.q value) (cycling AR>0.1, BH.q<0.3)")+
        ggtitle("Differentially Rhythmic TF enriched in gene list")+
        theme_minimal()
      ggsave("DR_TF_enriched_in_list.png", p1, width = 6, height = 5, units = "in")
  }
  #Check genes cycling in CTL for TF:
  if (nrow(cycling_CTL_tfs)>0){
      p2 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_in_CTL_BHQ)))+
        geom_label_repel(data = cycling_CTL_tfs,max.overlaps = 20,
                         aes(x = -log10(FDR), y = -log10(cycling_in_CTL_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.05)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(Cycling in CTL BH.q)")+
        ggtitle("Cycling in CTL TF enriched in gene list")+
        theme_minimal()
      ggsave("cycling_CTL_TF_enriched_in_list.png", p2, width = 6, height = 5, units = "in")
  }
  #Check genes cycling in AD for TF:
  if (nrow(cycling_AD_tfs)>0){
      p3 =  ggplot(TF_file)+
          geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_in_AD_BHQ)))+
          geom_label_repel(data = cycling_AD_tfs,
                           aes(x = -log10(FDR), y = -log10(cycling_in_AD_BHQ),label = TF_NAME),
                           nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                           box.padding = 0.05)+
          xlab("-Log(Pscan/enrichR FDR)")+
          ylab("-Log(Cycling in AD BH.q)")+
          ggtitle("Cycling in AD TF enriched in gene list")+
          theme_minimal()
      ggsave("cycling_AD_TF_enriched_in_list.png", p3, width = 6, height = 5, units = "in")
  }
  #Check Diff Expressed genes for TF:
  if(nrow(DE_tfs)>0){
    p4 =  ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(edgeR_DE_BHQ)))+
        geom_label_repel(data = DE_tfs,
                         aes(x = -log10(FDR), y = -log10(edgeR_DE_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue", 
                         box.padding = 0.35)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(edgeR Diff. Expr. BH.q)")+
        ggtitle("EdgeR diff expr TF enriched in gene list")+
        theme_minimal()
    ggsave("DE_TF_enriched_in_list.png", p4, width = 6, height = 5, units = "in")
  }
  #Check genes with different mesor for TF:
  if(nrow(DM_tfs)>0){
    p5 =  ggplot(TF_file)+
      geom_point(mapping = aes(x = -log10(FDR), y = -log10(diff_mesor_BHQ)))+
      geom_label_repel(data = DM_tfs,
                       aes(x = -log10(FDR), y = -log10(diff_mesor_BHQ),label = TF_NAME),
                       nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                       box.padding = 0.35)+
      xlab("-Log(Pscan/enrichR FDR)")+
      ylab("-Log(Diff Mesor BH.q)")+
      ggtitle("Diff Mesor TF enriched in gene list")+
      theme_minimal()
    ggsave("DM_TF_enriched_in_list.png", p5, width = 6, height = 5, units = "in")
  }
  if(use_PRTN_data){
    #PROTEIN plots:
    if (nrow(DR_PRTN_tfs) > 0){
      p6 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(DR_PRTN_BHQ)))+
        geom_label_repel(data = DR_PRTN_tfs,
                         aes(x = -log10(FDR), y = -log10(DR_PRTN_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.35)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(DR BH.q value) (cycling AR>0.1, BH.q<0.3)")+
        ggtitle("Differentially Rhythmic TF PROTEINS enriched in gene list")+
        theme_minimal()
      ggsave("DR_PRTN_TF_enriched_in_list.png", p6, width = 6, height = 5, units = "in")
    }
    if (nrow(cycling_PRTN_CTL_tfs)>0){
      p7 = ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_PRTN_CTL_BHQ)))+
        geom_label_repel(data = cycling_PRTN_CTL_tfs, max.overlaps = 20,
                         aes(x = -log10(FDR), y = -log10(cycling_PRTN_CTL_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.05)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(Cycling in CTL BH.q)")+
        ggtitle("Cycling in CTL TF PROTEIN enriched in gene list")+
        theme_minimal()
      ggsave("cycling_PRTN_CTL_TF_enriched_in_list.png", p7, width = 6, height = 5, units = "in")
    }
    if (nrow(cycling_PRTN_AD_tfs)>0){
      p8 =  ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(cycling_PRTN_AD_BHQ)))+
        geom_label_repel(data = cycling_PRTN_AD_tfs,
                         aes(x = -log10(FDR), y = -log10(cycling_PRTN_AD_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.05)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(Cycling in AD BH.q)")+
        ggtitle("Cycling in AD TF PROTEIN enriched in gene list")+
        theme_minimal()
      ggsave("cycling_PRTN_AD_TF_enriched_in_list.png", p8, width = 6, height = 5, units = "in")
    }
    if(nrow(DM_PRTN_tfs)>0){
      p9 =  ggplot(TF_file)+
        geom_point(mapping = aes(x = -log10(FDR), y = -log10(DM_PRTN_BHQ)))+
        geom_label_repel(data = DM_PRTN_tfs,
                         aes(x = -log10(FDR), y = -log10(DM_PRTN_BHQ),label = TF_NAME),
                         nudge_x = 0.5, nudge_y = 0.5,color = "blue",
                         box.padding = 0.35)+
        xlab("-Log(Pscan/enrichR FDR)")+
        ylab("-Log(Diff Mesor BH.q)")+
        ggtitle("Diff Mesor TF PROTEIN enriched in gene list")+
        theme_minimal()
      ggsave("DM_PRTN_TF_enriched_in_list.png", p9, width = 6, height = 5, units = "in")
    }
  }
  setwd(current_dir)
}
