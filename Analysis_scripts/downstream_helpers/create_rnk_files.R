#File to write out .rnk files uses in GSEA, updated Nov 2, 2023 to 
#work with my differential_rhyth_generalized.R which replaces use of compareRhythms
library(tidyverse)


write_rnks = function(path_to_cyclops_ordering, isCyclingBHQCutoff_str){
  print("Creating rnk files for fGSEA.")
  ## is cycling in CTL
  CTL_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_CTL.csv"), show_col_types = FALSE)
  ranked_CTL_cyclers = arrange(CTL_cyclers, p_statistic)
  df1 = data.frame(genes = ranked_CTL_cyclers$Gene_Symbols, metric = -log(ranked_CTL_cyclers$p_statistic))
  write.table(df1, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/CTL_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  ## is cycling in AD
  AD_cyclers = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/cosinor_results_AD.csv"), show_col_types = FALSE)
  ranked_AD_cyclers = arrange(AD_cyclers, p_statistic)
  df2 = data.frame(genes = ranked_AD_cyclers$Gene_Symbols, metric = -log(ranked_AD_cyclers$p_statistic))
  write.table(df2, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### DR genes #####
  # Cycling with AR > .1 & BHQ < cutoff Cycling -> DR testing, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, Log_AD_CTL_ampRatio)
  df3 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = ranked_DR_genes_AR1$Log_AD_CTL_ampRatio)
  write.table(df3, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .1 & BHQ < cutoff Cycling -> DR testing, -log(p_val) ranked
  ranked_DR_genes_AR1 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio1.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR1 = arrange(ranked_DR_genes_AR1, p_val)
  df4 = data.frame(genes = ranked_DR_genes_AR1$Gene_Symbols, metric = -log(ranked_DR_genes_AR1$p_val))
  write.table(df4, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio1_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .20 & BHQ < cutoff Cycling -> DR testing, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR20 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio20.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR20 = arrange(ranked_DR_genes_AR20, Log_AD_CTL_ampRatio)
  df5 = data.frame(genes = ranked_DR_genes_AR20$Gene_Symbols, metric = ranked_DR_genes_AR20$Log_AD_CTL_ampRatio)
  write.table(df5, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio20_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling with AR > .20 & BHQ < cutoff Cycling -> DR testing, -log(p_val) ranked
  ranked_DR_genes_AR20 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio20.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR20 = arrange(ranked_DR_genes_AR20, p_val)
  df6 = data.frame(genes = ranked_DR_genes_AR20$Gene_Symbols, metric = -log(ranked_DR_genes_AR20$p_val))
  write.table(df6, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio20_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling via Method 2 (BHQ < cutoff & AR > 0.20) -> diff rhythms, -log(p_val value DR) RANKED
  ranked_DR_genes_AR20_method2 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio20.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR20_method2 = arrange(ranked_DR_genes_AR20_method2, p_val)
  df7 = data.frame(genes = ranked_DR_genes_AR20_method2$Gene_Symbols, metric = -log(ranked_DR_genes_AR20_method2$p_val))
  write.table(df7, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio20_minusLogPRanked_method2.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #Cycling via Method 2 (BHQ < cutoff & AR > 0.20) -> diff rhythms, log(AD_amp/CTL_amp) ranked
  ranked_DR_genes_AR20_method2 = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingBHQCutoff_str,"AmpRatio20.csv"), show_col_types = FALSE)
  ranked_DR_genes_AR20_method2 = arrange(ranked_DR_genes_AR20_method2, Log_AD_CTL_ampRatio)
  df8 = data.table(genes = ranked_DR_genes_AR20_method2$Gene_Symbols, metric = ranked_DR_genes_AR20_method2$Log_AD_CTL_ampRatio)
  write.table(df8, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/DRgenesAmpRatio20_Log(AD-CTL)ranked_method2.rnk"), sep = '\t', col.names = F, row.names = F)
  
  #### differential mesor p ranked
  ranked_DM_genes = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/differential_mesor_all_genes.csv"), show_col_types = FALSE)
  ranked_DM_genes = arrange(ranked_DM_genes, p_mesor)
  df9 = data.table(genes = ranked_DM_genes$Gene_Symbols, metric = -log(ranked_DM_genes$p_mesor))
  write.table(df9, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/differential_mesor_all_genes_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)

  ###### AD severity DR p ranked #######
  ranked_cerad_genes = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AD_severity_AR20.csv"), show_col_types = FALSE)
  ranked_cerad_genes = arrange(ranked_cerad_genes, p_ceradsc)
  df_cerad = data.table(genes = ranked_cerad_genes$Gene_Symbols, metric = -log(ranked_cerad_genes$p_ceradsc))
  write.table(df_cerad, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_severity_ceradsc_AR20_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F )
  
  ranked_cogdx_genes = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AD_severity_AR20.csv"), show_col_types = FALSE)
  ranked_cogdx_genes = arrange(ranked_cogdx_genes, p_cogdx)
  df_cogdx = data.table(genes = ranked_cogdx_genes$Gene_Symbols, metric = -log(ranked_cogdx_genes$p_cogdx))
  write.table(df_cogdx, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_severity_cogdx_AR20_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F )
  
  ranked_braak_genes = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/diff_rhythms_AD_severity_AR20.csv"), show_col_types = FALSE)
  ranked_braak_genes = arrange(ranked_braak_genes, p_braak)
  df_braak = data.table(genes = ranked_braak_genes$Gene_Symbols, metric = -log(ranked_braak_genes$p_braak))
  write.table(df_braak, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_severity_braak_AR20_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F )
  
  ##### AD severity DM p ranked #######
  ranked_cerad_mesor = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/differential_mesor_all_genes_cogdx_cerad_braak.csv"), show_col_types = FALSE)
  ranked_cerad_mesor_genes = arrange(ranked_cerad_mesor, p_mesor_cerad)
  df_cerad_mesor = data.table(genes = ranked_cerad_mesor_genes$Gene_Symbols, metric = -log(ranked_cerad_mesor_genes$p_mesor_cerad))
  write.table(df_cerad_mesor, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_severity_mesor_p_ceradsc_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F )
  
  ranked_braak_mesor = read_csv(paste0(path_to_cyclops_ordering, "downstream_output/differential_mesor_all_genes_cogdx_cerad_braak.csv"), show_col_types = FALSE)
  ranked_braak_mesor_genes = arrange(ranked_braak_mesor, p_mesor_braak)
  df_braak_mesor = data.table(genes = ranked_braak_mesor_genes$Gene_Symbols, metric = -log(ranked_braak_mesor_genes$p_mesor_braak))
  write.table(df_braak_mesor, paste0(path_to_cyclops_ordering, "downstream_output/fGSEA/rnk_files/AD_severity_mesor_p_braak_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F )
  
  }
