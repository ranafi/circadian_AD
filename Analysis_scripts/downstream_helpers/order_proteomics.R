
library(tidyverse)
library(doParallel)
library(progress)

plot_prot = function( tmt_filename, path_to_cyclops_ordering,genelist,path_to_rosmap_clin, percentile = 0.025){
  
  setwd("../synapse_downloads")
  # path_to_rosmap_clin = "../ROSMAP_metadata/cleaned_rosmap_meta_cogdxConds.csv"
  prot = read_csv(tmt_filename, show_col_types = F)
  
  colnames(prot)[1] = "Protein_Symbols"
  data_key1 = read_csv("rosmap_50batch_specimen_metadata_for_batch_correction.csv", show_col_types = F) #metadata for batch 1 syn32835854
  data_key1$exp = 1
  data_key2 = readxl::read_xlsx("ROSMAP_Round2_Traits_FINAL.xlsx") #metadata for batch 2 syn21266449
  data_key2 = as.data.frame(data_key2)
  data_key2$exp = 2
  data_key2$Batch =  paste0("b",data_key2$Plex %>% as.numeric %>% + 50) #create a "batch" variable akin to data_key1
  data_key2$SampleID = paste0(data_key2$Batch,".", data_key2$Channel) %>% stringr::str_replace_all("_", "")
  colnames(data_key2)[5] = "projid"
  data_key = rbind(dplyr::select(data_key1, SampleID, projid, exp), dplyr::select(data_key2, SampleID, projid, exp) ) #combine the metadatas
  
  #IMPORTANT: some subjects were done in both batches,
  #right now I keep only the first batch versions, change fromLast = T otherwise
  data_key = data_key[!duplicated(data_key$projid, fromLast = F), ]
  
  ros_clin = read_csv(path_to_rosmap_clin, show_col_types = F)
  data_key = merge(data_key, ros_clin, by = "projid")
  keep_cols = which(colnames(prot) %in% data_key$SampleID) #columns of prot for which I have metadata
  prot = prot[, c(1, keep_cols)] #keep cols and col 1 (Protein_Symbols)
  colnames(prot)[-1] = data_key$projid[match(colnames(prot)[-1], data_key$SampleID)] #rename colms of prot to rosmap projid
  
  cyc_pred_file = list.files(path = paste0(path_to_cyclops_ordering, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(path_to_cyclops_ordering, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = F)
  cyc_pred = cyc_pred[na.exclude(match(colnames(prot), cyc_pred$ID)),] #keep the samples in cyc_pred that I have protein data for
  cyc_pred$pmi = ros_clin$pmi[match(cyc_pred$ID, ros_clin$projid)]
  cyc_pred$sex = ros_clin$msex[match(cyc_pred$ID, ros_clin$projid)]
  keep_subs = c(1, which(colnames(prot) %in% cyc_pred$ID))
  prot= prot[, keep_subs]  #keep the protein samples I have ordering for
  all(colnames(prot)[-1] == cyc_pred$ID) #double check all protein emat columns equal cyc_pred columns
  
  tmm = as.data.frame(prot)
  tmm = rbind(c("Cond_D", cyc_pred$Covariate_D),c("pmi_C", cyc_pred$pmi),
              c("sex_D", paste0("cond_", cyc_pred$sex)), prot ) 
  
  colnames(tmm)[1] = "Gene_Symbols"
  plot_gene_trace(cyc_pred = cyc_pred, tmm = tmm, seedlist = genelist, savePlots = F, adjust_points_for_batches = F)
}


order_prot = function( tmt_filename, path_to_cyclops_ordering ,path_to_rosmap_clin, isCyclingBHQCutoff = 0.1, percentile = 0.025, plot_prot = T){
  
  setwd("../synapse_downloads")
  # path_to_rosmap_clin = "../ROSMAP_metadata/cleaned_rosmap_meta_cogdxConds.csv"
  prot = read_csv(tmt_filename, show_col_types = F)
  
  colnames(prot)[1] = "Protein_Symbols"
  data_key1 = read_csv("rosmap_50batch_specimen_metadata_for_batch_correction.csv", show_col_types = F) #metadata for batch 1 syn32835854
  data_key1$exp = 1
  data_key2 = readxl::read_xlsx("ROSMAP_Round2_Traits_FINAL.xlsx") #metadata for batch 2 syn21266449
  data_key2 = as.data.frame(data_key2)
  data_key2$exp = 2
  data_key2$Batch =  paste0("b",data_key2$Plex %>% as.numeric %>% + 50) #create a "batch" variable akin to data_key1
  data_key2$SampleID = paste0(data_key2$Batch,".", data_key2$Channel) %>% stringr::str_replace_all("_", "")
  colnames(data_key2)[5] = "projid"
  data_key = rbind(dplyr::select(data_key1, SampleID, projid, exp), dplyr::select(data_key2, SampleID, projid, exp) ) #combine the metadatas
  
  #IMPORTANT: some subjects were done in both batches,
  #right now I keep only the first batch versions, change fromLast = T otherwise
  data_key = data_key[!duplicated(data_key$projid, fromLast = F), ]
  
  ros_clin = read_csv(path_to_rosmap_clin, show_col_types = F)
  data_key = merge(data_key, ros_clin, by = "projid")
  keep_cols = which(colnames(prot) %in% data_key$SampleID) #columns of prot for which I have metadata
  prot = prot[, c(1, keep_cols)] #keep cols and col 1 (Protein_Symbols)
  colnames(prot)[-1] = data_key$projid[match(colnames(prot)[-1], data_key$SampleID)] #rename colms of prot to rosmap projid
  
  cyc_pred_file = list.files(path = paste0(path_to_cyclops_ordering, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(path_to_cyclops_ordering, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = F)
  cyc_pred = cyc_pred[na.exclude(match(colnames(prot), cyc_pred$ID)),] #keep the samples in cyc_pred that I have protein data for
  cyc_pred$pmi = ros_clin$pmi[match(cyc_pred$ID, ros_clin$projid)]
  cyc_pred$sex = ros_clin$msex[match(cyc_pred$ID, ros_clin$projid)]
  keep_subs = c(1, which(colnames(prot) %in% cyc_pred$ID))
  prot= prot[, keep_subs]  #keep the protein samples I have ordering for
  all(colnames(prot)[-1] == cyc_pred$ID) #double check all protein emat columns equal cyc_pred columns
  
  tmm = as.data.frame(prot)
  tmm = rbind(c("Cond_D", cyc_pred$Covariate_D),c("pmi_C", cyc_pred$pmi),
              c("sex_D", paste0("cond_", cyc_pred$sex)), prot ) #Add condition to protein expression
  
  ###pca##
  # complete_prot = prot[complete.cases(prot), ] %>% column_to_rownames(var = "Protein_Symbols")
  # complete_prot = apply(complete_prot, 1, blunt_outliers)
  # exp = data_key$exp[match(rownames(complete_prot), data_key$projid)]
  # pca_res = prcomp(complete_prot, scale. = T)
  # ggplot(as.data.frame(pca_res$x), aes(PC1, PC2, color = exp))+geom_point()
  
  #Is_cycling 
  pb = progress_bar$new(total = dim(prot)[1])
  cycling_in_CTL = is_cycling(cyc_pred, tmm, "cond_0", pb = pb, percentile = percentile)
  pb = progress_bar$new(total = dim(prot)[1])
  cycling_in_AD = is_cycling(cyc_pred, tmm, "cond_1", pb = pb, percentile = percentile)
  pb = progress_bar$new(total = dim(prot)[1])
  cycling_AD_CTL_method2 = is_cycling_method2(cyc_pred, tmm, pb = pb, percentile = percentile)
  
  #diff mesor for all proteins
  pb = progress_bar$new(total = dim(prot)[1])
  gene_list_mesor =  unlist(unname(tmm[!grepl("_D|_C", unlist(tmm[,1])), 1])) # TEST ALL genes for Mesor diff (not just cyclers)
  diff_mesor = mesor_differences(cyc_pred, tmm, gene_list_mesor, pb = pb, percentile = percentile)
  
  #Create cycling lists for diff_rhythms
  strong_cycling_in_CTL = dplyr::filter(cycling_in_CTL, as.numeric(BHQ) < isCyclingBHQCutoff) %>% arrange(as.numeric(Bonf))
  strong_cycling_in_AD = dplyr::filter(cycling_in_AD, as.numeric(BHQ) < isCyclingBHQCutoff) %>% arrange(as.numeric(Bonf))
  seedlist = union(strong_cycling_in_CTL$Gene_Symbols, strong_cycling_in_AD$Gene_Symbols)
  strong_cyclers_method2 = dplyr::filter(cycling_AD_CTL_method2, as.numeric(BHQ) < isCyclingBHQCutoff) %>% arrange(as.numeric(BHQ))
  
  diff_rhythms = diff_rhyth(cyc_pred, tmm, seedlist, percentile = percentile)
  diff_rhythms_mthd2 = diff_rhyth(cyc_pred, tmm, strong_cyclers_method2$Gene_Symbols, percentile = percentile)
  
  setwd(paste0(path_to_cyclops_ordering, "/proteomics"))
  # plotting 
  if (plot_prot){
    seedlist = c(strong_cycling_in_CTL$Gene_Symbols[1:3], strong_cycling_in_AD$Gene_Symbols[1:3])
    colnames(tmm)[1] = "Gene_Symbols"
    print(plot_gene_trace(cyc_pred = cyc_pred, tmm = tmm, seedlist = seedlist, savePlots = T, adjust_points_for_batches = F))
  }
  isCyclingSigCutoff_str = str_extract(as.character(isCyclingBHQCutoff), "(?<=\\.)\\d+")
  percentile_str = str_extract(as.character(percentile), "(?<=\\.)\\d+")
  # write out files 
  diff_rhythms_mthd2 = separate(diff_rhythms_mthd2, Gene_Symbols, sep = "\\|", into = c("Symbol","Uniprot"))
  write.table(diff_rhythms_mthd2, paste0( "diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,"_bluntingPercentile", percentile_str,".csv"), col.names = T, row.names = F, sep = ',')
  
  diff_rhythms_mthd1 = separate(diff_rhythms, Gene_Symbols, sep = "\\|", into = c("Symbol","Uniprot"))
  write.table(diff_rhythms_mthd1, paste0("diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"_bluntingPercentile", percentile_str,".csv"), col.names = T, row.names = F, sep = ',')
  
  cycling_in_CTL = separate(cycling_in_CTL, Gene_Symbols, sep = "\\|", into = c("Symbol","Uniprot"))
  write.table(cycling_in_CTL, paste0("cycling_in_CTL_CyclingBHQ",isCyclingSigCutoff_str,"_bluntingPercentile", percentile_str,".csv"), col.names = T, row.names = F, sep = ',')
  
  cycling_in_AD = separate(cycling_in_AD, Gene_Symbols, sep = "\\|", into = c("Symbol","Uniprot"))
  write.table(cycling_in_AD, paste0("cycling_in_AD_CyclingBHQ",isCyclingSigCutoff_str,"_bluntingPercentile", percentile_str,".csv"), col.names = T, row.names = F, sep = ',')

  diff_mesor = separate(diff_mesor, Gene_Symbols, sep = "\\|", into = c("Symbol","Uniprot"))
  write.table(diff_mesor, paste0("diff_mesor_CyclingBHQ",isCyclingSigCutoff_str,"_bluntingPercentile", percentile_str,".csv"), col.names = T, row.names = F, sep = ',')
  
}

write_prot_rnks = function(path_to_cyclops_ordering){
  print("Creating rnk files for fGSEA.")
  ## is cycling in CTL ranked by -log(p_val)
  CTL_cyclers_file = list.files(path = paste0(path_to_cyclops_ordering, "/proteomics"), pattern = "cycling_in_CTL_CyclingBHQ\\d+_blunting.*csv")
  if (!purrr::is_empty(CTL_cyclers_file)){
    CTL_cyclers = read_csv(paste0(path_to_cyclops_ordering, "/proteomics/", CTL_cyclers_file), show_col_types = FALSE)
    ranked_CTL_cyclers = arrange(CTL_cyclers, p_statistic)
    df1 = data.frame(genes = ranked_CTL_cyclers$Uniprot, metric = -log(ranked_CTL_cyclers$p_statistic))
    write.table(df1, paste0(path_to_cyclops_ordering, "/proteomics/fGSEA/rnk_files/CTL_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  } 
  ## is cycling in AD ranked by -log(p_val)
  AD_cyclers_file = list.files(path = paste0(path_to_cyclops_ordering, "/proteomics"), pattern = "cycling_in_AD_CyclingBHQ\\d+_blunting.*csv")
  if (!purrr::is_empty(AD_cyclers_file)){
    AD_cyclers = read_csv(paste0(path_to_cyclops_ordering, "/proteomics/", AD_cyclers_file), show_col_types = F)
    ranked_AD_cyclers = arrange(AD_cyclers, p_statistic)
    df2 = data.frame(genes = ranked_AD_cyclers$Uniprot, metric = -log(ranked_AD_cyclers$p_statistic))
    write.table(df2, paste0(path_to_cyclops_ordering, "/proteomics/fGSEA/rnk_files/AD_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  }
  ## DR cyclers ranked by -log(p_val)
  DR_cyclers_file = list.files(path = paste0(path_to_cyclops_ordering, "/proteomics"), pattern = "diff_rhythms_Cycling.*.csv")
  if (!purrr::is_empty(DR_cyclers_file)){
    DR_cyclers = read_csv(paste0(path_to_cyclops_ordering, "/proteomics/", DR_cyclers_file), show_col_types = F)
    ranked_DR_cyclers = arrange(DR_cyclers, p_val)
    df3 = data.frame(genes = ranked_DR_cyclers$Uniprot, metric = -log(ranked_DR_cyclers$p_val))
    write.table(df3, paste0(path_to_cyclops_ordering, "/proteomics/fGSEA/rnk_files/DR_cyclers_minusLogPRanked.rnk"), sep = '\t', col.names = F, row.names = F)
  }
  
  ## DR cyclers ranked by log(ADamp/CTLamp)
  DR_cyclers_file = list.files(path = paste0(path_to_cyclops_ordering, "/proteomics"), pattern = "diff_rhythms_Cycling.*\\.csv")
  if (!purrr::is_empty(DR_cyclers_file)){
    DR_cyclers = read_csv(paste0(path_to_cyclops_ordering, "/proteomics/", DR_cyclers_file), show_col_types = F)
    ranked_DR_cyclers = arrange(DR_cyclers, Log_AD_CTL_ampRatio)
    df4 = data.frame(genes = ranked_DR_cyclers$Uniprot, metric = ranked_DR_cyclers$Log_AD_CTL_ampRatio)
    write.table(df4, paste0(path_to_cyclops_ordering, "/proteomics/fGSEA/rnk_files/DR_cyclers_Log(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  }
  
  ## Diff Mesor ranked by log(ADmesor/CTLmesor)
  diff_mesor_file = list.files(path = paste0(path_to_cyclops_ordering, "/proteomics"), pattern = "diff_mesor.*\\.csv")
  if (!purrr::is_empty(diff_mesor_file)){
    diff_mesor = read_csv(paste0(path_to_cyclops_ordering, "/proteomics/", diff_mesor_file), show_col_types = F)
    diff_mesor$AD_minus_CTL_mes = as.numeric(diff_mesor$mesor_AD) - as.numeric(diff_mesor$mesor_CTL)
    ranked_diff_mesor = arrange(diff_mesor, AD_minus_CTL_mes)
    df5 = data.frame(genes = ranked_diff_mesor$Uniprot, metric = ranked_diff_mesor$AD_minus_CTL_mes)
    write.table(df5, paste0(path_to_cyclops_ordering, "/proteomics/fGSEA/rnk_files/diff_mesor_(AD-CTL)ranked.rnk"), sep = '\t', col.names = F, row.names = F)
  }
}
