#functions for using cyclops order on metabolon data
library(tidyverse)
library(doParallel)

order_metabolon = function(metabolon_filename, metabolon_datakey, path_to_rosmap_clin, path_to_cyclops_ordering, BHQ_cutoff = 0.1, plot_mets = T, percentile = 0.){
  metabolon = read_csv(metabolon_filename, show_col_types = FALSE) #read in metabolon assay
  metabolon_data_key = read_csv(metabolon_datakey, show_col_types = FALSE) #read in metadata file with ID's for metabolites
  ros_clin = read_csv(path_to_rosmap_clin, show_col_types = FALSE)   #read rosmap clinical metadata
  metabolon = filter(metabolon, individualID %in% ros_clin$individualID)  #keep the metabolon subjects from my rosmap data
  metabolon$projid =  ros_clin$projid[match(metabolon$individualID, ros_clin$individualID)] #add projectID key to metabolon
  metab_data = dplyr::select(metabolon, projid | matches("^\\d+$")) #keeps columns for projid and metabolites (which are numbered)
  
  
  emat = as.data.frame(t(column_to_rownames(metab_data, var = "projid"))) #make subjects cols, metabs rows
  emat = log2(emat+1)
  emat$Gene_Symbols = rownames(emat) 
  emat = dplyr::select(emat, Gene_Symbols, everything())
  
  #read in cyclops predicted ordering
  cyc_pred_file = list.files(path = paste0(path_to_cyclops_ordering, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(path_to_cyclops_ordering, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)
  cyc_pred = cyc_pred[na.exclude(match(colnames(emat), cyc_pred$ID)),] #keep the subjects in cyc_pred that I have metabolites for
  cyc_pred$pmi = ros_clin$pmi[match( cyc_pred$ID, ros_clin$projid)]
  cyc_pred$sex = ros_clin$msex[match(cyc_pred$ID, ros_clin$projid)]
  keep_subs = c(1, which(colnames(emat) %in% cyc_pred$ID))
  emat= emat[, keep_subs]  #keep the emat samples I have ordering for
  all(colnames(emat)[-1] == cyc_pred$ID)
  
  #Add ratios to check for cycling
  ratio1 = as.numeric(emat["1310",-1]) / as.numeric(emat["428",-1])
  NAD_NADH <- c(GeneName = "NAD+/NADH", ratio1)
  ratio2 = as.numeric(emat["496",-1]) / as.numeric(emat["448",-1])
  GSH_GSSG <- c(GeneName = "GSH/GSSG", ratio2)
  ratio3 = as.numeric(emat["100000265",-1]) / as.numeric(emat["565",-1])
  KYN_TRY <- c(GeneName = "KYN/TRY", ratio3)
  ratio4 = as.numeric(emat["460",-1]) / as.numeric(emat["815",-1])
  PHE_TYR <- c(GeneName = "PHE/TYR", ratio4)
  ratio5 = as.numeric(emat["1263",-1]) / as.numeric(emat["197",-1])
  SAM_SAH <- c(GeneName = "SAM/SAH", ratio5)
  BCAAsum = as.numeric(emat["376",-1]) + as.numeric(emat["397",-1]) + as.numeric(emat["566",-1])
  BCAAs <- c(GeneName = "BCAAs", BCAAsum)
  AAAsum = as.numeric(emat["565",-1]) + as.numeric(emat["460",-1]) + as.numeric(emat["815",-1])
  AAAs <- c(GeneName = "AAAs", AAAsum)
  ratio6 = BCAAsum/AAAsum
  BCAAs_AAAs = c(GeneName = "BCAAs/AAAs", ratio6)
  ratio7 = BCAAsum/as.numeric(emat["811",-1]) 
  BCAA_Ala = c(GeneName = "BCAAs/Ala", ratio7)
  ratio8 = as.numeric(emat["823",-1]) / as.numeric(emat["482",-1])
  PYRUV_LAC <- c(GeneName = "PYRUV/LAC", ratio8)
  ratio9 = as.numeric(emat["1830",-1]) / as.numeric(emat["270",-1])
  ACoA_CoA <- c(GeneName = "ACoA/CoA", ratio9)
  
  emat = rbind(c("Cond_D", cyc_pred$Covariate_D),c("sex_D",   paste0("cond_", cyc_pred$sex)),
               c("pmi_C", cyc_pred$pmi),
               emat, NAD_NADH, GSH_GSSG, KYN_TRY, PHE_TYR, SAM_SAH,
               BCAAs, AAAs, BCAAs_AAAs, BCAA_Ala,PYRUV_LAC,ACoA_CoA )

  #run regression on emat
  cycling_in_CTL = is_cycling(cyc_pred, emat, "cond_0", percentile = percentile)
  cycling_in_AD = is_cycling(cyc_pred, emat, "cond_1", percentile = percentile)
  cycling_in_AD_CTL_mthd2 = is_cycling_method2(cyc_pred, emat, useBatch = F, percentile = percentile)
  
  #select the metabolites with rhythm BHq-value < cutoff
  strong_cycling_in_CTL = dplyr::filter(cycling_in_CTL, as.numeric(BHQ) < BHQ_cutoff) %>% arrange(as.numeric(BHQ))
  strong_cycling_in_AD = dplyr::filter(cycling_in_AD, as.numeric(BHQ) < BHQ_cutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_method2 = dplyr::filter(cycling_in_AD_CTL_mthd2, as.numeric(BHQ) < BHQ_cutoff) %>% arrange(as.numeric(BHQ))
  
  #test those genes for differential rhythms
  genelist = union(strong_cycling_in_CTL$Gene_Symbols, strong_cycling_in_AD$Gene_Symbols)
  diff_rhythms = diff_rhyth(cyc_pred, emat, genelist, percentile = percentile)
  diff_rhythms_mthd2 = diff_rhyth(cyc_pred, emat, strong_cyclers_method2$Gene_Symbols, percentile = percentile)
  
  # write out files 
  out_cycling_CTL = merge(dplyr::select(metabolon_data_key, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, SHORT_NAME, CHEMSPIDER, HMDB, PUBCHEM, KEGG),
                          cycling_in_CTL, by.x = "CHEM_ID", by.y = "Gene_Symbols", all.y = T)
  write.table(out_cycling_CTL, paste0(path_to_cyclops_ordering, "/metabolon/cycling_in_CTL.csv"), col.names = T, row.names = F, sep = ',')
  
  out_cycling_AD = merge(dplyr::select(metabolon_data_key, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, SHORT_NAME, CHEMSPIDER, HMDB, PUBCHEM, KEGG),
                         cycling_in_AD, by.x = "CHEM_ID", by.y = "Gene_Symbols", all.y = T)
  write.table(out_cycling_AD, paste0(path_to_cyclops_ordering, "/metabolon/cycling_in_AD.csv"), col.names = T, row.names = F, sep = ',')
  
  isCyclingSigCutoff_str = str_extract(as.character(BHQ_cutoff), "(?<=\\.)\\d+")
  
  diff_rhythms = merge(dplyr::select(metabolon_data_key, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, SHORT_NAME, CHEMSPIDER, HMDB, PUBCHEM, KEGG),
                       diff_rhythms, by.x = "CHEM_ID", by.y = "Gene_Symbols", all.y = T)
  write.table(diff_rhythms, paste0(path_to_cyclops_ordering, "metabolon/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,".csv"), col.names = T, row.names = F, sep = ',')
  
  diff_rhythms_mthd2 = merge(dplyr::select(metabolon_data_key, CHEM_ID, SUPER_PATHWAY, SUB_PATHWAY, CHEMICAL_NAME, SHORT_NAME, CHEMSPIDER, HMDB, PUBCHEM, KEGG),
                           diff_rhythms_mthd2, by.x = "CHEM_ID", by.y = "Gene_Symbols", all.y = T)
  write.table(diff_rhythms_mthd2, paste0(path_to_cyclops_ordering, "metabolon/diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,".csv"), col.names = T, row.names = F, sep = ',')
  if (plot_mets){
    ############
    # plotting #
    ############
    setwd(paste0(path_to_cyclops_ordering, "/metabolon"))
    top_dr_metabs = diff_rhythms %>% filter(as.numeric(p_val)<0.05) %>% dplyr::select(CHEM_ID) %>% unname %>% unlist
    seedlist = c(strong_cycling_in_CTL$Gene_Symbols[1:3], strong_cycling_in_AD$Gene_Symbols[1:3], top_dr_metabs)
    seedlist = unique(seedlist)
    rownames(emat) = NULL
    print(plot_gene_trace(cyc_pred = cyc_pred, tmm = emat, seedlist = seedlist, savePlots = T, percentile = percentile, adjust_points_for_batches = F))
  }
  
}


