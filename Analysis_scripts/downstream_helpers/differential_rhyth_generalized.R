library(tidyverse)
library(doParallel)
library(progress)

blunt_outliers = function(vec, percentile =0.){
  num =length(which(!is.na(vec)))
  blunt_n_points = round(percentile * num, 0)
  ord = sort(vec)
  upper_val = ord[num -blunt_n_points]
  lower_val = ord[blunt_n_points+1]

  vec[which(vec > upper_val)] = upper_val
  vec[which(vec < lower_val)] = lower_val
  return(vec)
}

#test which genes are cycling from cyclops subject phase prediction
is_cycling = function(cyc_pred, tmm, cond_subset, pb = NULL, useBatch = F, percentile = 0.){
  cat(paste("\nRunning is_cycling() on cond_subset:", cond_subset))
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  #test significant in the following genes, here that all of them.
  genelist = unlist(unname(tmm[!grepl("_D|_C", unlist(tmm[,1])), 1])) #ASSUMES FIRST COL is names

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch, sex, pmi) %>% filter(Covariate_D == cond_subset) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, sex, pmi) %>% filter(Covariate_D == cond_subset) %>% arrange(Phase)
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # since genelist is all genes, "gene" will be tmm without gene_names
  gene = gene %>% mutate_all(as.numeric)
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  # get the transpose, subjects x genes and put in order of CYCLOPS order
  colnames(gene1) =  unname(unlist(tmm[which(unname(unlist(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1

  #below 2 lines I use "match" in case I am given phases for subjects not in the data
  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) #the batch variable
  }
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #the given phase of each subject
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)])
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)])
  
  #loop:
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    s1 = s
    p1 = p
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }
      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      
      if (useBatch){
          partial_model = lm(gexp1 ~ b1 + p1 + s1)
          full_model = lm(gexp1 ~ sin(times1) + cos(times1)+ b1 + p1 + s1)
          design_matrix <- model.matrix(gexp1 ~ sin(times1) + cos(times1)+ b1 + p1 + s1)
          
        }else{
          partial_model = lm(gexp1 ~ p1 + s1)
          full_model = lm(gexp1 ~ sin(times1) + cos(times1) + p1 + s1)
          design_matrix <- model.matrix(gexp1 ~ sin(times1) + cos(times1)+ p1 + s1)

      }

      anova_results = anova(partial_model, full_model)
      sin_coff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      acrophase = atan2(sin_coff, cos_coeff) %% (2*pi)

      p_statistic = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      amplitude = sqrt(sin_coff^2 + cos_coeff^2)
        # batch_weighted_mesor = (sum(p1*full_model[["coefficients"]][["p1"]]) +(     full_model[["coefficients"]][["(Intercept)"]] * sum(b1 == levels(b1)[1] & s1 == levels(s1)[1]) + # num_B1S1* intercept
        #                             (full_model[["coefficients"]][["b1cond_1"]] + full_model[["coefficients"]][["(Intercept)"]])* sum(b1 == levels(b1)[2] & s1 == levels(s1)[1]) + # num_B2S1 * (int + B_offset)
        #                              (full_model[["coefficients"]][["s1cond_1"]] + full_model[["coefficients"]][["(Intercept)"]])*sum(b1 == levels(b1)[1] & s1 == levels(s1)[2]) + #num_B1S2 * (int + S_offset)
        #                               (full_model[["coefficients"]][["b1cond_1"]] + full_model[["coefficients"]][["s1cond_1"]] + full_model[["coefficients"]][["(Intercept)"]])* sum(b1 == levels(b1)[2] & s1 == levels(s1)[2])  ))/ #num_B2S2 * (int + B_offset + S_offset)
        #                           length(b1)
        
        batch_weighted_mesor = mean(design_matrix[,-c(2, 3)] %*% full_model[["coefficients"]][-c(2, 3)]) #tested, gives same results as above.
        
        amp_ratio = amplitude / batch_weighted_mesor

     
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }
      gene_summary = cbind( Gene_Symbols, acrophase,amplitude, p_statistic, amp_ratio, sin_coff, cos_coeff)

      return(gene_summary)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA))

  }
  all_genes = as_tibble(all_genes) %>% drop_na
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_statistic), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_statistic), "bonferroni")
  return(all_genes)

}

is_cycling_method2 = function(cyc_pred, tmm, pb = NULL, useBatch = F, percentile = 0.){
  cat("\nRunning is_cycling_method2() (analogue of compareRhythms) on all subjects")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  genelist = unlist(unname(tmm[!grepl("_D|_C", unlist(tmm[,1])), 1])) #ASSUMES FIRST COL is names

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch, pmi, sex) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, pmi, sex) %>% arrange(Phase)
  }


  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # since genelist is all genes, "gene" will be tmm without gene_names
  gene = gene %>% mutate_all(as.numeric)
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  # get the transpose, subjects x genes and put in order of CYCLOPS order
  colnames(gene1) =  unname(unlist(tmm[which(unname(unlist(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1

  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)])      #batch factor
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)]) # CTL or AD factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])  #in the case that I have CYCLOPS preds for subs not in tmm...
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I1 = I
    s1 = s
    p1 = p
    if(useBatch){b1 = b}
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)

      if (useBatch){
        partial_model = lm(gexp1 ~ I1 + b1 + p1 + s1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1 + p1 + s1)
        design_matrix <- model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1 + p1 + s1)
      }else{
        partial_model = lm(gexp1 ~ I1 + p1 + s1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1 )
        design_matrix <- model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1)
      }

      anova_results = anova(partial_model, full_model)

      p_val = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I1cond_1:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I1cond_1:cos(times1)"]] + cos_coeff
      acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
      acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
      amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))

      
      rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
      mesor_AD = mean(subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs]) #tested, gives same results as above.
      mesor_CTL = mean(subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs]) #tested, gives same results as above.
     
      amp_ratio_CTL = amplitude_CTL/ mesor_CTL
      amp_ratio_AD = amplitude_AD/ mesor_AD

      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }

      info = cbind( Gene_Symbols, p_val, acrophase_AD, acrophase_CTL, amplitude_AD, amplitude_CTL, amp_ratio_CTL, amp_ratio_AD, mesor_CTL, mesor_AD)
      return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA, NA, NA, NA))
  }


  all_genes = as_tibble(all_genes)
  all_genes$p_val = as.numeric(all_genes$p_val)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_val), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_val), "bonferroni")
  all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)


}

diff_rhyth = function(cyc_pred, tmm, genelist,  pb = NULL, useBatch = F, percentile = 0.){
  cat(paste("\nRunning diff_rhyth() on genelist of size:", length(genelist)))
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch, pmi, sex) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, pmi, sex) %>% arrange(Phase)
  }


  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # "gene" is tmm with only genelist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1

  if (useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)])      #batch factor
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])  # condtion factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I1 = I
    s1 = s
    p1 = p
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I1 = I1[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)

      if (useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + b1 + p1 + s1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1 + p1 + s1)
        design_matrix <- model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1 + p1 + s1)
      }else{
        # partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + 0)
        # full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + 0)
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + p1 + s1 )
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1)
        design_matrix <- model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1 )
      }
      anova_results = anova(partial_model, full_model)

      p_val = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I1cond_1:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I1cond_1:cos(times1)"]] + cos_coeff
      acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
      acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
      amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }

      info = cbind( Gene_Symbols, p_val, acrophase_AD, acrophase_CTL, amplitude_AD, amplitude_CTL)
      return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA))
  }


  all_genes = as_tibble(all_genes)
  all_genes$p_val = as.numeric(all_genes$p_val)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_val), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_val), "bonferroni")
  all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)

}

diff_rhyth_AD_severity = function(cyc_pred, tmm, genelist, rosmap_clin_path,  pb = NULL, useBatch = F, percentile = 0.){
  cat("\nRunning diff_rhyth_AD_severity()")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}
  ##### read in ROSMAP clin ####
  rosmap_clin = read_csv(rosmap_clin_path, show_col_types = FALSE)
  rosmap_clin = rosmap_clin[ na.exclude(match(cyc_pred$ID, rosmap_clin$projid)),]
  rosmap_clin = rosmap_clin %>%
    mutate(braaksc_bin = cut(braaksc, c(0, 3, 5, 7), right = F))
  rosmap_clin = rosmap_clin %>%
    mutate(ceradsc_bin = cut(ceradsc, c(1, 3, 5), right = F))
  # rosmap_clin$apoe_ordinal  = 1
  # rosmap_clin$apoe_ordinal[rosmap_clin$apoe_genotype == 34 | rosmap_clin$apoe_genotype == 24] = 2
  # rosmap_clin$apoe_ordinal[rosmap_clin$apoe_genotype == 44 ] = 3
  # rosmap_clin$apoe_ordinal[is.na(rosmap_clin$apoe_genotype) ] = NA
  ###############

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    cyc_pred_merged = merge(cyc_pred, dplyr::select(rosmap_clin, !pmi), by.x = "ID", by.y = "projid", y.keep = F)
    # preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin,  braaksc_bin, batch, pmi, sex) %>% arrange(Phase)
    preds = cyc_pred_merged %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, braaksc, braaksc_bin, batch, pmi, sex) %>% arrange(Phase)
    
  }else{
    cyc_pred_merged = merge(cyc_pred, dplyr::select(rosmap_clin, !pmi), by.x = "ID", by.y = "projid", y.keep = F)
    # preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, braaksc_bin, pmi, sex) %>% arrange(Phase)
    preds = cyc_pred_merged %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, braaksc, braaksc_bin, pmi, sex) %>% arrange(Phase)
    
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # "gene" is tmm with only genelist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1

  cog = as.factor(preds$cogdx[match(rownames(gene1), preds$ID)])      # cogdx score 1, 2, 4 or 5
  cerad = as.factor(preds$ceradsc_bin[match(rownames(gene1), preds$ID)]) #cerad score 2 bins
  braak = as.factor(preds$braaksc_bin[match(rownames(gene1), preds$ID)]) #Tau score 3 bins
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  
  if (useBatch){b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) }

  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I_local_cog = cog
    I_local_cerad = cerad
    I_local_braak = braak
    s1 = s
    p1 = p
    if(useBatch){b1 = b}

    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I_local_cog = cog[-rm_NA]
        I_local_cerad = cerad[-rm_NA]
        I_local_braak = braak[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }

      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      
      ##below is code for testing the cognitive variable##
    
      if(useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + b1 + p1 +s1)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1 + p1 + s1)
        # design_matrix <- model.matrix(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1 + p1 + s1)
      
      }else{
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + p1 + s1)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + p1 + s1)
        # design_matrix <- model.matrix(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + p1 + s1)
        
      }

      anova_results = anova(partial_model, full_model)
      p_cog = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]

      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I_local_cog2:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I_local_cog2:cos(times1)"]] + cos_coeff
      sin_coeff4 = full_model[["coefficients"]][["I_local_cog4:sin(times1)"]] + sin_coeff
      cos_coeff4 = full_model[["coefficients"]][["I_local_cog4:cos(times1)"]] + cos_coeff
      sin_coeff5 = full_model[["coefficients"]][["I_local_cog5:sin(times1)"]] + sin_coeff
      cos_coeff5 = full_model[["coefficients"]][["I_local_cog5:cos(times1)"]] + cos_coeff
      
      acrophase_cog1 = atan2(sin_coeff, cos_coeff) %% (2*pi)
      amplitude_cog1= sqrt((sin_coeff^2) + (cos_coeff^2))
      acrophase_cog2 = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_cog2= sqrt((sin_coeff2^2) + (cos_coeff2^2))
      acrophase_cog4 = atan2(sin_coeff4, cos_coeff4) %% (2*pi)
      amplitude_cog4= sqrt((sin_coeff4^2) + (cos_coeff4^2))
      acrophase_cog5 = atan2(sin_coeff5, cos_coeff5) %% (2*pi)
      amplitude_cog5= sqrt((sin_coeff5^2) + (cos_coeff5^2))
      ####### ceradsc_binned #########
      if(useBatch){
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + b1 + p1 + s1)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + b1 + p1 + s1)
      }else{
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + p1 + s1)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + p1 + s1)
      }
      anova_results_cerad = anova(partial_model_cerad, full_model_cerad)
      p_cerad = anova_results_cerad$`Pr(>F)`[2]

      sin_coeff_cerad = full_model_cerad[["coefficients"]][["sin(times1)"]]
      cos_coeff_cerad = full_model_cerad[["coefficients"]][["cos(times1)"]]
      sin_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):sin(times1)"]] + sin_coeff_cerad
      cos_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):cos(times1)"]] + cos_coeff_cerad
      acrophase_cerad1to2 = atan2(sin_coeff_cerad, cos_coeff_cerad) %% (2*pi)
      amplitude_cerad1to2 = sqrt((sin_coeff_cerad^2) + (cos_coeff_cerad^2))
      acrophase_cerad3to5 = atan2(sin_coeff2_cerad, cos_coeff2_cerad) %% (2*pi)
      amplitude_cerad3to5= sqrt((sin_coeff2_cerad^2) + (cos_coeff2_cerad^2))
      
      ######Braaksc binned #####
      if(useBatch){
        partial_model_braak = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_braak + b1 + p1 + s1)
        full_model_braak = lm(gexp1 ~ I_local_braak*sin(times1) + I_local_braak*cos(times1) + I_local_braak + b1 + p1 + s1)
      }else{
        partial_model_braak = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_braak + p1 + s1)
        full_model_braak = lm(gexp1 ~ I_local_braak*sin(times1) + I_local_braak*cos(times1) + I_local_braak + p1 + s1)
      }
      anova_results_braak = anova(partial_model_braak, full_model_braak)
      p_braak = anova_results_braak$`Pr(>F)`[2]
      
      sin_coeff_braak = full_model_braak[["coefficients"]][["sin(times1)"]] #the base class is [0,3)
      cos_coeff_braak = full_model_braak[["coefficients"]][["cos(times1)"]]
      sin_coeff2_braak = full_model_braak[["coefficients"]][["I_local_braak[3,5):sin(times1)"]] + sin_coeff_braak # [3,5) class
      cos_coeff2_braak = full_model_braak[["coefficients"]][["I_local_braak[3,5):cos(times1)"]] + cos_coeff_braak
      sin_coeff3_braak = full_model_braak[["coefficients"]][["I_local_braak[5,7):sin(times1)"]] + sin_coeff_braak # [5,7) class
      cos_coeff3_braak = full_model_braak[["coefficients"]][["I_local_braak[5,7):cos(times1)"]] + cos_coeff_braak
      
      acrophase_braak0to3 = atan2(sin_coeff_braak, cos_coeff_braak) %% (2*pi)
      amp_braak0to3 = sqrt((sin_coeff_braak^2) + (cos_coeff_braak^2))
      acrophase_braak3to5 = atan2(sin_coeff2_braak, cos_coeff2_braak) %% (2*pi)
      amp_braak3to5 = sqrt((sin_coeff2_braak^2) + (cos_coeff2_braak^2))
      acrophase_braak5to7 = atan2(sin_coeff3_braak, cos_coeff3_braak) %% (2*pi)
      amp_braak5to7 = sqrt((sin_coeff3_braak^2) + (cos_coeff3_braak^2))
      
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }

      info = c( Gene_Symbols, p_cog, p_cerad, p_braak, acrophase_cog1, acrophase_cog2, acrophase_cog4, acrophase_cog5,
                acrophase_cerad1to2, acrophase_cerad3to5, acrophase_braak0to3, acrophase_braak3to5, acrophase_braak5to7, amplitude_cog1, 
                amplitude_cog2, amplitude_cog4,amplitude_cog5, amplitude_cerad1to2, amplitude_cerad3to5, amp_braak0to3, amp_braak3to5, amp_braak5to7)
      return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  }


  all_genes = as_tibble(all_genes)
  colnames(all_genes) = c("Gene_Symbols", "p_cogdx","p_ceradsc","p_braak","acrophase_cog1", "acrophase_cog2", "acrophase_cog4", "acrophase_cog5",
                          "acrophase_cerad1to2", "acrophase_cerad3to5", "acrophase_braak0to3", "acrophase_braak3to5", "acrophase_braak5to7" ,
                          "amplitude_cog1", "amplitude_cog2" ,"amplitude_cog4",
                          "amplitude_cog5", "amplitude_cerad1to2", "amplitude_cerad3to5", "amplitude_braak0to3", "amplitude_braak3to5", "amplitude_braak5to7")
  all_genes$BHQ_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "BH")
  all_genes$Bonf_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "bonferroni")
  all_genes$BHQ_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "BH")
  all_genes$Bonf_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "bonferroni")
  all_genes$BHQ_braak = p.adjust(as.numeric(all_genes$p_braak), "BH")
  all_genes$Bonf_braak = p.adjust(as.numeric(all_genes$p_braak), "bonferroni")
  #all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)

}
diff_rhyth_AD_severity_AD_only = function(cyc_pred, tmm, genelist, rosmap_clin_path,  pb = NULL, useBatch = F, percentile = 0.){
  cat("\nRunning diff_rhyth_AD_severity() on just the cogdx > 3")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}
  ##### read in ROSMAP clin ####
  rosmap_clin = read_csv(rosmap_clin_path, show_col_types = FALSE)
  rosmap_clin = rosmap_clin[ na.exclude(match(cyc_pred$ID, rosmap_clin$projid)),]
  rosmap_clin = rosmap_clin %>%
    mutate(braaksc_bin = cut(braaksc, c(0, 3, 5, 7), right = F))
  rosmap_clin = rosmap_clin %>%
    mutate(ceradsc_bin = cut(ceradsc, c(1, 3, 5), right = F))
  
  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    cyc_pred_merged = merge(cyc_pred, dplyr::select(rosmap_clin, !pmi), by.x = "ID", by.y = "projid", y.keep = F)
    preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, braaksc, braaksc_bin, batch, pmi, sex) %>% arrange(Phase)

  }else{
    cyc_pred_merged = merge(cyc_pred, dplyr::select(rosmap_clin, !pmi), by.x = "ID", by.y = "projid", y.keep = F)
    preds = cyc_pred_merged %>% dplyr::filter(Covariate_D == "cond_1") %>% dplyr::select(ID, Phase, cogdx, ceradsc_bin, braaksc, braaksc_bin, pmi, sex) %>% arrange(Phase)

  }
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # "gene" is tmm with only genelist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1
  
  cog = as.factor(preds$cogdx[match(rownames(gene1), preds$ID)])      # cogdx score 1, 2, 4 or 5
  cerad = as.factor(preds$ceradsc_bin[match(rownames(gene1), preds$ID)]) #cerad score 2 bins
  braak = as.factor(preds$braaksc_bin[match(rownames(gene1), preds$ID)]) #Tau score 3 bins
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)]) #in the case that I have CYCLOPS preds for subs not in tmm...
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  
  if (useBatch){b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) }
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    I_local_cog = cog
    I_local_cerad = cerad
    I_local_braak = braak
    s1 = s
    p1 = p
    if(useBatch){b1 = b}
    
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        I_local_cog = cog[-rm_NA]
        I_local_cerad = cerad[-rm_NA]
        I_local_braak = braak[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }
      
      gexp1 = blunt_outliers(gexp1, percentile = percentile)
      
      ##below is code for testing the cognitive variable##
      
      if(useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + b1 + p1 +s1)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1 + p1 + s1)
        # design_matrix <- model.matrix(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + b1 + p1 + s1)
        
      }else{
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cog + p1 + s1)
        full_model = lm(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + p1 + s1)
        # design_matrix <- model.matrix(gexp1 ~ I_local_cog*sin(times1) + I_local_cog*cos(times1) + I_local_cog + p1 + s1)
        
      }
      
      anova_results = anova(partial_model, full_model)
      p_cog = anova_results$`Pr(>F)`[2]
      Gene_Symbols = colnames(gene1)[gene_i]
      
      sin_coeff = full_model[["coefficients"]][["sin(times1)"]]
      cos_coeff = full_model[["coefficients"]][["cos(times1)"]]
      sin_coeff2 = full_model[["coefficients"]][["I_local_cog5:sin(times1)"]] + sin_coeff
      cos_coeff2 = full_model[["coefficients"]][["I_local_cog5:cos(times1)"]] + cos_coeff
      acrophase_cog4 = atan2(sin_coeff, cos_coeff) %% (2*pi)
      amplitude_cog4= sqrt((sin_coeff^2) + (cos_coeff^2))
      acrophase_cog5 = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
      amplitude_cog5= sqrt((sin_coeff2^2) + (cos_coeff2^2))
      ####### ceradsc_binned #########
      if(useBatch){
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + b1 + p1 + s1)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + b1 + p1 + s1)
      }else{
        partial_model_cerad = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_cerad + p1 + s1)
        full_model_cerad = lm(gexp1 ~ I_local_cerad*sin(times1) + I_local_cerad*cos(times1) + I_local_cerad + p1 + s1)
      }
      anova_results_cerad = anova(partial_model_cerad, full_model_cerad)
      p_cerad = anova_results_cerad$`Pr(>F)`[2]
      
      sin_coeff_cerad = full_model_cerad[["coefficients"]][["sin(times1)"]]
      cos_coeff_cerad = full_model_cerad[["coefficients"]][["cos(times1)"]]
      sin_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):sin(times1)"]] + sin_coeff_cerad
      cos_coeff2_cerad = full_model_cerad[["coefficients"]][["I_local_cerad[3,5):cos(times1)"]] + cos_coeff_cerad
      acrophase_cerad1to2 = atan2(sin_coeff_cerad, cos_coeff_cerad) %% (2*pi)
      amplitude_cerad1to2 = sqrt((sin_coeff_cerad^2) + (cos_coeff_cerad^2))
      acrophase_cerad3to5 = atan2(sin_coeff2_cerad, cos_coeff2_cerad) %% (2*pi)
      amplitude_cerad3to5= sqrt((sin_coeff2_cerad^2) + (cos_coeff2_cerad^2))
      
      ######Braaksc binned #####
      if(useBatch){
        partial_model_braak = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_braak + b1 + p1 + s1)
        full_model_braak = lm(gexp1 ~ I_local_braak*sin(times1) + I_local_braak*cos(times1) + I_local_braak + b1 + p1 + s1)
      }else{
        partial_model_braak = lm(gexp1 ~ sin(times1) + cos(times1) + I_local_braak + p1 + s1)
        full_model_braak = lm(gexp1 ~ I_local_braak*sin(times1) + I_local_braak*cos(times1) + I_local_braak + p1 + s1)
      }
      anova_results_braak = anova(partial_model_braak, full_model_braak)
      p_braak = anova_results_braak$`Pr(>F)`[2]
      
      sin_coeff_braak = full_model_braak[["coefficients"]][["sin(times1)"]] #the base class is [0,3)
      cos_coeff_braak = full_model_braak[["coefficients"]][["cos(times1)"]]
      sin_coeff2_braak = full_model_braak[["coefficients"]][["I_local_braak[3,5):sin(times1)"]] + sin_coeff_braak # [3,5) class
      cos_coeff2_braak = full_model_braak[["coefficients"]][["I_local_braak[3,5):cos(times1)"]] + cos_coeff_braak
      sin_coeff3_braak = full_model_braak[["coefficients"]][["I_local_braak[5,7):sin(times1)"]] + sin_coeff_braak # [5,7) class
      cos_coeff3_braak = full_model_braak[["coefficients"]][["I_local_braak[5,7):cos(times1)"]] + cos_coeff_braak
      
      acrophase_braak0to3 = atan2(sin_coeff_braak, cos_coeff_braak) %% (2*pi)
      amp_braak0to3 = sqrt((sin_coeff_braak^2) + (cos_coeff_braak^2))
      acrophase_braak3to5 = atan2(sin_coeff2_braak, cos_coeff2_braak) %% (2*pi)
      amp_braak3to5 = sqrt((sin_coeff2_braak^2) + (cos_coeff2_braak^2))
      acrophase_braak5to7 = atan2(sin_coeff3_braak, cos_coeff3_braak) %% (2*pi)
      amp_braak5to7 = sqrt((sin_coeff3_braak^2) + (cos_coeff3_braak^2))
      
      if (!is.null(pb)){
        if(!pb$finished){
          pb$tick()
        }
      }
      
      info = c( Gene_Symbols, p_cog, p_cerad, p_braak, acrophase_cog4, acrophase_cog5,
                acrophase_cerad1to2, acrophase_cerad3to5, acrophase_braak0to3, acrophase_braak3to5, acrophase_braak5to7, 
                amplitude_cog4,amplitude_cog5, amplitude_cerad1to2, amplitude_cerad3to5, amp_braak0to3, amp_braak3to5, amp_braak5to7)
      return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
  }
  
  
  all_genes = as_tibble(all_genes)
  colnames(all_genes) = c("Gene_Symbols", "p_cogdx","p_ceradsc","p_braak","acrophase_cog1", "acrophase_cog2", "acrophase_cog4", "acrophase_cog5",
                          "acrophase_cerad1to2", "acrophase_cerad3to5", "acrophase_braak0to3", "acrophase_braak3to5", "acrophase_braak5to7" ,
                          "amplitude_cog1", "amplitude_cog2" ,"amplitude_cog4",
                          "amplitude_cog5", "amplitude_cerad1to2", "amplitude_cerad3to5", "amplitude_braak0to3", "amplitude_braak3to5", "amplitude_braak5to7")
  all_genes$BHQ_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "BH")
  all_genes$Bonf_cogdx = p.adjust(as.numeric(all_genes$p_cogdx), "bonferroni")
  all_genes$BHQ_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "BH")
  all_genes$Bonf_cerad = p.adjust(as.numeric(all_genes$p_ceradsc), "bonferroni")
  all_genes$BHQ_braak = p.adjust(as.numeric(all_genes$p_braak), "BH")
  all_genes$Bonf_braak = p.adjust(as.numeric(all_genes$p_braak), "bonferroni")
  #all_genes$Log_AD_CTL_ampRatio = log(as.numeric(all_genes$amplitude_AD) / as.numeric(all_genes$amplitude_CTL))
  return(all_genes)
  
}
mesor_differences = function(cyc_pred, tmm, genelist, pb = NULL, useBatch = F, percentile = 0.){ ##
  cat("\nRunning Mesor_differences()")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}

  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch, pmi, sex) %>% arrange(Phase)

  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, pmi, sex) %>% arrange(Phase)
  }

  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # "gene" is tmm with only genelist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1

  if(useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) # batch
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])  # condition factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  

  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    times1 = times
    s1 = s
    p1 = p
    I1 = I
    if(useBatch){b1 = b}
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(gene1))){ #only proceed if >70% of data are not NA
      
      if(!is_empty(rm_NA)){
        gexp1 = gexp1[-rm_NA]
        times1 = times1[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        I1 = I1[-rm_NA]
        if(useBatch){b1 = b1[-rm_NA]}
      }
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
    if(useBatch){
      partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + b1 + p1 + s1)
      full_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + b1 + p1 + s1)
      design_mat = model.matrix(gexp1 ~ sin(times1) + cos(times1) + I1 + b1 + p1 + s1)
    }else{
      partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + p1 + s1)
      full_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + p1 + s1)
      design_mat = model.matrix(gexp1 ~ sin(times1) + cos(times1) + I1 + p1 + s1)
    }

    anova_results = anova(partial_model, full_model)
    # wilcox_test = wilcox.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    # p_wilcox = wilcox_test$p.value

    # t_test = t.test(gexp1[I == levels(I)[1]], gexp1[I == levels(I)[2]])
    # p_ttest = t_test$p.value

    p_mesor = anova_results$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    
    rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
    mesor_AD = mean(subset(design_mat[,-rm_coeffs], design_mat[, "I1cond_1"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    mesor_CTL = mean(subset(design_mat[,-rm_coeffs], design_mat[, "I1cond_1"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    
    # sin_coeff = full_model1[["coefficients"]][["sin(times)"]]
    # cos_coeff = full_model1[["coefficients"]][["cos(times)"]]
    # sin_coeff2 = full_model1[["coefficients"]][["Icond_1:sin(times)"]] + sin_coeff
    # cos_coeff2 = full_model1[["coefficients"]][["Icond_1:cos(times)"]] + cos_coeff
    # acrophase_CTL = atan2(sin_coeff, cos_coeff) %% (2*pi)
    # acrophase_AD = atan2(sin_coeff2, cos_coeff2) %% (2*pi)
    # amplitude_CTL = sqrt((sin_coeff^2) + (cos_coeff^2))
    # amplitude_AD = sqrt((sin_coeff2^2) + (cos_coeff2^2))
    if (!is.null(pb)){
      if(!pb$finished){
        pb$tick()
      }
    }

    info = cbind( Gene_Symbols, p_mesor, mesor_CTL, mesor_AD)
    return(info)
    }
    return(cbind( colnames(gene1)[gene_i], NA,NA, NA))
    
  }
  all_genes = as_tibble(all_genes)
  all_genes$BHQ = p.adjust(as.numeric(all_genes$p_mesor), "BH")
  all_genes$Bonf = p.adjust(as.numeric(all_genes$p_mesor), "bonferroni")
  # all_genes$BHQ_wilcox = p.adjust(as.numeric(all_genes$p_wilcox), "BH")
  # all_genes$BHQ_ttest = p.adjust(as.numeric(all_genes$p_ttest), "BH")

  return(all_genes)
}

diff_mesor_AD_severity = function(cyc_pred, tmm, genelist, pb = NULL, useBatch = F, percentile = 0.){ ##
  cat("\nRunning diff_mesor_AD_severity()")
  if(useBatch){cat("\nNOTE: Using batches in regression.")}
  
  cond_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "cond_d")
  cyc_pred$Covariate_D = tmm[cond_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  sex_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "sex_d")
  cyc_pred$sex = tmm[sex_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  pmi_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "pmi_c")
  cyc_pred$pmi = tmm[pmi_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  cerad_bin_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "ceradsc_d")
  cyc_pred$ceradsc =  tmm[cerad_bin_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  braak_bin_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "braaksc_d")
  cyc_pred$braaksc =  tmm[braak_bin_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
  
  if (useBatch){
    batch_row_of_tmm = which(tolower(unlist(tmm[, 1])) == "batch_d")
    cyc_pred$batch = tmm[batch_row_of_tmm, na.exclude(match(cyc_pred$ID, colnames(tmm)))] %>% unname %>% unlist
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, batch, pmi, sex, ceradsc, braaksc) %>% arrange(Phase)
    
  }else{
    preds= dplyr::select(cyc_pred, ID, Covariate_D, Phase, pmi, sex, ceradsc, braaksc) %>% arrange(Phase)
  }
  
  gene = tmm[which(unlist(unname(tmm[,1])) %in% genelist), -1] # "gene" is tmm with only genelist subset
  gene1 = t(gene[,na.exclude(match(preds$ID, colnames(gene)))])  #the transpose, subjects x genes for tidyverse purposes
  colnames(gene1) =  unname(unlist(tmm[which(unlist(unname(tmm[,1])) %in% genelist), 1]))  #add the gene names to the columns of gene1
  
  if(useBatch){
    b = as.factor(preds$batch[match(rownames(gene1), preds$ID)]) # batch
  }
  I = as.factor(preds$Covariate_D[match(rownames(gene1), preds$ID)])  # condition factor
  times = as.numeric(preds$Phase[match(rownames(gene1), preds$ID)])
  s = as.factor(preds$sex[match(rownames(gene1), preds$ID)]) #sex of each subject
  p = as.numeric(preds$pmi[match(rownames(gene1), preds$ID)]) #pmi of each subject
  cerad = as.factor(preds$ceradsc[match(rownames(gene1), preds$ID)])
  braak = as.factor(preds$braaksc[match(rownames(gene1), preds$ID)])
  
  all_genes = foreach (gene_i = 1:ncol(gene1), .combine = rbind) %do%{
    gexp1 = as.numeric(unlist(gene1[,gene_i]))
    gexp1 = blunt_outliers(gexp1, percentile = percentile)
    
    #cogdx stratification ----------------------
    if(useBatch){
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + b + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + I + b + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + I + b + p + s)
    }else{
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + I + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + I + p + s)
    }
    
    anova_results = anova(partial_model, full_model)
    p_mesor_cogdx = anova_results$`Pr(>F)`[2]
    Gene_Symbols = colnames(gene1)[gene_i]
    
    rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
    #Find the individuals with AD and the mean of the non-circadian variance (fitted full_model without circadian terms)
    mesor_AD_cogdx = mean(subset(design_mat[,-rm_coeffs], design_mat[, "Icond_1"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    #Find the individuals with CTL and the mean of the non-circadian variance (fitted full_model without circadian terms)
    mesor_CTL_cogdx = mean(subset(design_mat[,-rm_coeffs], design_mat[, "Icond_1"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
 
    #Cerad stratification -------------------------
    if(useBatch){
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + b + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + cerad + b + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + cerad + b + p + s)
    }else{
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + cerad + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + cerad + p + s)
    }
    
    anova_results = anova(partial_model, full_model)
    p_mesor_cerad = anova_results$`Pr(>F)`[2]
    rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
    mesor_cerad_bin1 = mean(subset(design_mat[,-rm_coeffs], design_mat[, "ceradcond_[3,5)"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    mesor_cerad_bin2 = mean(subset(design_mat[,-rm_coeffs], design_mat[, "ceradcond_[3,5)"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    
    #Braak stratification ------------------------
    if(useBatch){
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + b + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + braak + b + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + braak + b + p + s)
    }else{
      partial_model = lm(gexp1 ~ sin(times) + cos(times) + p + s)
      full_model = lm(gexp1 ~ sin(times) + cos(times) + braak + p + s)
      design_mat = model.matrix(gexp1 ~ sin(times) + cos(times) + braak + p + s)
    }
    
    anova_results = anova(partial_model, full_model)
    p_mesor_braak = anova_results$`Pr(>F)`[2]
    rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
    
    mesor_braak_bin1 = mean(subset(design_mat[,-rm_coeffs], (design_mat[, "braakcond_[3,5)"]== 0 & design_mat[, "braakcond_[5,7)"]== 0) ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    mesor_braak_bin2 = mean(subset(design_mat[,-rm_coeffs], (design_mat[, "braakcond_[3,5)"]== 1 & design_mat[, "braakcond_[5,7)"]== 0) ) %*% full_model[["coefficients"]][-rm_coeffs]) 
    mesor_braak_bin3 = mean(subset(design_mat[,-rm_coeffs], (design_mat[, "braakcond_[3,5)"]== 0 & design_mat[, "braakcond_[5,7)"]== 1) ) %*% full_model[["coefficients"]][-rm_coeffs]) 
      
    if (!is.null(pb)){
      if(!pb$finished){
        pb$tick()
      }
    }
    
    info = cbind( Gene_Symbols, p_mesor_cogdx, mesor_CTL_cogdx, mesor_AD_cogdx,
                  p_mesor_cerad, mesor_cerad_bin1, mesor_cerad_bin2, p_mesor_braak,
                  mesor_braak_bin1, mesor_braak_bin2, mesor_braak_bin3)
    return(info)
  }
  all_genes = as_tibble(all_genes)
  all_genes$BHQ_cogdx = p.adjust(as.numeric(all_genes$p_mesor_cogdx), "BH")
  all_genes$Bonf_cogdx = p.adjust(as.numeric(all_genes$p_mesor_cogdx), "bonferroni")
  all_genes$BHQ_cerad = p.adjust(as.numeric(all_genes$p_mesor_cerad), "BH")
  all_genes$BHQ_braak = p.adjust(as.numeric(all_genes$p_mesor_braak), "BH")
  
  return(all_genes)
}


##### main function #####

run_cycling_and_dr_analysis = function(order_path, tmm_path, rosmap_clin_path, isCyclingSigCutoff = 0.05, useBatch = F, percentile = 0.){
  tmm = read_csv(tmm_path, show_col_types = FALSE)      #read expression data, unordered
  colnames(tmm)[1] = "gene_names" #set first column name bc sometimes they are different

  #get cyclops subject predictions
  cyc_pred_file = list.files(path = paste0(order_path, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(order_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)

  #perform is_cycling (method 1) nested regression on CTL data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_CTL = is_cycling(cyc_pred, tmm, cond_subset = "cond_0", useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  # Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_CTL$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # cycling_in_CTL = cbind(Ensembl, cycling_in_CTL)
  #record strong cyclers in CTL, for several different amplitude cutoffs
  strong_cyclers_CTL_AR20 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.20 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_CTL_AR33 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.33 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_CTL_AR1 = dplyr::filter(cycling_in_CTL, as.numeric(amp_ratio) >=0.1 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  #perform is_cycling (method 1) nested regression on AD data
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_AD = is_cycling(cyc_pred, tmm, cond_subset = "cond_1", useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  # Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_AD$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # cycling_in_AD = cbind(Ensembl, cycling_in_AD)
  #record strong cyclers in AD, for several different amplitude cutoffs
  strong_cyclers_AD_AR20 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.20 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_AD_AR33 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.33 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_AD_AR1 = dplyr::filter(cycling_in_AD, as.numeric(amp_ratio) >=0.1 & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  #perform is_cycling (method 2) nested regression on all data ( all data tested together, same to compareRhythms)
  pb <- progress_bar$new(total = dim(tmm)[1])
  cycling_in_either_cond = is_cycling_method2(cyc_pred, tmm, useBatch = useBatch, pb = pb, percentile = percentile)
  #Add Ensemble genenames to output
  # Ensembl = Ensembl_dict$ENSEMBL[match(cycling_in_either_cond$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # cycling_in_either_cond = cbind(Ensembl, cycling_in_either_cond)
  #record strong cyclers in either condition, for several different amplitude cutoffs
  strong_cyclers_method2_AR20 = dplyr::filter(cycling_in_either_cond, ( (as.numeric(amp_ratio_CTL) >=0.20) | (as.numeric(amp_ratio_AD) >=0.20)) & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))
  strong_cyclers_method2_AR1 = dplyr::filter(cycling_in_either_cond, ( (as.numeric(amp_ratio_CTL) >=0.1) | (as.numeric(amp_ratio_AD) >=0.1)) & as.numeric(BHQ) < isCyclingSigCutoff) %>% arrange(as.numeric(BHQ))

  # We only test for diff rhythmicity if a gene cycles in AD OR CTL. Here I create those unions
  genelist_AR20 = union(strong_cyclers_AD_AR20$Gene_Symbols, strong_cyclers_CTL_AR20$Gene_Symbols)
  genelist_AR33 = union(strong_cyclers_AD_AR33$Gene_Symbols, strong_cyclers_CTL_AR33$Gene_Symbols)
  genelist_AR1 = union(strong_cyclers_AD_AR1$Gene_Symbols, strong_cyclers_CTL_AR1$Gene_Symbols)
  genelist_method2_AR20 = strong_cyclers_method2_AR20$Gene_Symbols
  genelist_method2_AR1 = strong_cyclers_method2_AR1$Gene_Symbols

  DR_genelist_list = list(genelist_AR20, genelist_AR33, genelist_AR1, genelist_method2_AR20, genelist_method2_AR1)

  ##### mesor differences ######
  gene_list_mesor =  unlist(unname(tmm[!grepl("_D|_C", unlist(tmm[,1])), 1])) # TEST ALL genes for Mesor diff (not just cyclers)
  pb <- progress_bar$new(total = length(gene_list_mesor))
  differential_mesor = mesor_differences(cyc_pred, tmm, gene_list_mesor,useBatch = useBatch, pb = pb, percentile = percentile)
  # Ensembl = Ensembl_dict$ENSEMBL[match(differential_mesor$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # differential_mesor = cbind(Ensembl, differential_mesor)
  
  pb <- progress_bar$new(total = length(gene_list_mesor))
  mesor_ad_severity = diff_mesor_AD_severity(cyc_pred, tmm, gene_list_mesor,useBatch = useBatch, pb = pb, percentile = percentile)

  ####### differential rhtyhms with continuous cerad covs####
  if(length(genelist_AR20) != 0){
  diff_rhythms_AD_severity = diff_rhyth_AD_severity(cyc_pred, tmm,
                                                    # unname(unlist(strong_cyclers_AD_AR20$Gene_Symbols)),
                                                    genelist_AR20,
                                                    rosmap_clin_path = rosmap_clin_path,
                                                    percentile = percentile, useBatch = useBatch)
  diff_rhythms_AD_severity_AD_only = diff_rhyth_AD_severity_AD_only(cyc_pred,
                                 tmm,
                                 genelist_AR20,
                                 rosmap_clin_path = rosmap_clin_path,
                                 percentile = percentile,
                                 useBatch = useBatch)
  }else{
    print("Not running AD_severity functions, genelist_AR20 is empty")
  }
  ####### differential rhythms #####
  DR_results = list()
  for(genelist in DR_genelist_list){
    if (length(genelist) == 0) {
      # Append NULL for empty genelist to preserve alignment
      DR_results <- c(DR_results, list(NULL))
    } else {
    pb <- progress_bar$new(total = length(genelist))
    diff_rhythms_results = diff_rhyth(cyc_pred, tmm, genelist, useBatch = useBatch, pb = pb, percentile = percentile)
    # Ensembl = Ensembl_dict$ENSEMBL[match(diff_rhythms_results$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
    # diff_rhythms_results = cbind(Ensembl, diff_rhythms_results)
    DR_results = c(DR_results, list(diff_rhythms_results))
    }
  }

  diff_rhythms20 = DR_results[[1]]
  diff_rhythms33 = DR_results[[2]]
  diff_rhythms1 = DR_results[[3]]
  diff_rhythms_mthd2_AR20 = DR_results[[4]]
  diff_rhythms_mthd2_AR1 = DR_results[[5]]

  #Create list of strong cyclers (AR 0.20 or 0.33 and BHQ < BHQcutoff) in CTL subjects
  CTL_cyclers_AR20BHQCutoff = dplyr::select(strong_cyclers_CTL_AR20 , Gene_Symbols )
  CTL_cyclers_AR33BHQCutoff = dplyr::select(strong_cyclers_CTL_AR33 , Gene_Symbols )

  #Create list of strong cyclers (AR 0.20 or 0.33 and BHQ < BHQcutoff) in AD subjects
  AD_cyclers_AR20BHQCutoff = dplyr::select(strong_cyclers_AD_AR20 , Gene_Symbols )
  AD_cyclers_AR33BHQCutoff = dplyr::select(strong_cyclers_AD_AR33 , Gene_Symbols )

  #Create list of strong cyclers (AR 0.20 or 0.33 and BHQ < BHQcutoff) in Either (method 2)
  mthd2_cyclers_AR20BHQCutoff = dplyr::select(strong_cyclers_method2_AR20, Gene_Symbols)
  mthd2_cyclers_AR1BHQCutoff = dplyr::select(strong_cyclers_method2_AR1, Gene_Symbols)

  # All genes expressed in CTL and AD:
  EnrichR_background = dplyr::select(cycling_in_CTL, Gene_Symbols)

  null_DR_elements <- sapply(DR_results, is.null)
  if(any(!(null_DR_elements))){
    #Create list of strong DR genes (AR 0.33 or 0.20 or 0.1, and BHQ< 0.3 
    DR_cyclers_AR33_DRBHQ3 = dplyr::filter(diff_rhythms33,  as.numeric(BHQ) < 0.3) %>% arrange(as.numeric(BHQ))
    DR_cyclers_AR20_DRBHQ3 = dplyr::filter(diff_rhythms20,  as.numeric(BHQ) < 0.3) %>% arrange(as.numeric(BHQ))
    DR_cyclers_AR1_DRBHQ3 = dplyr::filter(diff_rhythms1,  as.numeric(BHQ) < 0.3) %>% arrange(as.numeric(BHQ))
    DR_cyclers_mthd2_AR20_DRBHQ3 = dplyr::filter(diff_rhythms_mthd2_AR20,  as.numeric(BHQ) < 0.3) %>% arrange(as.numeric(BHQ))
    DR_cyclers_mthd2_AR1_DRBHQ3 = dplyr::filter(diff_rhythms_mthd2_AR1,  as.numeric(BHQ) < 0.3) %>% arrange(as.numeric(BHQ))
  

    #Create list of lost cycling DR genes:
    DR_lostAmpAD_AR33BHQ3 = filter(diff_rhythms33, BHQ < 0.3, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select( Gene_Symbols)
    DR_lostAmpAD_AR20BHQ3 = filter(diff_rhythms20, BHQ < 0.3, Log_AD_CTL_ampRatio < 0)  %>% dplyr::select( Gene_Symbols)
    DR_lostAmpAD_AR1BHQ3 = filter(diff_rhythms1, BHQ < 0.3, Log_AD_CTL_ampRatio < 0) %>% dplyr::select( Gene_Symbols)
    DR_mthd2_lostAmpAD_AR20BHQ3 = filter(diff_rhythms_mthd2_AR20, BHQ < 0.3, Log_AD_CTL_ampRatio < 0) %>% dplyr::select( Gene_Symbols)
    DR_mthd2_lostAmpAD_AR1BHQ3 = filter(diff_rhythms_mthd2_AR1, BHQ < 0.3, Log_AD_CTL_ampRatio < 0) %>% dplyr::select( Gene_Symbols)
  
    #Create lists of gained cycling DR genes
    DR_gainAmpAD_AR33BHQ3= filter(diff_rhythms33, BHQ < 0.3, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select( Gene_Symbols)
    DR_gainAmpAD_AR20BHQ3= filter(diff_rhythms20, BHQ < 0.3, Log_AD_CTL_ampRatio > 0)  %>% dplyr::select( Gene_Symbols)
    DR_gainAmpAD_AR1BHQ3 = filter(diff_rhythms1, BHQ < 0.3, Log_AD_CTL_ampRatio > 0) %>% dplyr::select( Gene_Symbols)
    DR_mthd2_gainAmpAD_AR20BHQ3 = filter(diff_rhythms_mthd2_AR20, BHQ < 0.3, Log_AD_CTL_ampRatio > 0) %>% dplyr::select( Gene_Symbols)
    DR_mthd2_gainAmpAD_AR1BHQ3 = filter(diff_rhythms_mthd2_AR1, BHQ < 0.3, Log_AD_CTL_ampRatio > 0) %>% dplyr::select( Gene_Symbols)

  }
  
  if (!(dir.exists(paste(order_path, "downstream_output", sep = '/')))){
    dir.create(paste(order_path, "downstream_output", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_results", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "enrichR_files", sep = '/'))
    dir.create(paste(order_path, "downstream_output", "PSEA_files", sep = '/'))

  }
  #Create string for isCyclingSigCutoff.  E.g. 0.05 -> "05"
  isCyclingSigCutoff_str = str_extract(as.character(isCyclingSigCutoff), "(?<=\\.)\\d+")
  blunting_percentile_str = str_extract(as.character(percentile), "(?<=\\.)\\d+")

  #write out all results of cycling and DR analysis
  if(any(!(null_DR_elements))){
    write.table(diff_rhythms20, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio20.csv"), sep = ',', row.names = F, col.names = T)
    write.table(diff_rhythms33, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio33.csv"), sep = ',', row.names = F, col.names = T)
    write.table(diff_rhythms1, paste0(order_path, "/downstream_output/diff_rhythms_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio1.csv"), sep = ',', row.names = F, col.names = T)
    write.table(diff_rhythms_mthd2_AR20, paste0(order_path, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio20.csv"), sep = ',', row.names = F, col.names = T)
    write.table(diff_rhythms_mthd2_AR1, paste0(order_path, "/downstream_output/diff_rhythms_method2_CyclingBHQ",isCyclingSigCutoff_str,"AmpRatio1.csv"), sep = ',', row.names = F, col.names = T)
  }
  write.table(cycling_in_CTL, paste(order_path, "downstream_output","cosinor_results_CTL.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_AD, paste(order_path, "downstream_output","cosinor_results_AD.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  write.table(cycling_in_either_cond, paste(order_path, "downstream_output","cosinor_results_method2_cyclingInEither.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #gene lists for enrichR
  write.table(CTL_cyclers_AR20BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR20BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(CTL_cyclers_AR33BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/CTL_cyclers_AR33BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR20BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR20BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(AD_cyclers_AR33BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_cyclers_AR33BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)
  write.table(mthd2_cyclers_AR20BHQCutoff, paste0(order_path, "/downstream_output/enrichR_files/AD_CTL_mthd2_cyclers_AR20BHQ", isCyclingSigCutoff_str, ".csv"), sep = ',', row.names = F, col.names = T)

  #background genes for enrichR
  write.table(EnrichR_background, paste(order_path, "downstream_output", "enrichR_files","EnrichR_background.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
  #DR genes for enrichR
  if(any(!(null_DR_elements))){
    write.table(DR_cyclers_AR33_DRBHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_cyclers_AR20_DRBHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_cyclers_AR1_DRBHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_cyclers_mthd2_AR20_DRBHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_mthd2_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_cyclers_mthd2_AR1_DRBHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_cyclers_mthd2_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
  
    write.table(DR_lostAmpAD_AR33BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_lostAmpAD_AR20BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_lostAmpAD_AR1BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_mthd2_lostAmpAD_AR20BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_mthd2_lostAmpAD_AR1BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_lostAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
  
    write.table(DR_gainAmpAD_AR33BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR33_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_gainAmpAD_AR20BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_gainAmpAD_AR1BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_mthd2_gainAmpAD_AR20BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR20_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
    write.table(DR_mthd2_gainAmpAD_AR1BHQ3, paste0(order_path, "/downstream_output/enrichR_files/DR_mthd2_gainAmpAD_CyclingBHQ",isCyclingSigCutoff_str,"AR1_DRBHQ3.csv"), sep = ',', row.names = F, col.names = T)
  }
  
  #Mesor differences
  write.table(differential_mesor, paste(order_path, "downstream_output", "differential_mesor_all_genes.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  sig_diff_mesor = filter(differential_mesor, as.numeric(BHQ) < 0.05 ) %>% dplyr::select(  Gene_Symbols)
  write.table(sig_diff_mesor, paste(order_path, "downstream_output","enrichR_files","diff_mesor_all_genes_BHQ05.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  write.table(mesor_ad_severity, paste(order_path, "downstream_output", "differential_mesor_all_genes_cogdx_cerad_braak.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  
  #Continuous AD differences
  if(length(genelist_AR20) != 0){ #If AD_severity was tested...
  write.table(diff_rhythms_AD_severity_AD_only,paste(order_path, "downstream_output", "diff_rhythms_AD_severity_AR20_AD_only.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  
  write.table(diff_rhythms_AD_severity, paste(order_path, "downstream_output", "diff_rhythms_AD_severity_AR20.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  strong_cogdx_diffs = filter(diff_rhythms_AD_severity, BHQ_cogdx< 0.1) %>% dplyr::select(Gene_Symbols)
  # Ensembl = Ensembl_dict$ENSEMBL[match(strong_cogdx_diffs$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # strong_cogdx_diffs = cbind(Ensembl, strong_cogdx_diffs) %>% dplyr::select(Ensembl , Gene_Symbols )
  write.table(strong_cogdx_diffs, paste(order_path, "downstream_output", "enrichR_files", "strong_cogdx_diffs_AR20.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  strong_braak_diffs = filter(diff_rhythms_AD_severity, BHQ_braak < 0.1) %>% dplyr::select(Gene_Symbols)
  # Ensembl = Ensembl_dict$ENSEMBL[match(strong_braak_diffs$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # strong_braak_diffs = cbind(Ensembl, strong_braak_diffs) %>% dplyr::select(Ensembl , Gene_Symbols )
  write.table(strong_braak_diffs, paste(order_path, "downstream_output", "enrichR_files", "strong_braak_diffs_AR20.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  strong_cerad_diffs = filter(diff_rhythms_AD_severity, BHQ_cerad< 0.1) %>% dplyr::select(Gene_Symbols)
  # Ensembl = Ensembl_dict$ENSEMBL[match(strong_cerad_diffs$Gene_Symbols, Ensembl_dict$Gene_Symbol)]
  # strong_cerad_diffs = cbind(Ensembl, strong_cerad_diffs) %>% dplyr::select(Ensembl , Gene_Symbols )
  write.table(strong_cerad_diffs, paste(order_path, "downstream_output", "enrichR_files", "strong_cerad_diffs_AR20.csv", sep = "/"), sep = ',', row.names = F, col.names = T)
  }
  
  #create lists of genes for PSEA
  PSEA_CTL_cyclers_AR20BHQCutoff = strong_cyclers_CTL_AR20 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_CTL_cyclers_AR20BHQCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_CTL_cyclers_AR20BHQ", isCyclingSigCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  PSEA_AD_cyclers_AR20BHQCutoff = strong_cyclers_AD_AR20 %>% dplyr::select(Gene_Symbols, acrophase ) %>% mutate(acrophase = as.numeric(acrophase) * 12 / pi)
  write.table(PSEA_AD_cyclers_AR20BHQCutoff, paste0(order_path, "/downstream_output/PSEA_files/PSEA_AD_cyclers_AR20BHQ", isCyclingSigCutoff_str, ".txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  if(any(!(null_DR_elements))){
     PSEA_DR_AR20_acrodiffs = diff_rhythms20 %>% mutate(acro_diff = (as.numeric(acrophase_AD) - as.numeric(acrophase_CTL))*12/pi ) %>%
    dplyr::select(Gene_Symbols, acro_diff)
    write.table(PSEA_DR_AR20_acrodiffs, paste0(order_path, "/downstream_output/PSEA_files/PSEA_DR_AR20_acrodiffs.txt"), sep = '\t', row.names = F, col.names = F, quote = F)
  }
  #write out nice summary of cycling and DR genes
  if(any(!(null_DR_elements))){
  summary = data.frame(List = c("TMM name","Using batch in Regression", "isCyclingBHQCutoff", "Blunting_Percentile" ,paste0("CTL_cyclers_AR20BHQ", isCyclingSigCutoff_str), paste0("CTL_cyclers_AR33BHQ", isCyclingSigCutoff_str),
      paste0("AD_cyclers_AR20BHQ", isCyclingSigCutoff_str), paste0("AD_cyclers_AR33BHQ", isCyclingSigCutoff_str), paste0("AD_CTL_mthd2_cyclers_AR20BHQ", isCyclingSigCutoff_str),
      "DR_cyclers_AR1BHQ3", "DR_cyclers_AR20BHQ3", "DR_cyclers_AR33BHQ3", "DR_cyclers_mthd2_AR1_DRBHQ3", "DR_cyclers_mthd2_AR20_DRBHQ3",
      "DR_lostAmpAD_AR1BHQ3", "DR_lostAmpAD_AR20BHQ3", "DR_lostAmpAD_AR33BHQ3", "DR_mthd2_lostAmpAD_AR1BHQ3", "DR_mthd2_lostAmpAD_AR20BHQ3",
      "DR_gainAmpAD_AR1BHQ3", "DR_gainAmpAD_AR20BHQ3", "DR_gainAmpAD_AR33BHQ3", "DR_mthd2_gainAmpAD_AR1BHQ3", "DR_mthd2_gainAmpAD_AR20BHQ3"),
      Num_genes = c(tmm_path, useBatch, isCyclingSigCutoff, percentile, dim(CTL_cyclers_AR20BHQCutoff)[1], dim(CTL_cyclers_AR33BHQCutoff)[1],
      dim(AD_cyclers_AR20BHQCutoff)[1],dim(AD_cyclers_AR33BHQCutoff)[1],
      dim(mthd2_cyclers_AR20BHQCutoff)[1],
      dim(DR_cyclers_AR1_DRBHQ3)[1],dim(DR_cyclers_AR20_DRBHQ3)[1], dim(DR_cyclers_AR33_DRBHQ3)[1], dim(DR_cyclers_mthd2_AR1_DRBHQ3)[1], dim(DR_cyclers_mthd2_AR20_DRBHQ3)[1],
      dim(DR_lostAmpAD_AR1BHQ3)[1], dim(DR_lostAmpAD_AR20BHQ3)[1],dim(DR_lostAmpAD_AR33BHQ3)[1], dim(DR_mthd2_lostAmpAD_AR1BHQ3)[1], dim(DR_mthd2_lostAmpAD_AR20BHQ3)[1],
      dim(DR_gainAmpAD_AR1BHQ3)[1], dim(DR_gainAmpAD_AR20BHQ3)[1], dim(DR_gainAmpAD_AR33BHQ3)[1], dim(DR_mthd2_gainAmpAD_AR1BHQ3)[1], dim(DR_mthd2_gainAmpAD_AR20BHQ3)[1])
  )
  write.table(summary, paste(order_path, "downstream_output","cosinor_DR_summary.csv", sep = '/'), sep = ',', row.names = F, col.names = T)
}

}
