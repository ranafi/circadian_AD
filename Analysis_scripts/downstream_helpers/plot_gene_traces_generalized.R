library(tidyverse)
library(doParallel)
library(gridExtra)
library(progress)

blunt_outliers = function(vec, percentile = 0.){
  num =length(which(!is.na(vec)))
  blunt_n_points = round(percentile * num, 0)
  ord = sort(vec)
  upper_val = ord[num -blunt_n_points]
  lower_val = ord[blunt_n_points+1]
  
  vec[which(vec > upper_val)] = upper_val
  vec[which(vec < lower_val)] = lower_val
  return(vec)
}


plot_gene_trace = function(cyc_pred, tmm, seedlist,  useBatch = F,
                           percentile = 0., savePlots = F, split_cond_plots = T, adjust_points_for_batches = T, save_pdf = F){
  if(useBatch){print("NOTE: Using batches in plotting")}
  
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
  
   df = column_to_rownames(tmm, var = "Gene_Symbols") %>%
     t %>% as.data.frame %>% rownames_to_column( var = "Subjid") %>%
     merge(., cyc_pred, by.x = "Subjid", by.y = "ID") %>% arrange(Phase)

  if (useBatch){
    b = as.factor(df$Batch_D)      #batch factor
  }
  I = as.factor(df$Cond_D)  # condtion factor
  times = as.numeric(df$Phase) #in the case that I have CYCLOPS preds for subs not in tmm...
  s = as.factor(df$sex_D) #sex of each subject
  p = as.numeric(df$pmi_C) #pmi of each subject
  
  all_genes = foreach (i = 1:length(seedlist)) %do%{
    if (!(seedlist[i] %in% colnames(df))) {print(paste(seedlist[i],"not found")); return()}
    gexp1 = as.numeric(unlist(df[,seedlist[i]]))
    times1 = df$Phase
    I1 = I
    s1 = s
    p1 = p
    my_df = df
    if(useBatch){b1 = b}
    
    rm_NA = which(is.na(gexp1))
    if (length(rm_NA) <= floor(.7*nrow(df))){ #only proceed if >70% of data are not NA
      if(!is_empty(rm_NA)){ #if there are NA data points...
        gexp1 = gexp1[-rm_NA] #remove them from the expression, times, and condition vectors
        times1 = times1[-rm_NA]
        I1 = I[-rm_NA]
        s1 = s1[-rm_NA]
        p1 = p1[-rm_NA]
        my_df = my_df[-rm_NA,]
        if(useBatch){b1 = b1[-rm_NA]} #and batch vector if using batch
      }
      
      #blunt x percentile in each condition separately
      gexp1[I1==levels(I1)[1]] = blunt_outliers(gexp1[I1==levels(I1)[1]], percentile = percentile)
      gexp1[I1==levels(I1)[2]] = blunt_outliers(gexp1[I1==levels(I1)[2]], percentile = percentile)
      
      if (useBatch){
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + b1 + p1 + s1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1+ p1 + s1)
        design_matrix = model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + b1+ p1 + s1)
      }else{
        partial_model = lm(gexp1 ~ sin(times1) + cos(times1) + I1 + p1 + s1)
        full_model = lm(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1)
        design_matrix = model.matrix(gexp1 ~ I1*sin(times1) + I1*cos(times1) + I1 + p1 + s1)
      }
      anova_results = anova(partial_model, full_model)
      
      rm_coeffs = grep("sin|cos",names(full_model[["coefficients"]]))
      mesor_AD = mean(subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs])
      mesor_CTL = mean(subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs])
      
      #From the fitted values for all points we are removing the non-circadian variance. We do this separately for the AD points and CTL points and add back in the mesor for both AD and CTL:
      full_model$fitted.values[I1 == levels(I1)[1]] = full_model$fitted.values[I1 == levels(I1)[1]] - subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 0 ) %*% full_model[["coefficients"]][-rm_coeffs] + mesor_CTL 
      full_model$fitted.values[I1 == levels(I1)[2]] = full_model$fitted.values[I1 == levels(I1)[2]] - subset(design_matrix[,-rm_coeffs], design_matrix[, "I1cond_1"]== 1 ) %*% full_model[["coefficients"]][-rm_coeffs] + mesor_AD

      
      fitted_ctl = full_model$fitted.values[I1 == levels(I1)[1]]
      fitted_ad = full_model$fitted.values[I1 == levels(I1)[2]]
      
      if(adjust_points_for_batches){
        batch_offset_correct = design_matrix[,"b1cond_1"] * full_model$coefficients[["b1cond_1"]]
        gexp1 = gexp1 - batch_offset_correct
      }
      
      my_df[,seedlist[i]] =  gexp1
      y_sd = sd(my_df[,seedlist[i]])
      y_mean = mean(my_df[,seedlist[i]])
  
      # ylim_max = max(my_df[,seedlist[i]])
      # ylim_min = min(my_df[,seedlist[i]])
      df_AD = my_df %>% dplyr::filter(Cond_D == "cond_1")
      df_AD$fitted_values = fitted_ad
      df_CTL = my_df %>% dplyr::filter(Cond_D == "cond_0")
      df_CTL$fitted_values = fitted_ctl
      if(split_cond_plots){
        
        plot1 = ggplot(df_CTL , aes(x = Phase , y = df_CTL[,seedlist[i]])) +
          # ylim(ylim_min, ylim_max)+
          ylim(max(0,y_mean - 3 * y_sd), y_mean + 3 * y_sd)+
          geom_point(aes(color = "CTL")) +
          geom_line(mapping=aes(x=Phase, y=fitted_values, color = "CTL"), linetype = "solid",linewidth = 2) +
          labs(title = paste0(seedlist[i], " in CTL"), x = "Circadian Phase", y = "Expression")+
          #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
          #scale_shape_manual(values=c(2, 16))+
          scale_colour_manual(values = c("#0091ff"))+
          scale_x_continuous(breaks = seq(0, 2 * pi, by = pi/2),
                             labels = c("0", expression(pi/2), expression(pi),
                             expression(3*pi/2), expression(2*pi)))+
          theme(
            plot.title = element_text(size = 14),      # Title font size
            axis.title = element_text(size = 14),     # Axis title font size
            axis.text = element_text(size = 14),      # Axis text font size
            legend.text = element_text(size = 12),    # Legend text font size
            legend.title = element_text(size = 14)    # Legend title font size
          )
        
       
        plot2 = ggplot(df_AD , aes(x = Phase , y = df_AD[,seedlist[i]])) +
          # ylim(ylim_min, ylim_max)+
          ylim(max(0,y_mean - 3 * y_sd), y_mean + 3 * y_sd)+
          geom_point(aes(color = "AD")) +
          geom_line(mapping=aes(x=Phase, y=fitted_values, color = "AD"), linetype = "solid",linewidth = 2) +
          labs(title = paste0(seedlist[i], " in AD"), x = "Circadian Phase", y = "Expression")+
          #annotate("text", x=min(plot_df2$Phase)+.5,y=lims[1], label = paste("DR FDR", DR_FDR))+
          #scale_shape_manual(values=c(2, 16))+
          scale_colour_manual(values = c("red"))+
          scale_x_continuous(breaks = seq(0, 2 * pi, by = pi/2),
                             labels = c("0", expression(pi/2), expression(pi),
                                        expression(3*pi/2), expression(2*pi)))+
          theme(
            plot.title = element_text(size = 14),      # Title font size
            axis.title = element_text(size = 14),     # Axis title font size
            axis.text = element_text(size = 14),      # Axis text font size
            legend.text = element_text(size = 12),    # Legend text font size
            legend.title = element_text(size = 14)    # Legend title font size
          )
        plot = grid.arrange(plot1, plot2, nrow = 1)
      }else{
        plot = ggplot(df_AD, aes(x = Phase , y = df_AD[,seedlist[i]], color = "AD")) +
          geom_point(aes(color = "AD"),  alpha = 0.4) +
          geom_line(mapping=aes(x=Phase, y=fitted_values, color = "AD"), linetype = "solid", linewidth = 2) +
          geom_point(data = df_CTL, mapping = aes(x = Phase , y = df_CTL[,seedlist[i]], color = "CTL"),  alpha = 0.4) +
          geom_line(data=df_CTL, mapping=aes(x=Phase, y=fitted_values, color = "CTL"), linetype = "solid", linewidth = 2) +
          labs(title = paste0(seedlist[i]), x = "Predicted Phase", y = "Expression")+
          scale_colour_manual(values = c("red", "#0091ff"))+
          scale_x_continuous(breaks = seq(0, 2 * pi, by = pi/2),
                             labels = c("0", expression(pi/2), expression(pi),
                                        expression(3*pi/2), expression(2*pi)))+
          ylim(max(0,y_mean - 3 * y_sd), y_mean + 3 * y_sd)+
          theme(
            plot.title = element_text(size = 14),      # Title font size
            axis.title = element_text(size = 14),     # Axis title font size
            axis.text = element_text(size = 14),      # Axis text font size
            legend.text = element_text(size = 12),    # Legend text font size
            legend.title = element_text(size = 14)    # Legend title font size
          )
        
      }
      
      
      #print(p)
      if(savePlots){
        gene_name = str_replace(seedlist[i], "/", "_")
        if(split_cond_plots == F){gene_name = paste0(gene_name, "_combined_plot")}
        ggsave(paste0("plots/", gene_name, "_gene_trace.png"),plot)
        if(save_pdf){
          if(split_cond_plots){
            ggsave(paste0("plots/", gene_name, "_gene_trace.pdf"),plot= plot, width = 7.2, height = 4.5, units = "in")
          }else{
            ggsave(paste0("plots/", gene_name, "_gene_trace.pdf"),plot= plot, width = 3.6, height = 4.5, units = "in")
          }
        }
      }
    }
  }

}


plot_genes = function(tmm_path, cyclops_path, genes_to_plot = NULL,force_align = NULL, force_flipped = F, useBatch = useBatch,
                                 percentile = percentile, split_cond_plots = T, save_pdf = F){
  current_dir = getwd()
  tmm = read_csv(tmm_path, show_col_types = FALSE)
  cyc_pred_file = list.files(path = paste0(cyclops_path, "/Fits/"), pattern = '*Fit_Output_*')
  cyc_pred = read_csv(paste(cyclops_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)
  
  if(!is.null(force_align)){
    cyc_pred$Phase = (cyc_pred$Phase - force_align) %% (2*pi)
    print("Shifting phases using `force align`")
  }
  if (force_flipped){
    cyc_pred$Phase  = -cyc_pred$Phase %% (2*pi)
    print("flipping acrophases.")
  }
  setwd(paste0(cyclops_path, "/downstream_output"))
  if (!(dir.exists("plots"))){
    dir.create("plots")
  }
  
  genelist = c("ARNTL", "NPAS2", "CLOCK", "CRY1", "CRY2", "NR1D1", "NR1D2", "PER1", "PER2", "PER3", "DBP", "TEF")
  if(is.null(genes_to_plot)){
    genes_to_plot = genelist
  }
  plot_gene_trace(cyc_pred, tmm, genes_to_plot, useBatch = useBatch,
                  savePlots = T, percentile = percentile, split_cond_plots = split_cond_plots, save_pdf = save_pdf)
  setwd(current_dir)
  
}

plot_subject_histogram = function(cyclops_path, cond_subset){
  
  cyc_pred_file = list.files(path = paste0(cyclops_path, "/Fits/"), pattern = '*Fit_Output_*')
  fits = read_csv(paste(cyclops_path, "Fits", cyc_pred_file[1], sep = '/'), show_col_types = FALSE)
  setwd(cyclops_path)
  
  fits = filter(fits, Covariate_D == cond_subset)
  if (nrow(fits) > 0){
    bin = 2*pi/8
    title = ifelse(cond_subset == "cond_0", "CTL Subjects Phase Distribution", "AD Subjects Phase Distribution") 
    fill_color = ifelse(cond_subset == "cond_0", "#0091ff", "red") 
    jpeg(file=paste0("downstream_output/plots/", cond_subset, "_phase_histogram.jpg"))
    hist(fits$Phase, breaks = seq(0, 2*pi, by = pi/6), xaxt = "n", main = title, xlab = "Phase Prediction", col = fill_color)
    axis(side=1, at=c(0,pi, 2*pi),
         labels=c("0",expression(pi),expression(2*pi)))
    dev.off()
  }else{
    print(paste("No subjects matching cond_subset:", cond_subset))
  }
  
}


