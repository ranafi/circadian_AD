library(GEOquery)
library(tidyverse)
library(doParallel)

# load series and platform data from GEO
gset <- getGEO("GSE71620", GSEMatrix =TRUE, getGPL=T)
gene_assign = gset[["GSE71620_series_matrix.txt.gz"]]@featureData@data %>%
  # filter(gene_assignment !="---") %>%
  separate_wider_delim("gene_assignment",names = c("GeneBank", "Gene_Symbols", "desc"),  delim = " // ", too_many = "drop", too_few = "align_start", cols_remove = F) %>%
  dplyr::select(ID, Gene_Symbols)
gene_assign$ID = as.character(gene_assign$ID)

#which subjects are BA11?
BA11_subs = grep("11", gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[["source_name_ch1"]])
emat = gset[["GSE71620_series_matrix.txt.gz"]]@assayData[["exprs"]][,BA11_subs] #only the BA11 subjects
#Grab subject times
tod = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[["tod:ch1"]][BA11_subs] #only the BA11 subjects
tod = (as.numeric(tod) %% 24 ) *pi / 12
subs_wo_tod = which(is.na(tod))
tod = tod[-subs_wo_tod]

#expression data of 146 subjects, logged
emat = emat[, -subs_wo_tod]

#Resolve probe duplicates in same manner as paper:
# "If a gene was represented by multiple probe sets, the one with the largest intensity interquartile region was selected to represent that gene."

# ---  Prepare Data ---
expression_data <- as.data.frame(emat) %>%
    rownames_to_column(var = "ID") 

# Identify sample columns ( all columns EXCEPT the ID column)
sample_cols <- setdiff(colnames(expression_data), "ID")

# --- Calculate IQR per Probe ---
# Using apply (base R approach within the pipe)
probe_iqrs <- expression_data %>%
  # Select only sample columns for IQR calculation
  select(all_of(sample_cols)) %>% 
  # Calculate IQR for each row (probe). na.rm=TRUE is important if you have missing values
  apply(1, IQR, na.rm = TRUE) 

# Add IQRs back to the data frame with probe IDs
probe_iqr_df <- data.frame(
  ID = expression_data$ID,
  IQR = probe_iqrs
)

# ---  Merge Annotation ---
probe_iqr_annotated <- probe_iqr_df %>%
  left_join(gene_assign %>% select(ID, Gene_Symbols), by = "ID")

# --- Filter Missing/Empty Gene Symbols ---
probe_iqr_filtered <- probe_iqr_annotated %>%
  filter(!is.na(Gene_Symbols) & Gene_Symbols != "")

# --- 5. Group by Gene Symbol & 6. Find Max IQR within Group ---
selected_probes <- probe_iqr_filtered %>%
  group_by(Gene_Symbols) %>%
  # Arrange by IQR descending within each group and take the first one
   arrange(desc(IQR)) %>% 
  slice(1) %>% # Takes the top row after arranging
  ungroup() # Ungroup for further steps

# --- View the selected probes and their IQRs ---
print("Selected probes per gene based on max IQR:")
print(selected_probes)

# ---  Filter Original Matrix/Data Frame ---
# Get the vector of probe IDs to keep
probes_to_keep <- selected_probes$ID

# Filter the original matrix (assuming probe IDs are row names)
filtered_expression_matrix <- emat[probes_to_keep, ]

# --- Display the dimensions of the filtered matrix ---
print(paste("Original matrix dimensions:", paste(dim(emat), collapse = " x ")))
print(paste("Filtered matrix dimensions:", paste(dim(filtered_expression_matrix), collapse = " x "))) 
# Check that the number of rows matches the number of unique selected gene symbols
print(paste("Number of unique selected gene symbols:", length(unique(selected_probes$Gene_Symbols))))



emat_filt_logged <- cbind(Gene_Symbols = selected_probes$Gene_Symbols, 
                                  as.data.frame(filtered_expression_matrix)) 
emat_filt <- cbind(Gene_Symbols = selected_probes$Gene_Symbols, 
                   as.data.frame(2^(filtered_expression_matrix)))

pmi = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[colnames(filtered_expression_matrix), "characteristics_ch1.2"] %>%
  str_remove("pmi: ") %>% as.numeric 
print(mean(pmi))
age = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[colnames(filtered_expression_matrix), "age:ch1"] %>%
  as.numeric 
print(mean(age))
tod = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[colnames(filtered_expression_matrix), "tod:ch1"] %>%
  as.numeric 
rin = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[colnames(filtered_expression_matrix), "rin:ch1"] %>%
  as.numeric 
print(mean(rin))
sex = gset[["GSE71620_series_matrix.txt.gz"]]@phenoData@data[colnames(filtered_expression_matrix), "Sex:ch1"] 
table(sex)
# Create the mapping
lookup <- c("M" = 0, "F" = 1)
# Use the character vector to look up values in the named vector
sex_vec <- paste0("cond_", unname(lookup[sex])) # unname() removes the "M"/"F" names from the result
todRad = (tod %% 24)*pi/12
out = rbind(c("TOD_C", tod), c("TODRAD_C", todRad),  c("age_C", age), c("pmi_C",pmi), c("RIN_C", rin), c("sex_D", sex_vec), emat_filt)

write.table(out, file = "~/Desktop/BA11.csv", col.names = T, row.names = F, sep = ",")
