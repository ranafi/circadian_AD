#Henry Hollis Oct 22 2024
#This script takes the seurat objects from syn53366818
#and creates a subcluster pseudobulk matrix and metadata for
#excitatory, inhibitory, microglia, and astrocytes.
#The result is meant to be fed to select_best_subclusters.qmd
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(Matrix)

files = c("~/Downloads/microglia.h5Seurat", "~/Downloads/astrocytes.h5Seurat",
          "~/Downloads/inhibitory.h5Seurat", "~/Downloads/cux2-.h5Seurat", "~/Downloads/cux2+.h5Seurat",
          "~/Downloads/oligodendroglia.h5Seurat", "~/Downloads/vascular.niche.h5Seurat")
# Function to outer join an arbitrary number of matrices
#Need this b/c gene names will not match b/w cell types
#but we dont just want the intersect of the genes, we want the union
outer_join_matrices <- function(matrices) {
  # Get unique row names from all matrices
  all_rows <- unique(unlist(lapply(matrices, rownames)))
  # Calculate total number of columns
  total_cols <- sum(sapply(matrices, ncol))
  
  # Initialize the combined matrix
  combined_mat <- Matrix(0, nrow = length(all_rows), ncol = total_cols, sparse = T)
  rownames(combined_mat) <- all_rows
  
  # Track the current column index in the combined matrix
  current_col <- 1
  
  # Fill in values from each matrix
  for (mat in matrices) {
    # Get the number of columns in the current matrix
    ncol_mat <- ncol(mat)
    
    # Get the row names current matrix
    mat_row_names <- rownames(mat)
    
    # Fill in the combined matrix with values from the current matrix
    combined_mat[match(mat_row_names, all_rows), current_col:(current_col + ncol_mat - 1)] <- mat
    
    # Update the current column index
    current_col <- current_col + ncol_mat
  }
  
  # Set column names to reflect the original matrices
  colnames(combined_mat) <- unlist(lapply(matrices, colnames))
  
  return(combined_mat)
}
###############################
### Test outer join function: #
# # Example usage with multiple dgCMatrix matrices
# # Create example matrices
# mat1 <- Matrix(c(1, 2, 0, 0, 3, 4), nrow = 3, ncol = 2, sparse = TRUE, 
#                dimnames = list(c("A", "B", "C"), c("Value1", "Value2")))
# mat2 <- Matrix(c(0, 0, 7, 8, 0, 0, 9, 10), nrow = 4, ncol = 2, sparse = TRUE, 
#                dimnames = list(c("B", "C", "D", "E"), c("Value3", "Value4")))
# mat3 <- Matrix(c(0, 0, 0, 0, 13, 14, 15, 16), nrow = 4, ncol = 2, sparse = TRUE, 
#                dimnames = list(c("C", "D", "F", "X"), c("Value5", "Value6")))
# 
# # Put matrices in a list
# mat_list <- list(mat1, mat2, mat3)
# result = outer_join_matrices(mat_list)
#################################


# Create empty list of matrices
matrices = list()
# Create an empty dataframe to start
cell_info = data.frame()
# pseudobulk_meta that we will build:
pseudobulk_meta = data.frame()

#Subjects with more than one batch:
subjects_of_interest = c("R2144127","R2880377","R4817881","R4996203","R5447358","R5693901",
                        "R6280004","R6679530","R7702934","R8724814","R8760165","R8781115")


i = 1
for (seurat_file in files){
  print(seurat_file)
  #read in seurat_data
  seurat_object = LoadH5Seurat(seurat_file,assays = c("RNA"), reductions = F, graphs = F, neighbors = F, images = F)
  
  #separate the 12 subjects with 2 batches as multiple individuals:
  indiv_ids_of_dub_subs = seurat_object$individualID[which(seurat_object$individualID %in% subjects_of_interest)]
  batches_of_dub_subs = seurat_object$batch[which(seurat_object$individualID %in% subjects_of_interest)]
  batches_of_dub_subs_substring = str_extract(batches_of_dub_subs, "B[0-9]+")
  seurat_object$individualID[which(seurat_object$individualID %in% subjects_of_interest)] = paste0(indiv_ids_of_dub_subs,"-", batches_of_dub_subs_substring)
  
  cell_info_rows  = data.frame(individualID = seurat_object$individualID, nCount_RNA = seurat_object$nCount_RNA,
                        nFeature_RNA = seurat_object$nFeature_RNA, batch = seurat_object$batch,
                        cluster = seurat_object$subset, subcluster = seurat_object$state)
  
  # cell_info_row = c(seurat_file, length(seurat_object@meta.data$subset))
  cell_info = rbind(cell_info, cell_info_rows)
  
  #Sum counts per individual per subcluster "state"
  seurat_pseudobulk= AggregateExpression(seurat_object,
                                         assay = "RNA",
                                         group.by = c("individualID", "state"), 
                                         return.seurat = T)
 
 
  

  #Retrieve cluster information to add to pseudobulk_meta
  cluster_info = seurat_object@meta.data[, c("individualID", "state", "subset")]
  cluster_info = cluster_info %>% unite(col = "united_info", individualID, state, sep = "_") %>% distinct
  #retrieve broadclass info
  cluster_info_sorted = cluster_info[match(seurat_pseudobulk@meta.data$orig.ident, cluster_info$united_info),]
  #assert metadata is in same order as pseudobulk metadata:
  stopifnot(all(cluster_info_sorted$united_info == seurat_pseudobulk$orig.ident))
  cluster_info_sorted = cluster_info_sorted %>% separate_wider_delim(united_info, "_", 
                                                                     too_few = "error", too_many = "error", names = c("individualID", "state") )
  #Add metadata to pseudobulk_meta
  pseudobulk_meta = rbind(pseudobulk_meta, cluster_info_sorted)
  #Add counts matrix to list
  matrices[[i]] = seurat_pseudobulk@assays$RNA$counts
  
  i = i +1
}

#Add column names to dataframes
# colnames(cell_info) = c("file", "num_cells")
colnames(pseudobulk_meta) = c("individualID", "state", "broadclass")


# Perform the outer join on pseudobulk mats
result = outer_join_matrices(matrices)

##############################
# #test matrix for code below:
# mat_data <- matrix(c(1, 0, 0, 0, 0, 0,
#                      2, 1, 0, 0, 0, 0,
#                      0, 0, 1, 1, 0, 0,
#                      0, 0, 2, 1, 0, 0), nrow = 4, byrow = TRUE)
# # Convert the matrix to a sparse dgCMatrix format
# sparse_mat <- Matrix(mat_data, sparse = TRUE)
# 
# # Assign duplicate column names
# colnames(sparse_mat) <- c("A", "B", "A", "C", "D", "B")
#############################

##################SLow version of removing dups ################
# # Get the unique column names and duplicated column names
# col_name_groups = split(seq_along(colnames(result)), colnames(result))
# # Identify groups where there are duplicate column names
# duplicate_groups = col_name_groups[sapply(col_name_groups, length) > 1]
# 
# # Create a new dgCMatrix with the number of unique columns (to store the results)
# unique_colnames <- names(col_name_groups)
# n_rows <- nrow(result)
# n_cols <- length(unique_colnames)
# 
# # Initialize an empty sparse matrix to hold combined results
# combined_result <- Matrix(0, nrow = n_rows, ncol = n_cols, sparse = TRUE)
# colnames(combined_result) <- unique_colnames
# 
# 
#  # For duplicate groups, use matrix multiplication to sum the columns
#  for (i in seq_along(unique_colnames)) {
#    col_name <- unique_colnames[i]
#    group <- col_name_groups[[col_name]]
#    
#    # Sum the columns that correspond to the same name and assign to combined_result
#    combined_result[, i] <- result[, group, drop = FALSE] %*% rep(1, length(group))
#    # Disable buffering to ensure immediate output
#    if(i %%10 == 0){
#      flush.console()
#      print(paste0(i,"/",length(unique_colnames) ," unique names done"))
#      flush.console()
#    }
#  }
#
# # Remove the corresponding rows from the metadata_df
# pseudobulk_meta = pseudobulk_meta[columns_to_keep, , drop = FALSE]  # Drop the rows corresponding to removed columns
#######################################

# Function to sum duplicated column names and return filtered metadata
sum_duplicate_columns_with_metadata <- function(mat, metadata_df) {
  # Step 1: Identify column name groups (duplicated and non-duplicated)
  col_name_groups <- split(seq_along(colnames(mat)), colnames(mat))
  unique_colnames <- names(col_name_groups)
  n_cols <- length(unique_colnames)
  
  # Step 2: Create the indicator matrix for group summation
  indicator_matrix <- Matrix(0, ncol(mat), n_cols, sparse = TRUE)
  for (i in seq_along(unique_colnames)) {
    cols_in_group <- col_name_groups[[unique_colnames[i]]]
    indicator_matrix[cols_in_group, i] <- 1
  }
  print("made it here")
  # Step 3: Perform matrix multiplication to sum all duplicate columns
  combined_result <- mat %*% indicator_matrix
  
  # Step 4: Set the correct column names in the resulting matrix
  colnames(combined_result) <- unique_colnames
  
  # Step 5: Filter the metadata
  # Keep only the first occurrence of each column in the metadata
  first_occurrence_indices <- sapply(col_name_groups, `[`, 1)
  filtered_metadata <- metadata_df[first_occurrence_indices, , drop = FALSE]
  
  # Return the combined matrix and the filtered metadata
  return(list(combined_result = combined_result, filtered_metadata = filtered_metadata))
}

# Call the function to sum duplicate columns and filter the metadata
output <- sum_duplicate_columns_with_metadata(result, pseudobulk_meta)
counts= output[[1]]
pseudobulk_meta = output[[2]]

# Ensure no duplicates remain
print(any(duplicated(colnames(counts))))  # This checks if any duplicates remain


save(counts, pseudobulk_meta, cell_info, file = "rosmap437_pseudobulk_by_cluster_all_celltypes.rda")
