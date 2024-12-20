#Henry Hollis Oct 30 2024
#This script takes the seurat objects from syn53366818
#and counts the number of exc.3 and exc.5 cells for each individual
library(tidyverse)
library(Seurat)
library(SeuratDisk)

#This cell popluation has all the Exc.3 and Exc.5 neurons:
filename = c("cux2+.h5Seurat")

#Subjects with more than one batch:
subjects_of_interest = c("R2144127","R2880377","R4817881","R4996203","R5447358","R5693901",
                         "R6280004","R6679530","R7702934","R8724814","R8760165","R8781115")

#read in seurat_data
seurat_object = LoadH5Seurat(filename, assays = c("RNA"), reductions = F, graphs = F, neighbors = F, images = F)

#separate the 12 subjects with 2 batches as multiple individuals:
indiv_ids_of_dub_subs = seurat_object$individualID[which(seurat_object$individualID %in% subjects_of_interest)]
batches_of_dub_subs = seurat_object$batch[which(seurat_object$individualID %in% subjects_of_interest)]
batches_of_dub_subs_substring = str_extract(batches_of_dub_subs, "B[0-9]+")
seurat_object$individualID[which(seurat_object$individualID %in% subjects_of_interest)] = paste0(indiv_ids_of_dub_subs,"-", batches_of_dub_subs_substring)

Exc35_cell_counts_per_subject = table(seurat_object$state, seurat_object$individualID)
write.table(t(Exc35_cell_counts_per_subject), "Exc35_cell_counts_per_subject.csv", sep = ",", col.names = NA, row.names = T)