### These are the analysis scripts for performing downstream analysis on data after CYCLOPS algorithm is run.

### The most important file in this directory is "*Perform_downstream_analysis.qmd*." This notebook runs all the functions (whose code is in downstream_helpers directory).

The other qmd notebooks in this directory perform smaller tasks:
-   *create_subcluster_psudobulk* takes single cell files from AD Knowledge portal and processes them into rda object for select_best_scROSMAPsubcluster.qmd

-   *select_best_scROSMAPsubclusters.qmd* processes rda object from ROSMAP into TMMs that go into the "normed_counts" directory. In addition, it is responsible for finding the subclusters of the cell types with the "purest circadian signal," as assessed by a correlation analysis (deltaccd library). It also performs differential expression analysis on the normalized data, which is stored in the "edgeR_diff_expression" directory.

-   *downstream_on_diff_expression.qmd* takes the differential expression results (which is not the primary focus of this work. We are interested in the differential rhythmicity analysis, which is similar but incorporates circadian variation in genes) and leverages the downstream analysis scripts in downstream_helpers (i.e. fGSEA and enrichR).

-   *MetaboAnalystsHelper.qmd* contains simple code to match the Metabolon data to Metaboanalyst names, which is used for MetaboAnalyst pathway analysis and so on.

-   *PSEA_accross_celltypes.qmd* uses the PSEA downstream functionality to look at differences in phase of genes in different cell types.

-   *Anchor_phases_to_excitatory* takes the predicted acrophase of ARNTL in the trusted excitatory neuron subset and subtracts the Phase column of the cyclops fit file to anchor the phases.

-   *CountExc35CellsPerSubject* simply counts number of excitatory 3 & 5 neurons in each subject

### 
