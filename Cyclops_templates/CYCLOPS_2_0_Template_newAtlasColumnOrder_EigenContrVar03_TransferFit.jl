using DataFrames, Statistics, StatsBase, LinearAlgebra, MultivariateStats, PyPlot, Distributed, Random, CSV, Revise, Distributions, Dates, MultipleTesting

base_path = joinpath(homedir(), "Box Sync", "Henry_stuff", "AD_project", "scROSMAP", "scrosmap_newAtlas", "scROSMAP_project") # define the base path in which the data file, the seed gene list, and the sample collection time file are located
data_path = joinpath(base_path, "normed_counts") # path for data folder in base folder
output_path = base_path # real training run output folder in base folder

TMM = CSV.read(joinpath(data_path,  "ExcSubtypes35_FiltByExprDefault_TMM.csv"), DataFrame) # load data file in data folder

TMM = TMM[Not(3:6), :]
nth_seed_cutoff = 10000

seed_genes = unique(CSV.read(joinpath(base_path, "seed_gene_list_nonUnique.csv"), DataFrame, header = false)).Column1
seed_genes= String.(seed_genes)
# if ((length(sample_ids_with_collection_times)+length(sample_collection_times))>0) && (length(sample_ids_with_collection_times) != length(sample_collection_times))
#     error("ATTENTION REQUIRED! Number of sample ids provided (\'sample_ids_with_collection_times\') must match number of collection times (\'sample_collection_times\').")
# end

# make changes to training parameters, if required. Below are the defaults for the current version of cyclops.
training_parameters = Dict(:regex_cont => r".*_C",			# What is the regex match for continuous covariates in the data file
:regex_disc => r".*_D",							# What is the regex match for discontinuous covariates in the data file

:blunt_percent => 0.975, 						# What is the percentile cutoff below (lower) and above (upper) which values are capped

:seed_min_CV => 0.14, 							# The minimum coefficient of variation a gene of interest may have to be included in eigen gene transformation
:seed_max_CV => 0.7, 							# The maximum coefficient of a variation a gene of interest may have to be included in eigen gene transformation
:seed_mth_Gene => nth_seed_cutoff, 						# The minimum mean a gene of interest may have to be included in eigen gene transformation

:norm_gene_level => true, 						# Does mean normalization occur at the seed gene level
:norm_disc => false, 							# Does batch mean normalization occur at the seed gene level
:norm_disc_cov => 1, 							# Which discontinuous covariate is used to mean normalize seed level data

:eigen_reg => true, 							# Does regression again a covariate occur at the eigen gene level
:eigen_reg_disc_cov => 1, 						# Which discontinous covariate is used for regression
:eigen_reg_exclude => false,					# Are eigen genes with r squared greater than cutoff removed from final eigen data output
:eigen_reg_r_squared_cutoff => 0.6,				# This cutoff is used to determine whether an eigen gene is excluded from final eigen data used for training
:eigen_reg_remove_correct => false,				# Is the first eigen gene removed (true --> default) or it's contributed variance of the first eigne gene corrected by batch regression (false)

:eigen_first_var => false, 						# Is a captured variance cutoff on the first eigen gene used
:eigen_first_var_cutoff => 0.85, 				# Cutoff used on captured variance of first eigen gene

:eigen_total_var => 0.85, 						# Minimum amount of variance required to be captured by included dimensions of eigen gene data
:eigen_contr_var => 0.03, 						# Minimum amount of variance required to be captured by a single dimension of eigen gene data
:eigen_var_override => false,					# Is the minimum amount of contributed variance ignored
:eigen_max => 6, 								# Maximum number of dimensions allowed to be kept in eigen gene data

:out_covariates => true, 						# Are covariates included in eigen gene data
:out_use_disc_cov => true,						# Are discontinuous covariates included in eigen gene data
:out_all_disc_cov => true, 						# Are all discontinuous covariates included if included in eigen gene data
:out_disc_cov => 1,								# Which discontinuous covariates are included at the bottom of the eigen gene data, if not all discontinuous covariates
:out_use_cont_cov => false,						# Are continuous covariates included in eigen data
:out_all_cont_cov => true,						# Are all continuous covariates included in eigen gene data
:out_use_norm_cont_cov => false,				# Are continuous covariates Normalized
:out_all_norm_cont_cov => true,					# Are all continuous covariates normalized
:out_cont_cov => 1,								# Which continuous covariates are included at the bottom of the eigen gene data, if not all continuous covariates, or which continuous covariates are normalized if not all
:out_norm_cont_cov => 1,						# Which continuous covariates are normalized if not all continuous covariates are included, and only specific ones are included

:init_scale_change => true,						# Are scales changed
:init_scale_1 => false,							# Are all scales initialized such that the model sees them all as having scale 1
                                                # Or they'll be initilized halfway between 1 and their regression estimate.

:train_n_models => 80, 							# How many models are being trained
:train_μA => 0.001, 							# Learning rate of ADAM optimizer
:train_β => (0.9, 0.999), 						# β parameter for ADAM optimizer
:train_min_steps => 1500, 						# Minimum number of training steps per model
:train_max_steps => 2050, 						# Maximum number of training steps per model
:train_μA_scale_lim => 1000, 					# Factor used to divide learning rate to establish smallest the learning rate may shrink to
:train_circular => false,						# Train symmetrically
:train_collection_times => false,						# Train using known times
:train_collection_time_balance => 1.0,					# How is the true time loss rescaled
# :train_sample_id => sample_ids_with_collection_times,
# :train_sample_phase => sample_collection_times,

:cosine_shift_iterations => 192,				# How many different shifts are tried to find the ideal shift
:cosine_covariate_offset => true,				# Are offsets calculated by covariates

:align_p_cutoff => 0.05,						# When aligning the acrophases, what genes are included according to the specified p-cutoff
:align_base => "radians",						# What is the base of the list (:align_acrophases or :align_phases)? "radians" or "hours"
:align_disc => false,							# Is a discontinuous covariate used to align (true or false)
:align_disc_cov => 1,							# Which discontinuous covariate is used to choose samples to separately align (is an integer)
:align_other_covariates => false,				# Are other covariates included
:align_batch_only => false,
# :align_samples => sample_ids_with_collection_times,
# :align_phases => sample_collection_times,
# :align_genes => Array{String, 1},				# A string array of genes used to align CYCLOPS fit output. Goes together with :align_acrophases
# :align_acrophases => Array{<: Number, 1}, 	# A number array of acrophases for each gene used to align CYCLOPS fit output. Goes together with :align_genes

:X_Val_k => 10,									# How many folds used in cross validation.
:X_Val_omit_size => 0.1,						# What is the fraction of samples left out per fold

:plot_use_o_cov => true,
:plot_correct_batches => true,
:plot_disc => false,
:plot_disc_cov => 1,
:plot_separate => false,
:plot_color => ["b", "orange", "g", "r", "m", "y", "k"],
:plot_only_color => true,
:plot_p_cutoff => 0.05)

Distributed.addprocs(10)
@everywhere include(joinpath(homedir(), "Desktop", "CYCLOPS-2.0", "CYCLOPS.jl"))

subset_subjects = TMM[:, Vector(TMM[1,:])  .!= "cond_1"]

#Apply the pca from normals to all subjects:
eigen_data, metricDataframe, correlationDataframe, best_model, options = CYCLOPS.TransferFit_d1(subset_subjects, TMM, seed_genes, training_parameters)
file_name, ~ = CYCLOPS.Align(TMM, metricDataframe, correlationDataframe, best_model, options, output_path)

#for reapplying fit for just normal subjects
eigen_data2, metricDataframe2, correlationDataframe2, model2, options2  = CYCLOPS.ReApplyFit_d1(best_model, subset_subjects, TMM, subset_subjects, seed_genes, training_parameters)
CYCLOPS.Align(subset_subjects, metricDataframe2, correlationDataframe2, model2, options2, joinpath(output_path, file_name))

smoothness_norms = CYCLOPS.smoothness_metric(eigen_data2, model2, metricDataframe2; covariates = true)
smoothness_all = CYCLOPS.smoothness_metric(eigen_data, best_model, metricDataframe; covariates = true)
touch(joinpath(output_path,file_name, "smoothness.txt"))
file = open(joinpath(output_path,file_name, "smoothness.txt"), "w")
write(file, "eigendata on all subs, sqrt smoothness: $smoothness_all \neigendata on norm subs, sqrt smoothness: $smoothness_norms ")
close(file)

#Running regular CYCLOPS covariates (not transfer fit) on just the controls. This fit putput should match the statErr_unshuffled_data_best_circular_fit.csv in the next step:
eigen_data_ctl_run, metricDataframe_ctl_run, correlationDataframe_ctl_run, best_model_ctl_run, options_ctl_run = CYCLOPS.Fit(subset_subjects, seed_genes, training_parameters)
CYCLOPS.Align(subset_subjects, metricDataframe_ctl_run, correlationDataframe_ctl_run, best_model_ctl_run, options_ctl_run, output_path)

smoothness = CYCLOPS.smoothness_metric(eigen_data_ctl_run,best_model_ctl_run , metricDataframe_ctl_run; covariates = true)
touch(joinpath(output_path,"smoothness_ctl_model.txt"))
file = open(joinpath(output_path,"smoothness_ctl_model.txt"), "w")
write(file, "eigendata on all subs, sqrt smoothness: $smoothness")
close(file)


#Run stat_err on the control subjects:
stat_err = CYCLOPS.stat_err(subset_subjects, seed_genes, training_parameters, output_path; perms = 200, eigen_shuffle = true, freeze_covs = true, covs = true)

