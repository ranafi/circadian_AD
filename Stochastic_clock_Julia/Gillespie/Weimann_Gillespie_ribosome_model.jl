
using Gillespie
using Plots # Or Gadfly, or another plotting package
using Random
using Plots
using Peaks # Import the package
using Statistics # For mean and std

# If V = 1 pL = 10⁻¹² L:
# correct = (6.022 x 10^(23) molecules/mole) * (10⁻¹² L) / (10⁹ nM L / mole)
# correct = 6.022 x 10^(23 - 12 - 9) molecules/nM

#----parameters----
correct = 602.2
### This script implements a modified version of this model: https://pmc.ncbi.nlm.nih.gov/articles/PMC1304775/

v1b = 9*correct #Maximal rate of Per2/Cry transcription
k1b = 1*correct #Michaelis constant of Per2/Cry transcription
k1i = 0.56*correct #Inhibition constant of Per2/Cry transcription
c = 0.01*correct #Concentration of constitutive activator
p_hill = 8.0 #Hill coefficient of inhibition of Per2/Cry transcription
k1d = 0.12 #Degradation rate of Per2/Cry mRNA
k2b = 43 #/correct #Complex formation rate of PER2/CRY
q = 2.0  #No. of PER2/CRY complex forming subunits
k2d = 0.05 #Degradation rate of the cytoplasmatic PER2/CRY
k2t = 0.24 #Nuclear import rate of the PER2/CRY complex
k3t = 0.02 #Nuclear export rate of the PER2/CRY complex
k3d = 0.12 #Degradation rate of the nuclear PER2/CRY complex
v4b = 3.6*correct  #Maximal rate of Bmal1 transcription
k4b = 2.16*correct #Michaelis constant of Bmal1 transcription
r = 3.0      #Hill coefficient of activation of Bmal1 transcription
k4d = 0.75 #Degradation rate of Bmal1 mRNA
k5b = 43#0.24 #Translation rate of BMAL1
k5d = 0.06 #Degradation rate of cytoplasmatic BMAL1
k5t = 0.45 #Nuclear import rate of BMAL1
k6t = 0.06 #Nuclear export rate of BMAL1
k6d = 0.12 #Degradation rate of nuclear BMAL1
k6a = 0.09 #Activation rate of nuclear BMAL1
k7a = 0.003 #Deactivation rate of nuclear BMAL1*
k7d = 0.09 #Degradation rate of nuclear BMAL1*
K_bind = 1
K_release = 43 #h-1 from (T_transl) = L / Rate_elong = 500 aa / (21,600 aa/h) ≈ 0.023 hours

#params
p = [v1b, k1b, k1i, c, p_hill, k1d, k2b, q, k2d, k2t, k3t, k3d, v4b, k4b, r, k4d, k5b, k5d, k5t, k6t, k6d, k6a, k7a, k7d, K_bind, K_release]

# --- Stoichiometry Matrix ---
      #R1 R2 R3 R4 R5 R6 R7 R8 R9 10 11 12 13 1415 16 17 18 19 20 21 22 23 24 25 26  27
nu = [ 1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  1  0;  # perCry mRNA
       0  0  1  1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  # PER/CRY protein
       0  0  0  0  0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0;  # nucleur PER/CRY protein
       0  0  0  0  0  0  0  0  0  1 -1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0  1;  # bmal mRNA
       0  0  0  0  0  0  0  0  0  0  0  1  1 -1 -1  0  0  0  0  0  0  0  0  0  0  0  0;  # BMAL1 protein
       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1 -1 -1 -1  0  0  0  0  0  0  0;  # nucleur BMAL1 protein
       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  0  0  0;  # active BMAL1 protein
       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1  1  1;  #free ribosome
       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1  0;  #per bound ribosome
       0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0 -1]' #bmal bound ribosome 

# --- Propensity Function  ---
# Takes only state (u) and parameters (p)
function calculate_propensities(u, p)
    y = u
    R_free, R_bound_per, R_bound_bmal = y[8], y[9], y[10]    
    v1b, k1b, k1i, c, p_hill, k1d, k2b, q, k2d, k2t, k3t, k3d, v4b, k4b, r, k4d, k5b, k5d, k5t, k6t, k6d, k6a, k7a, k7d, K_bind, K_release = p

    a1 = ((v1b*(y[7]+c)) / (k1b*(1+(y[3]/k1i)^p_hill) +(y[7]+c))) #create perCRY mRNA
    a2 = k1d*y[1] #degrade perCRY mRNA

    a3 = k2b * R_bound_per #k2b*y[1]^q #create PER2/CRY protein
    a4 = k3t * y[3] #create PER2/CRY protein
    a5 = k2d*y[2] #degrade PER2/CRY protein
    a6 = k2t * y[2] #degrade PER2/CRY protein

    a7 =  k2t*y[2] #create nuclear PER2/CRY protein
    a8 =  k3t*y[3] #degrade nuclear PER2/CRY protein
    a9 = k3d*y[3] #degrade nuclear PER2/CRY protein

    a10 = ((v4b*y[3]^r)/(k4b^r + y[3]^r)) #create bmal mRNA
    a11 = k4d*y[4] #degrade bmal mRNA

    a12 = k5b * R_bound_bmal #k5b*y[4]  #Create BMAL1 protein
    a13 = k6t*y[6] #create BMAL1 protein
    a14 = k5d*y[5] #degrade BMAL1 protein
    a15 = k5t*y[5] #degrade BMAL1 protein

    a16 = k5t*y[5] #create nuclear BMAL1 protein
    a17 = k7a*y[7] #create nuclear BMAL1 protein
    a18 = k6t*y[6] #degrade nuclear BMAL1 protein
    a19 =  k6d*y[6] #degrade nuclear BMAL1 protein
    a20 =  k6a*y[6] #degrade nuclear BMAL1 protein

    a21 = k6a*y[6] #create active BMAL1 protein
    a22 = k7a*y[7] #degrade active BMAL1 protein
    a23 = k7d*y[7] #degrade active BMAL1 protein

    #Translational Propensities
    a24 = K_bind * y[1] * R_free   #Binding of ribosome to perCry mRNA
    a25 = K_bind * y[4] * R_free   #Binding of ribosome to bmal mRNA

    a26 = K_release * R_bound_per #Release of ribosome from perCry mRNA
    a27 = K_release * R_bound_bmal# Release of ribosome from bmal mRNA
    return [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27]
end



# Assuming you have:
# result_t: Vector of time points from ssa
# result_u: Vector of state vectors from ssa
# states = hcat(result_u...)' # Combine results into a matrix (time x species)
# times = result_t

function find_library_peaks_filtered(times, species_data, t_start; min_prominence=nothing, min_height=nothing, window_size=5)
    original_peak_indices = Int[]
    peak_times = Float64[]

    # Find index where steady state starts
    start_idx = findfirst(t -> t >= t_start, times)
    if isnothing(start_idx) || start_idx > length(times)
         println("Warning: t_start is beyond the simulation time.")
         return original_peak_indices, peak_times
    end
     if start_idx == 1
        # Consider adjusting start_idx to 2 if you need a point before the first potential peak
        println("Info: t_start is at or before the first time point. Considering data from index $start_idx.")
    end

    # Extract data from the steady state onwards
    steady_data = species_data[start_idx:end]
    steady_times = times[start_idx:end]

    # Check if there's enough data for peak finding
     if length(steady_data) < 3 # Need at least y[i-1], y[i], y[i+1] for local maxima
        println("Warning: Not enough data points after t_start for peak finding.")
        return original_peak_indices, peak_times
    end

    # --- Step 1: Find all local maxima indices in the steady-state data ---
    # findmaxima returns a tuple, first element is the indices
    peak_indices_rel = findmaxima(steady_data, window_size) # Get only the indices

    if isempty(peak_indices_rel.indices)
        println("Info: No local maxima found after t_start.")
        return original_peak_indices, peak_times
    end

    # --- Step 2: Calculate Prominences (if filtering by prominence) ---
    # We need prominences to filter. Calculate them for all peaks found.
    # peakproms needs the indices and the data. Returns NamedTuple with indices, proms, etc.
    prom_results = peakproms(peak_indices_rel)

    # --- Step 3: Filter based on criteria ---
    filtered_indices_rel = Int[] # Store indices that pass filters
    for i in 1:length(prom_results.indices)
        idx_rel = prom_results.indices[i] # Index relative to steady_data
        prom = prom_results.proms[i]      # Prominence of this peak
        height = steady_data[idx_rel]     # Height of this peak

        passes_prom = isnothing(min_prominence) || prom >= min_prominence
        passes_height = isnothing(min_height) || height >= min_height

        if passes_prom && passes_height
            push!(filtered_indices_rel, idx_rel)
        end
    end

    if isempty(filtered_indices_rel)
        println("Info: No peaks passed the filtering criteria (min_prominence=$(min_prominence), min_height=$(min_height)).")
        return original_peak_indices, peak_times
    end

    # --- Step 4: Convert filtered relative indices back to original indices ---
    original_peak_indices = filtered_indices_rel .+ (start_idx - 1)
    peak_times = times[original_peak_indices]

    return original_peak_indices, peak_times
end




# --- Assuming previous code ran and you have: ---
# times:         Vector of time points from ssa
# states:        Matrix of states (time x species)
# species_data:  Vector for the species you analyzed
# indices_filtered: Indices of peaks passing initial filters (from find_library_peaks_filtered)
# peak_t_filtered: Times of peaks passing initial filters

function consolidate_peaks_and_period(peak_indices, peak_times, species_data, min_separation; outlier_factor=1.5)
    if length(peak_indices) < 2
        println("Need at least 2 peaks to calculate a period.")
        return Int[], Float64[], Float64[], NaN, NaN # Return empty/NaN if not enough peaks
    end

    # Ensure peaks are sorted by time (should be, but good practice)
    sort_order = sortperm(peak_times)
    sorted_indices = peak_indices[sort_order]
    sorted_times = peak_times[sort_order]

    final_peak_indices = Int[]
    final_peak_times = Float64[]

    current_cluster_indices = [sorted_indices[1]]
    current_cluster_times = [sorted_times[1]]

    for i = 2:length(sorted_times)
        # Check time difference to the *previous* peak in the sorted list
        time_diff = sorted_times[i] - sorted_times[i-1]

        if time_diff < min_separation
            # Add to the current cluster
            push!(current_cluster_indices, sorted_indices[i])
            push!(current_cluster_times, sorted_times[i])
        else
            # End of the previous cluster, process it
            # Find the index within the cluster with the highest value
            cluster_values = species_data[current_cluster_indices]
            max_val_idx_in_cluster = argmax(cluster_values)
            representative_index = current_cluster_indices[max_val_idx_in_cluster]
            representative_time = current_cluster_times[max_val_idx_in_cluster] # Use time corresponding to max value index

            push!(final_peak_indices, representative_index)
            push!(final_peak_times, representative_time)

            # Start a new cluster
            current_cluster_indices = [sorted_indices[i]]
            current_cluster_times = [sorted_times[i]]
        end
    end

    # Process the last cluster
    if !isempty(current_cluster_indices)
        cluster_values = species_data[current_cluster_indices]
        max_val_idx_in_cluster = argmax(cluster_values)
        representative_index = current_cluster_indices[max_val_idx_in_cluster]
        representative_time = current_cluster_times[max_val_idx_in_cluster] # Use time corresponding to max value index

        push!(final_peak_indices, representative_index)
        push!(final_peak_times, representative_time)
    end

    # --- Calculate Period ---
    if length(final_peak_times) < 2
         println("Only found 0 or 1 consolidated peak clusters. Cannot calculate period.")
         return final_peak_indices, final_peak_times, Float64[], NaN, NaN
    end

    periods = diff(final_peak_times)
    
    median_period_est = median(periods) # Use median as a robust estimate

    # Define bounds for acceptable periods
    lower_bound = median_period_est / outlier_factor
    upper_bound = median_period_est * outlier_factor

    # Filter periods
    filtered_periods = periods[lower_bound .<= periods .<= upper_bound]
    outlier_indices = findall(.!(lower_bound .<= periods .<= upper_bound)) # Indices of outliers *within the periods array*

    robust_period_est = mean(filtered_periods)
    std_dev_period = std(filtered_periods)
    

    return final_peak_indices, final_peak_times, periods, robust_period_est, std_dev_period, outlier_indices
end


window_size = 50 # Adjust window size for peak finding
min_prom = 1.0 # Example: Require peak to be 
species_to_analyze = 4
min_peak_sep = 5.0 
ribo_test = [30, 25, 20, 18, 16, 14, 12, 10, 8, 6] # Ribosome counts to test
iterate_lst = repeat(ribo_test, inner=(75,))
t_start_steady = 1000.0
out = zeros(size(iterate_lst, 1), 8) # Initialize as an empty array of type Float64
mult = 3 # Multiplier for initial conditions
Random.seed!(1234) # For reproducibility
for (i, ribo_count) in enumerate(iterate_lst)
  
    # Update initial conditions
    u0 = [ 45*mult,70 * mult,90 *mult,70 *mult,25*mult,30*mult,40*mult, ribo_count, 0, 0]
    println("Ribosome count: $(u0[8])")

    # Run SSA again with updated initial conditions
    tf = 2000.0 # Single float value

    # --- Run SSA ---
    result= ssa(u0, calculate_propensities, nu, p, tf)
    

    times = result.time


    start_idx = findfirst(t -> t >= t_start_steady, times)
    median_bmalmrna = median(result.data[start_idx:end, 4])
    median_bmal = median(result.data[start_idx:end, 5])
    median_permrna = median(result.data[start_idx:end, 1])
    median_fr = median(result.data[start_idx:end, 8])

    min_h = median_bmalmrna

    data_vector = result.data[:, species_to_analyze]
    indices_filtered, peak_t_filtered = find_library_peaks_filtered(times, data_vector, t_start_steady, min_prominence=min_prom, min_height=min_h, window_size=window_size)

    # println("Found $(length(indices_filtered)) peaks using Peaks.jl (filtered) at times: $peak_t_filtered")

    # # --- Plotting Verification ---
    # plot(times, data_vector, label="Species $species_to_analyze", xlabel="Time", ylabel="Count", legend=:outertopright)
    # if !isempty(indices_filtered)
    #     scatter!(peak_t_filtered, data_vector[indices_filtered], label="Detected Peaks (Filtered)", markersize=5, color=:red)
    # end

    

    # 2. Consolidate and calculate period
    final_indices, final_times, periods, mean_p, std_p, outlier_p_idx = consolidate_peaks_and_period(
        indices_filtered,
        peak_t_filtered,
        data_vector, # Need the original data vector here
        min_peak_sep
    )
    first_non_outlier = periods[findfirst(i -> !(i in outlier_p_idx), eachindex(periods))]



    println("Found $(length(final_indices)) consolidated peaks.")
    println("Consolidated peak times: $final_times")
    println("Calculated periods: $periods")
    println("Mean Period: $mean_p +/- $std_p")

    # --- Plotting Verification ---
    # plot(times, data_vector, label="Species $species_to_analyze", xlabel="Time", ylabel="Count", legend=:outertopright, title="Peak Consolidation")

    # if !isempty(final_indices)
    #     scatter!(final_times, data_vector[final_indices], label="Consolidated Peaks", markersize=6, color=:red)
    # end
    out[i, :] = [u0[8], mean_p, std_p, first_non_outlier, median_bmalmrna, median_bmal, median_permrna, median_fr] 

end

scatter(out[:, 1], out[:, 2], label="Mean Period", xlabel="Ribosome Count", ylabel="Mean Period", title="Mean Period vs Ribosome Count")
scatter(out[:, 1], out[:, 3], label="Standard Deviation of Period", xlabel="Ribosome Count", ylabel="Standard Deviation of Period", title="Standard Deviation of Period vs Ribosome Count")


using Plots, StatsPlots, DataFrames
n_per_group = 20
distinct_ribo_counts = unique(out[:, 1]) 


df = DataFrame(out, :auto)
rename!(df,:x1 => :RibosomeCount, :x2 => :MeanPeriod, :x3 => :StdDevPeriod)

plt = @df df boxplot(
    :RibosomeCount,    # Grouping variable for x-axis
    :StdDevPeriod,     # Values for y-axis
    xlabel = "Ribosome Count",
    ylabel = "Standard Deviation of Period",
    title = "Period Std. Dev. by Ribosome Count",
    legend = false,      # The boxplot itself doesn't need a legend entry per box
    linecolor = :black,  # Color of boxplot lines
    fillalpha = 0.75,    # Transparency of box fill
    # To ensure boxes are plotted for each unique RibosomeCount value
    # and not treated as a continuous variable for a single box:
    # group = :RibosomeCount # Often implicit if x is suitable
    xticks = (unique(sort(df.RibosomeCount)), string.(unique(sort(df.RibosomeCount)))) # Ensure all distinct counts are ticks
)

# 2. Overlay the individual data points (scatter plot)
# We use the mutating version scatter! to add to the existing plot 'plt'
@df df scatter!(
    plt, # Add to the plot 'plt'
    :RibosomeCount,
    :StdDevPeriod,
    markercolor = :black,
    markersize = 3,
    markeralpha = 0.5,
    label = "" # No separate legend entry for these points
)

# 3. Calculate means for each group and plot a line connecting them
# Group by RibosomeCount and calculate the median of StdDevPeriod
grouped_meds = combine(groupby(df, :RibosomeCount), :StdDevPeriod => median => :MedStdDevPeriod)

# Sort by RibosomeCount to ensure the line connects points in the correct order
sort!(grouped_meds, :RibosomeCount)

# Plot the line connecting the means
# Using @df again for the new 'grouped_meds' DataFrame
@df grouped_meds plot!(
    plt, # Add to the plot 'plt'
    :RibosomeCount,
    :MedStdDevPeriod,
    seriestype = :line,       # Specify line plot
    linecolor = :red,
    linewidth = 2,
    linestyle = :dash,
    marker = :circle,        # Add markers at each mean point
    markercolor = :red,
    label = "Mean Trend"     # Label for the legend
)

# Display the plot
display(plt)

savefig(plt, joinpath(homedir(), "Desktop", "ribosome_stddev_boxplot.pdf"))




# 1. Create the grouped violin plot
plt = @df df violin( 
    :RibosomeCount,
    :StdDevPeriod,
    xlabel = "Total ribosome count",
    ylabel = "Standard deviation of period (hrs)",
    title = "Period std. dev. by ribosome count",
    legend = false,
    linecolor = :black, # Outline of the violin
    fillalpha = 0.75,
    # side = :both, # Default, can also be :right or :left for half-violins
    # points = false, # Default, we are adding points separately with scatter!
    xticks = (unique(sort(df.RibosomeCount)), string.(unique(sort(df.RibosomeCount))))
)

# 2. Overlay the individual data points (scatter plot) - NO CHANGE HERE
@df df scatter!(
    plt,
    :RibosomeCount,
    :StdDevPeriod,
    markercolor = :black,
    markersize = 3,
    markeralpha = 0.5,
    label = ""
)
savefig(plt, joinpath(homedir(), "Desktop", "ribosome_stddev_violinplot.pdf"))

#####Plots for mean period:

plt = @df df boxplot(
    :RibosomeCount,    # Grouping variable for x-axis
    :MeanPeriod,     # Values for y-axis
    xlabel = "Ribosome Count",
    ylabel = "Mean Period (hrs)",
    title = "Period length by Ribosome Count",
    legend = false,      # The boxplot itself doesn't need a legend entry per box
    linecolor = :black,  # Color of boxplot lines
    fillalpha = 0.75,    # Transparency of box fill
    # To ensure boxes are plotted for each unique RibosomeCount value
    # and not treated as a continuous variable for a single box:
    # group = :RibosomeCount # Often implicit if x is suitable
    xticks = (unique(sort(df.RibosomeCount)), string.(unique(sort(df.RibosomeCount)))) # Ensure all distinct counts are ticks
)

# 2. Overlay the individual data points (scatter plot)
# We use the mutating version scatter! to add to the existing plot 'plt'
@df df scatter!(
    plt, # Add to the plot 'plt'
    :RibosomeCount,
    :MeanPeriod,
    markercolor = :black,
    markersize = 3,
    markeralpha = 0.5,
    label = "" # No separate legend entry for these points
)

# 3. Calculate means for each group and plot a line connecting them
# Group by RibosomeCount and calculate the median of StdDevPeriod
grouped_meds = combine(groupby(df, :RibosomeCount), :MeanPeriod => median => :Median_meanPeriod)

# Sort by RibosomeCount to ensure the line connects points in the correct order
sort!(grouped_meds, :RibosomeCount)

# Plot the line connecting the means
# Using @df again for the new 'grouped_meds' DataFrame
@df grouped_meds plot!(
    plt, # Add to the plot 'plt'
    :RibosomeCount,
    :Median_meanPeriod,
    seriestype = :line,       # Specify line plot
    linecolor = :red,
    linewidth = 2,
    linestyle = :dash,
    marker = :circle,        # Add markers at each mean point
    markercolor = :red
)

# Display the plot
display(plt)

savefig(plt, joinpath(homedir(), "Desktop", "ribosome_period_boxplot.pdf"))




# 1. Create the grouped violin plot
plt = @df df violin( 
    :RibosomeCount,
    :MeanPeriod,
    xlabel = "Total ribosome count",
    ylabel = "Mean period (hrs)",
    title = "Period length by ribosome count",
    legend = false,
    linecolor = :black, # Outline of the violin
    fillalpha = 0.75,
    # side = :both, # Default, can also be :right or :left for half-violins
    # points = false, # Default, we are adding points separately with scatter!
    xticks = (unique(sort(df.RibosomeCount)), string.(unique(sort(df.RibosomeCount))))
)

# 2. Overlay the individual data points (scatter plot) - NO CHANGE HERE
@df df scatter!(
    plt,
    :RibosomeCount,
    :MeanPeriod,
    markercolor = :black,
    markersize = 3,
    markeralpha = 0.5,
    label = ""
)
savefig(plt, joinpath(homedir(), "Desktop", "ribosome_period_violinplot.pdf"))



# --- Simulation End Time ---
tf = 1000.0 # Single float value
u0 = [ 45,70,90,70,25,30,40, 6, 0, 0]
# --- Run SSA ---
Random.seed!(1234) # For reproducibility
result= ssa(u0, calculate_propensities, nu, p, tf)

# --- Plotting ---
# Process the output (list of state vectors) into a matrix for plotting
times = result.time
states = result.data

plot(times, result.data[:,4])