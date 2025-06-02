using DifferentialEquations
using Plots
using Statistics
using StatsBase
using BenchmarkTools
using Random
using GLM
using DataFrames
using CSV
### This script implements a modified version of this model: https://pmc.ncbi.nlm.nih.gov/articles/PMC1304775/

v1b = 9 #Maximal rate of Per2/Cry transcription
k1b = 1 #Michaelis constant of Per2/Cry transcription
k1i = 0.56 #Inhibition constant of Per2/Cry transcription
c = 0.01 #Concentration of constitutive activator
p_hill = 8.0 #Hill coefficient of inhibition of Per2/Cry transcription
k1d = 0.12 #Degradation rate of Per2/Cry mRNA
k2b = 0.3 #Complex formation rate of PER2/CRY
q = 2.0  #No. of PER2/CRY complex forming subunits
k2d = 0.05 #Degradation rate of the cytoplasmatic PER2/CRY
k2t = 0.24 #Nuclear import rate of the PER2/CRY complex
k3t = 0.02 #Nuclear export rate of the PER2/CRY complex
k3d = 0.12 #Degradation rate of the nuclear PER2/CRY complex
v4b = 3.6  #Maximal rate of Bmal1 transcription
k4b = 2.16 #Michaelis constant of Bmal1 transcription
r = 3.0      #Hill coefficient of activation of Bmal1 transcription
k4d = 0.75 #Degradation rate of Bmal1 mRNA
k5b = 0.24 #Translation rate of BMAL1
k5d = 0.06 #Degradation rate of cytoplasmatic BMAL1
k5t = 0.45 #Nuclear import rate of BMAL1
k6t = 0.06 #Nuclear export rate of BMAL1
k6d = 0.12 #Degradation rate of nuclear BMAL1
k6a = 0.09 #Activation rate of nuclear BMAL1
k7a = 0.003 #Deactivation rate of nuclear BMAL1*
k7d = 0.09 #Degradation rate of nuclear BMAL1*
#ribo_slow=1
#ribo_noise=0.15


function stock_oscillator!(dy, y, p, t)
    ribo_slow = p[1]
    ribo_noise = p[2]

    dy[1] = ((v1b*(y[7]+c)) / (k1b*(1+(y[3]/k1i)^p_hill) + (y[7]+c))) - k1d*y[1]
    dy[2] = ribo_slow * k2b*y[1]^q - k2d*y[2] - k2t * y[2] + k3t * y[3] 
    dy[3] =  k2t*y[2] - k3t*y[3]-k3d*y[3]
    dy[4] = ((v4b*y[3]^r)/(k4b^r + y[3]^r)) -k4d*y[4]
    dy[5] = ribo_slow * k5b*y[4] - k5d*y[5] - k5t*y[5]+ k6t*y[6]
    dy[6] = k5t*y[5] - k6t*y[6] - k6d*y[6] + k7a*y[7] - k6a*y[6]
    dy[7] = k6a*y[6] - k7a*y[7] - k7d*y[7]
    dy[8] = 0 # reporter
    
end


function σ_modified_oscillator!(dy, y, p, t)
    ribo_slow=p[1]
    ribo_noise=p[2]

    dy[1] = 0
    dy[2] =  ribo_noise * k2b * y[1]
    dy[3] = 0
    dy[4] = 0
    dy[5] = ribo_noise * k5b * y[4]
    dy[6] = 0
    dy[7] = 0
    dy[8] = 0 #reporter
end

############ Verify Model Works Correctly ##########
# init_conds = [ 0.90,1.41,1.76,1.37,0.51,0.65,0.88, 1] #Already on the stable limit cycle

# tspan = (0.0, 1000.0)
# prob = ODEProblem(stock_oscillator!, init_conds, tspan, [1., 0.])
# sol = solve(prob, Tsit5());
# plot(sol, tspan = (400, 472))

#################################

function bmal_extremum(y, t, integrator) 
    ((v4b*y[3]^r)/(k4b^r + y[3]^r)) -k4d*y[4]   #bmal mRNA
end

function counter!(integrator)
    integrator.u[8] = integrator.u[8]*-1
end


cb = ContinuousCallback(bmal_extremum, counter!)
#############################
#### this function converts the data structure that is output from the
#### solution to the differential equation into a matrix 
#############################
function convert_solution_to_matrix(local_sol)
    y1=map(x->x,local_sol.t)
    y2=map(x->x[1],local_sol.u)
    y3=map(x->x[2],local_sol.u)
    y4=map(x->x[3],local_sol.u)
    y5=map(x->x[4],local_sol.u)
    y6=map(x->x[5],local_sol.u)
    y7=map(x->x[6],local_sol.u)
    y8=map(x->x[7],local_sol.u)
    y9=map(x->x[8],local_sol.u)
    y=hcat(y1,y2,y3,y4,y5, y6, y7, y8, y9)
    return y
end

#############################
#### this function takes the matrix form of the solution
#### then finds the times where the reporter variable shifted from -1 to 1 (these are relative maxima)
#### as there is stochastic nosie i worried (and have seen) that there can be a few near points that are noisey local maximaa
#### thus this function will discount transitions less than an hour long
#############################
function find_robust_extrema(solmat)
    max_times=[]
    times=solmat[:,1]
    reporter=solmat[:,9]
    N_entries=length(reporter)

    for index in (1:N_entries-1)
        if ((times[index]==times[index+1]) && (reporter[index+1]==1.0))
            push!(max_times,times[index])
        end
    end

    N_max_times=length(max_times)
    robust_max_times=[max_times[1]]
    for index in (2:N_max_times)
        if ((max_times[index]-robust_max_times[end])>12.0)
            push!(robust_max_times,max_times[index])
        end
    end
    
    return robust_max_times

end

function find_closest_index(time::Vector{Float64}, value::Float64)
    index = argmin(abs.(time .- value))
    return index
end


##############################
function calculate_IV(time::Vector{Float64}, signal::Vector{Float64})
    # Ensure time and signal vectors are the same length
    if length(time) != length(signal)
        error("The 'time' and 'signal' vectors must have the same length.")
    end
    
    # Calculate the differences between consecutive signal values
    signal_diffs = diff(signal)
    
    # Compute the standard deviation of signal differences
    std_diffs = std(signal_diffs)
    
    # Compute the mean of the signal
    mean_signal = mean(signal)
    
    # Handle the case where the mean signal is zero
    if mean_signal == 0
        error("The mean of the signal is zero; IV cannot be computed.")
    end
    
    # Calculate IV
    IV = std_diffs / mean_signal
    return IV
end

function run_model_withparameters(ribo_slow,ribo_noise)
           y0 = [ 0.90,1.41,1.76,1.37,0.51,0.65,0.88, 1]
           prob_sde_weimann = SDEProblem(stock_oscillator!, σ_modified_oscillator!,y0 ,(0.0, 6000.0),[ribo_slow,ribo_noise])
           sol = solve(prob_sde_weimann,SRIW1(),callback = cb,abstol=reltol=1e-6)
        #    plot(sol)
           solution_matrix=convert_solution_to_matrix(sol)
           extremes=find_robust_extrema(solution_matrix)
           cycle_lengths=diff(extremes)
           cycle_lengths=cycle_lengths[20:end]  ###throw out first 20 as coming to quasi steady state
           start_idx = find_closest_index(solution_matrix[:, 1], extremes[20]) #What index does the 20th cycle start at?
           end_idx = find_closest_index(solution_matrix[:, 1], extremes[21])
           IV = calculate_IV(solution_matrix[start_idx:end,1], solution_matrix[start_idx:end, 5])
           return [std(cycle_lengths),mad(cycle_lengths),mean(cycle_lengths),median(cycle_lengths), IV, sol, solution_matrix, start_idx, end_idx]
end

function process_gene_traces(vect)
    ribo_param = map(x -> x[1], vect)
    print("Riboslow param: ", unique(ribo_param))
    try
        if length(unique(ribo_param)) > 1
            throw(ErrorException("The vector contains more than one unique element: $(unique(ribo_param))"))
        end
    catch e
        println("Error: ", e)
        return NaN, NaN, NaN
    end
    bmal_expression = map(x -> x[5], vect)
    per_expression = map(x -> x[6], vect)
    times = map(x -> x[7], vect)
    start_idx = map(x -> x[3], vect)
    end_idx = map(x -> x[4], vect)
    trimmed_bmal= [v[start_idx[i]:end_idx[i]] for (i, v) in enumerate(bmal_expression)]
    trimmed_per = [v[start_idx[i]:end_idx[i]] for (i, v) in enumerate(per_expression)]
    trimmed_times = [v[start_idx[i]:end_idx[i]] .- v[start_idx[i]] for (i, v) in enumerate(times)]

    return trimmed_times, trimmed_bmal, trimmed_per
end

function plot_processed_gene_traces(t1, t2, b1,b2, p1, p2)
    # Plot each trimmed vector
    P = plot(xlims = (0,24), ylims = (0,2), legend =true)
    for i in 1:20
        plot!(P, t1[i], b1[i], c = "green", label="")
        plot!(P, t1[i], p1[i], c = "orange", label="")
        plot!(P, t2[i], b2[i],  c = "green", linestyle =:dash, label="")
        plot!(P, t2[i], p2[i], c = "orange", linestyle =:dash, label="")
        
    end
    # Show the plot
    display(P)
end

function anova_discrete_exp_old(t1, t2, b1,b2, p1, p2)
    times = [0., 4., 8., 12.,16.,20.]
    idx_vec = map(y-> map(x -> find_closest_index(x, y), t1), times)
    # Initialize an empty matrix of zeros
    bmal_expression = zeros(Float64, 20, length(idx_vec))
    per_expression = zeros(Float64, 20, length(idx_vec))
    for (i, indices) in enumerate(idx_vec)
        sampled_bmal= [b1[j][v] for (j, v) in enumerate(indices)]
        sampled_per= [p1[j][v] for (j, v) in enumerate(indices)]
        bmal_expression[:, i] = sampled_bmal
        per_expression[:, i] = sampled_per
    end

    idx_vec2 = map(y-> map(x -> find_closest_index(x, y), t2), times)
    # Initialize an empty matrix of zeros
    bmal_expression2 = zeros(Float64, 20, length(idx_vec))
    per_expression2 = zeros(Float64, 20, length(idx_vec))
    for (i, indices) in enumerate(idx_vec2)
        sampled_bmal2= [b2[j][v] for (j, v) in enumerate(indices)]
        sampled_per2= [p2[j][v] for (j, v) in enumerate(indices)]
        bmal_expression2[:, i] = sampled_bmal2
        per_expression2[:, i] = sampled_per2
    end

  # Convert matrices to long format
    function reshape_expression_matrix(expression_matrix, group_name)
        n_individuals, n_times = size(expression_matrix)
        DataFrame(
            group = fill(group_name, n_individuals * n_times),  # Group identifier
            individual = repeat(1:n_individuals, outer=n_times),  # Individual IDs
            time = repeat(1:n_times, inner=n_individuals),  # Time points
            value = vec(expression_matrix)  # Flatten matrix to column
        )
    end
    # Reshape each group's expression matrix
    per_group1 = reshape_expression_matrix(per_expression, "Group1")
    per_group2 = reshape_expression_matrix(per_expression2, "Group2")
    # Combine both groups into a single DataFrame
    combined_df1 = vcat(per_group1, per_group2)
    CSV.write("/Users/henryhollis/Desktop/per_differences_old_method.csv", combined_df1)

    bmal_group1 = reshape_expression_matrix(bmal_expression, "Group1")
    bmal_group2 = reshape_expression_matrix(bmal_expression2, "Group2")
    # Combine both groups into a single DataFrame
    combined_df2 = vcat(bmal_group1, bmal_group2)
    CSV.write("/Users/henryhollis/Desktop/bmal_differences_old_method.csv", combined_df2)

    # Plot each row of the expression matrix
    P = plot(xlims = (0,20), ylims = (0,2),
        xtickfontsize=12,ytickfontsize=12,legendfontsize=15,
        legend =false)
    for i in 1:size(bmal_expression, 1)
        plot!(P, times, bmal_expression[i, :], marker=:circle, c = "green")
        plot!(P, times, per_expression[i, :], marker=:circle, c = "orange")
        plot!(P, times, bmal_expression2[i, :], marker=:utriangle, c = "green")
        plot!(P, times, per_expression2[i, :], marker=:utriangle, c = "orange")
    end
    display(P)
end

function anova_discrete_exp(t1, t2, b1,b2, p1, p2)
    #TODO can you look at half interval, quarter interval, etc.?
    function helper(time_vec)
        dur = time_vec[end] #total time of oscillation
        first = dur/6
        first_ind =  find_closest_index(time_vec, first)
        second = dur/3
        second_ind =  find_closest_index(time_vec, second)
        third = dur/2
        third_ind =  find_closest_index(time_vec, third)
        fourth = dur/(6/4)
        fourth_ind = find_closest_index(time_vec, fourth)
        fifth = dur/(6/5)
        fifth_ind = find_closest_index(time_vec, fifth)
        sixth = dur[end]
        sixth_ind = find_closest_index(time_vec, sixth)
        return first_ind, second_ind, third_ind, fourth_ind, fifth_ind, sixth_ind
    end

    idx_vec = map(x -> helper(x), t1)
    # Initialize an empty matrix of zeros
    bmal_expression = zeros(Float64, 20, 6)
    per_expression = zeros(Float64, 20, 6)
    for (i, indices) in enumerate(idx_vec)
        sampled_bmal= [b1[i][v] for v in indices]
        sampled_per= [p1[i][v] for v in indices]
        bmal_expression[i, :] = sampled_bmal
        per_expression[i, :] = sampled_per
    end

    idx_vec2 = map(x -> helper(x), t2)
    # Initialize an empty matrix of zeros
    bmal_expression2 = zeros(Float64, 20, 6)
    per_expression2 = zeros(Float64, 20, 6)
    for (i, indices) in enumerate(idx_vec2)
        sampled_bmal= [b2[i][v] for v in indices]
        sampled_per= [p2[i][v] for v in indices]
        bmal_expression2[i, :] = sampled_bmal
        per_expression2[i, :] = sampled_per
    end

  # Convert matrices to long format
    function reshape_expression_matrix(expression_matrix, group_name)
        n_individuals, n_times = size(expression_matrix)
        DataFrame(
            group = fill(group_name, n_individuals * n_times),  # Group identifier
            individual = repeat(1:n_individuals, outer=n_times),  # Individual IDs
            time = repeat(1:n_times, inner=n_individuals),  # Time points
            value = vec(expression_matrix)  # Flatten matrix to column
        )
    end
    # Reshape each group's expression matrix
    per_group1 = reshape_expression_matrix(per_expression, "Group1")
    per_group2 = reshape_expression_matrix(per_expression2, "Group2")
    # Combine both groups into a single DataFrame
    combined_df1 = vcat(per_group1, per_group2)
    CSV.write("/Users/henryhollis/Desktop/per_differences.csv", combined_df1)

    bmal_group1 = reshape_expression_matrix(bmal_expression, "Group1")
    bmal_group2 = reshape_expression_matrix(bmal_expression2, "Group2")
    # Combine both groups into a single DataFrame
    combined_df2 = vcat(bmal_group1, bmal_group2)
    CSV.write("/Users/henryhollis/Desktop/bmal_differences.csv", combined_df2)

    # Plot each row of the expression matrix
    P = plot(xlims = (0,1), ylims = (0,2),
        xtickfontsize=12,ytickfontsize=12,legendfontsize=15,
        legend =false)
    times = [1/6, 2/6, 3/6, 4/6, 5/6, 1]
    for i in 1:size(bmal_expression, 1)
        plot!(P, times, bmal_expression[i, :], marker=:circle, c = "green")
        plot!(P, times, per_expression[i, :], marker=:circle, c = "orange")
        plot!(P, times, bmal_expression2[i, :], marker=:utriangle, c = "green")
        plot!(P, times, per_expression2[i, :], marker=:utriangle, c = "orange")
    end
    display(P)
end


Random.seed!(1234)

Answers2=[1, 0,0,0,24,24, 0]
for i in LinRange(1, .7, 10)
    for j in LinRange(0.03, 0.06, 10)
        for cnt in 1:20
           ribo_noise = j
           ribo_slow =  i
           lans=run_model_withparameters(ribo_slow,ribo_noise)
           sd=lans[1]
           md=lans[2]
           mn=lans[3]
           mdn=lans[4]
           IV=lans[5]
           print([cnt,ribo_slow, ribo_noise])
           Answers2=hcat(Answers2,[ribo_slow, ribo_noise,sd,md,mn,mdn, IV])
       end
    end
end 
Answers2=transpose(Answers2)
Answers2=Answers2[2:end,:]


using Plots; pyplot();

# Group by (ribo_slow, ribo_noise) and average sd
xy = [(Answers2[i,1], Answers2[i, 2]) for i in 1:size(Answers2, 1)]
grouped = Dict()
for i in 1:length(xy)
    push!(get!(grouped, xy[i], []), Answers2[i, 3])
end


# Compute mean z for each group
z_avg = Dict(k => median(v) for (k, v) in grouped)

# Extract averaged data for plotting
x_avg = [k[1] for k in keys(z_avg)]
y_avg = [k[2] for k in keys(z_avg)]
z_avg_values = [v for v in values(z_avg)]


# plot(x_avg, y_avg, z_avg_values,
#     st=:surface,camera=(-30,30),
#     xtickfontsize=12,ytickfontsize=12,ztickfontsize=12,legendfontsize=15,
#     xflip = true)

# 1. Extract and sort the grid points
x_vals = sort(unique(first.(keys(z_avg))))
y_vals = sort(unique(last.(keys(z_avg))))

# 2. Build the Z matrix
#    Note the order: rows correspond to y, columns to x
Z = [ z_avg[(x,y)] for y in y_vals, x in x_vals ]

# 3. Plot the heatmap
p = heatmap(
    x_vals, y_vals, Z;
    xlabel         = "Translation Reduction",
    ylabel         = "Translation Noise",
    colorbar_title = " ",
    xtickfontsize  = 12,
    ytickfontsize  = 12,
    legendfontsize = 15,
    xflip          = true        # if you still want to flip X
)
savefig(p, "/Users/henryhollis/Desktop/heatmap.pdf")

out = DataFrame(Answers2, :auto)
rename!(out, ["ribo_slow", "ribo_noise", "sd", "md", "mn", "mdn", "IV"])  # Assign custom column names
CSV.write("/Users/henryhollis/Desktop/SD_vs_noise_vs_slow.csv", out)

 # Threshold for comparison
 threshold = 1e-6


# Fix ribo_slow at 1 and grab all ribo_noises
fixed_slow = Answers2[abs.(Answers2[:, 1] .- 1.0) .< threshold, :]
scatter(fixed_slow[:,2], fixed_slow[:, 3],
 xlabel="Translational noise",
 ylabel="Standard dev. oscillator period (hr)",
 title="Spread of oscillator period vs translation noise",
 legend=false,
 ylims=(0.1, .4),
 xtickfontsize=12,
 ytickfontsize=12,
 c="#0091ff",
 markersize=10,
 xflip = false)


 # Fix ribo_noise at 0.03 and grab all ribo_slows
fixed_noise = Answers2[abs.(Answers2[:, 2] .- 0.03) .< threshold, :]
scatter(fixed_noise[:,1], fixed_noise[:, 3],
 xlabel="Translational reduction",
 ylabel="Standard dev. oscillator period (hr)",
 title="Spread of oscillator period vs translation slow",
 legend=false,
 ylims=(0.1, .25),
 xtickfontsize=12,
 ytickfontsize=12,
 c="#0091ff",
 markersize=10,
 xflip = true)
 

Answers=[1,0,0,24,24, 0]
gene_tracings = Vector{Tuple{Float64, Int64,Int64, Int64, Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
for i in LinRange(0.03, 0.06, 10)
        for cnt in 1:20
           ribo_noise = i
           ribo_slow =  1 #0.85
           lans=run_model_withparameters(ribo_slow,ribo_noise)
           sd=lans[1]
           md=lans[2]
           mn=lans[3]
           mdn=lans[4]
           IV=lans[5]
           start_idx = lans[8]
           end_idx = lans[9]
           push!(gene_tracings, (i, cnt, start_idx, end_idx, lans[7][:,5], lans[7][:,2], lans[7][:,1]))
           print([cnt,ribo_slow])
           Answers=hcat(Answers,[ribo_slow,sd,md,mn,mdn, IV])
       end
       end
Answers=transpose(Answers)
Answers=Answers[2:end,:]

t1, b1, p1 = process_gene_traces(gene_tracings[1:20])
t2, b2, p2 = process_gene_traces(gene_tracings[181:200])
plot_processed_gene_traces(t1, t2, b1,b2, p1, p2)
# anova_discrete_exp_old(t1, t2, b1,b2, p1, p2)
anova_discrete_exp(t1, t2, b1,b2, p1, p2)


#What if I vary both ribo_slow and ribo_noise?

Random.seed!(1234)
Answers=[1,0,0,24,24, 0]
ribo_noise_list = [0.03, 0.06]
ribo_slow_list = [1.0, 0.66]
gene_tracings = Vector{Tuple{Float64, Int64,Int64, Int64, Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
for i in [1,2]
        for cnt in 1:20
           ribo_noise = ribo_noise_list[i]
           ribo_slow =  ribo_slow_list[i]
           lans=run_model_withparameters(ribo_slow,ribo_noise)
           sd=lans[1]
           md=lans[2]
           mn=lans[3]
           mdn=lans[4]
           IV=lans[5]
           start_idx = lans[8]
           end_idx = lans[9]
           push!(gene_tracings, (i, cnt, start_idx, end_idx, lans[7][:,5], lans[7][:,2], lans[7][:,1]))
           print([cnt,ribo_slow])
           Answers=hcat(Answers,[ribo_slow,sd,md,mn,mdn, IV])
       end
       end
Answers=transpose(Answers)
Answers=Answers[2:end,:]


t1, b1, p1 = process_gene_traces(gene_tracings[1:20])
t2, b2, p2 = process_gene_traces(gene_tracings[21:40])
plot_processed_gene_traces(t1, t2, b1,b2, p1, p2)
anova_discrete_exp(t1, t2, b1,b2, p1, p2)
