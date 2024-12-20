using PyPlot, Colors

# Check if the correct number of command-line arguments is provided
if length(ARGS) != 7
    println("Usage: julia plot_clock_face.jl <plotname> <Gene_Symbols> <Gene_Phases> <Ref_Gene_Symbols> <Ref_Phases> <BHQ> <amp_ratios>")
    exit(1)  # Exit with error status
end

plotname_arg = ARGS[1]

observed_symbol_names = ARGS[2]
GOI_arg = split(observed_symbol_names, ",")

observed_symbol_phase = ARGS[3]
GOI_Acrophases_arg = map(x -> parse(Float64, x), split(observed_symbol_phase, ","))

ref_symbol_phase = ARGS[4]
GOI_Ideal_arg = split(ref_symbol_phase, ",")

Ideal_Acrophases = ARGS[5]
Ideal_Acrophases_arg = map(x -> parse(Float64, x), split(Ideal_Acrophases, ","))

BHQ = ARGS[6]
signif_Values_arg = map(x -> parse(Float64, x), split(BHQ, ","))

AR = ARGS[7]
amp_ratio_arg = map(x -> parse(Float64, x), split(AR, ","))


function Circular_Mean(phases::Array{T,1} where T <: Union{Float64, Float32})
	if length(phases) < 1
		return NaN
	end
  sinterm=sum(sin.(phases))
  costerm=sum(cos.(phases))
  return mod(atan(sinterm, costerm), 2pi)
end
function ±(X, Y)
    X_plus_Y = X .+ Y
    X_minus_Y = X .- Y
    return X_plus_Y, X_minus_Y
end


function plot_clock_face(plotname, GOI,GOI_Acrophases, GOI_Ideal, Ideal_Acrophases, signif_Values, amp_ratio; space_factor=π/35.5, subplot_space_factor_ratio = 1.5, subplot_fontsize = 8)

    fig = figure(figsize = (10, 11))
    ax = PyPlot.axes(polar = true)
	ax.spines["polar"].set_visible(false)
    axis([0, 2π, 0, 1.1])
    xticks([π/2, π, 3π/2, 2π], [L"$0$", L"$\frac{3π}{2}$", L"$π$", L"$\frac{π}{2}$"], fontsize = 22)
    yticks([0, 0.5, 1], ["", "", ""])
    ax.yaxis.grid(true)
    ax.xaxis.grid(false)
	# println("\n\n")
	# show(stdout, P_Values)
    #comparable_indices = signif_Values .< signif_Cutoff
	# println("\n\n")
	# show(stdout, comparable_indices)
	# println("\n\n")
    # Comparable_GOI = GOI[comparable_indices]
    # Comparable_GOI_Ideal = GOI_Ideal[comparable_indices]
    # comparable_core_clock_acrophases = GOI_Acrophases[comparable_indices]
    # mouse_comparable_acrophases = Ideal_Acrophases[comparable_indices]
    plot_acrophases = mod.(-(GOI_Acrophases .- π/2), 2π)
    plot_ideals = mod.(-(Ideal_Acrophases .- π/2), 2π)
    # cmap = ColorMap("mycmap",[RGB(1,0,0),RGB(.5,.5,.5),RGB(0,1,0)])
    PyPlot.set_cmap("RdYlBu_r")
    minus_log_sig =  -log.(signif_Values)
    sc = PyPlot.scatter(plot_acrophases, ones(length(plot_acrophases)), alpha = 0.8, c =minus_log_sig, s = amp_ratio*1e3, label = "Estimated Acrophases")
    PyPlot.scatter(plot_ideals, ones(length(plot_acrophases)) .* 0.5, alpha = 0.8, s = 75, c = "g", label = "Ideal Acrophases")
    legend( loc = (0.45, -0.1875))
    clim(minimum(minus_log_sig), maximum(minus_log_sig))
    colorbar( orientation = "horizontal", shrink = .8, label = "-log(BHQ value)")

    significant_acrophase_mean = Circular_Mean(plot_acrophases)
    range_upper, range_lower = mod.(±(significant_acrophase_mean, π/2), 2π)
    significant_ideal_acrophase_mean = Circular_Mean(plot_ideals)
    range_upper_ideal, range_lower_ideal = mod.(±(significant_ideal_acrophase_mean, π/2), 2π)
    if range_upper > range_lower
        acrophases_in_range_logical = range_lower .< plot_acrophases .<= range_upper
    else
        acrophases_in_range_logical = .!(range_upper .< plot_acrophases .<= range_lower)
    end
    if range_upper_ideal > range_lower_ideal
        ideal_acrophases_in_range_logical = range_lower_ideal .< plot_ideals .<= range_upper_ideal
    else
        ideal_acrophases_in_range_logical = .!(range_upper_ideal .< plot_ideals .<= range_lower_ideal)
    end
    in_range_sig_acrophases_logical = acrophases_in_range_logical
    in_range_sig_ideal_acrophases_logical = ideal_acrophases_in_range_logical
    acrophase_in_range_mean = Circular_Mean(plot_acrophases[in_range_sig_acrophases_logical])
    ideal_acrophase_in_range_mean = Circular_Mean(plot_ideals[in_range_sig_ideal_acrophases_logical])
    closest_gene_index = findmin(acos.(cos.(acrophase_in_range_mean .- plot_acrophases)))[2]
    closest_ideal_gene_index = findmin(acos.(cos.(ideal_acrophase_in_range_mean .- plot_ideals)))[2]
    middle_gene_acrophase = plot_acrophases[closest_gene_index]
    ideal_middle_gene_acrophase = plot_ideals[closest_ideal_gene_index]
    annotate(GOI[closest_gene_index], xy = [middle_gene_acrophase, 1], xytext = [middle_gene_acrophase, 1.35], arrowprops=Dict("arrowstyle"=>"->", "facecolor"=>"grey"), fontsize = 12)
    annotate(GOI[closest_ideal_gene_index], xy = [ideal_middle_gene_acrophase, 0.5], xytext = [ideal_middle_gene_acrophase, 0.75], arrowprops=Dict("arrowstyle"=>"->", "facecolor"=>"grey"), fontsize = 12)
    distance_from_middle_annotation = middle_gene_acrophase .- plot_acrophases
    sig_phases_larger_middle_logical = distance_from_middle_annotation .< 0
    sig_phases_smaller_middle_logical = distance_from_middle_annotation .> 0
    distance_from_ideal_middle_annotation = ideal_middle_gene_acrophase .- plot_ideals
    sig_phases_larger_ideal_middle_logical = distance_from_ideal_middle_annotation .< 0
    sig_phases_smaller_ideal_middle_logical = distance_from_ideal_middle_annotation .> 0
    sig_larger_phases = plot_acrophases[sig_phases_larger_middle_logical]
    sig_smaller_phases = plot_acrophases[sig_phases_smaller_middle_logical]
    sig_larger_ideal_phases = plot_ideals[sig_phases_larger_ideal_middle_logical]
    sig_smaller_ideal_phases = plot_ideals[sig_phases_smaller_ideal_middle_logical]
    sorted_sig_larger_phases = sort(sig_larger_phases)
    sorted_sig_smaller_phases = sort(sig_smaller_phases, rev = true)
    sorted_sig_larger_ideal_phases = sort(sig_larger_ideal_phases)
    sorted_sig_smaller_ideal_phases = sort(sig_smaller_ideal_phases, rev = true)
    annotation_x_vals_larger = deepcopy(sorted_sig_larger_phases)
    nearest_neighbor = diff(vcat(middle_gene_acrophase, annotation_x_vals_larger))
    nearest_neighbor_too_close_logical = nearest_neighbor .< space_factor
    length(nearest_neighbor_too_close_logical) > 1 ? too_close = |(nearest_neighbor_too_close_logical...) : too_close = false
    while too_close
        annotation_x_vals_larger .+= space_factor .* nearest_neighbor_too_close_logical
        nearest_neighbor = diff(vcat(middle_gene_acrophase, annotation_x_vals_larger))
        nearest_neighbor_too_close_logical = nearest_neighbor .< space_factor
        too_close = |(nearest_neighbor_too_close_logical...)
    end
    annotation_x_vals_smaller = deepcopy(sorted_sig_smaller_phases)
    nearest_neighbor = diff(vcat(middle_gene_acrophase, annotation_x_vals_smaller))
    nearest_neighbor_too_close_logical = nearest_neighbor .> space_factor
    length(nearest_neighbor_too_close_logical) > 1 ? too_close = |(nearest_neighbor_too_close_logical...) : too_close = false
    while too_close
        annotation_x_vals_smaller .-= space_factor .* nearest_neighbor_too_close_logical
        nearest_neighbor = diff(vcat(middle_gene_acrophase, annotation_x_vals_smaller))
        nearest_neighbor_too_close_logical = nearest_neighbor .> -space_factor
        too_close = |(nearest_neighbor_too_close_logical...)
    end
    for mwm in 1:length(sorted_sig_larger_phases)
        desired_phases = sorted_sig_larger_phases[mwm]
        desired_annotation_phases = annotation_x_vals_larger[mwm]
        desired_gene_index = findall(in(desired_phases), plot_acrophases)[1]
        annotate(GOI[desired_gene_index], xy = [desired_phases, 1], xytext = [desired_annotation_phases, 1.35], arrowprops=Dict("arrowstyle"=>"->"), fontsize = 12)
    end
    for wmw in 1:length(sorted_sig_smaller_phases)
        desired_phases = sorted_sig_smaller_phases[wmw]
        desired_annotation_phases = annotation_x_vals_smaller[wmw]
        desired_gene_index = findall(in(desired_phases), plot_acrophases)[1]
        annotate(GOI[desired_gene_index], xy = [desired_phases, 1], xytext = [desired_annotation_phases, 1.35], arrowprops=Dict("arrowstyle"=>"->"), fontsize = 12)
    end
    ideal_annotation_x_vals_larger = deepcopy(sorted_sig_larger_ideal_phases)
    nearest_neighbor = diff(vcat(ideal_middle_gene_acrophase, ideal_annotation_x_vals_larger))
    nearest_neighbor_too_close_logical = nearest_neighbor .< (space_factor * subplot_space_factor_ratio)
    length(nearest_neighbor_too_close_logical) > 1 ? too_close = |(nearest_neighbor_too_close_logical...) : too_close = false
    while too_close
        ideal_annotation_x_vals_larger .+= (space_factor * subplot_space_factor_ratio) .* nearest_neighbor_too_close_logical
        nearest_neighbor = diff(vcat(ideal_middle_gene_acrophase, ideal_annotation_x_vals_larger))
        nearest_neighbor_too_close_logical = nearest_neighbor .< (space_factor * subplot_space_factor_ratio)
        too_close = |(nearest_neighbor_too_close_logical...)
    end
    ideal_annotation_x_vals_smaller = deepcopy(sorted_sig_smaller_ideal_phases)
    nearest_neighbor = diff(vcat(ideal_middle_gene_acrophase, ideal_annotation_x_vals_smaller))
    nearest_neighbor_too_close_logical = nearest_neighbor .> -(space_factor * subplot_space_factor_ratio)
    length(nearest_neighbor_too_close_logical) > 1 ? too_close = |(nearest_neighbor_too_close_logical...) : too_close = false
    while too_close
        ideal_annotation_x_vals_smaller .-= (space_factor * subplot_space_factor_ratio) .* nearest_neighbor_too_close_logical
        nearest_neighbor = diff(vcat(ideal_middle_gene_acrophase, ideal_annotation_x_vals_smaller))
        nearest_neighbor_too_close_logical = nearest_neighbor .> -(space_factor * subplot_space_factor_ratio)
        too_close = |(nearest_neighbor_too_close_logical...)
    end
    c = 1
    while c <= length(sorted_sig_larger_ideal_phases)
        desired_phases = sorted_sig_larger_ideal_phases[c]
        desired_gene_indices = findall(in(desired_phases), plot_ideals)
        sig_desired_gene_indices = desired_gene_indices
        # sig_desired_gene_indices = intersect(desired_gene_indices, findall(l_sig_p_values[kk]))
        d = 1
        while d <= length(sig_desired_gene_indices)
            desired_annotation_phases = ideal_annotation_x_vals_larger[c + d - 1]
            annotate(GOI_Ideal[sig_desired_gene_indices[d]], xy = [desired_phases, 0.5], xytext = [desired_annotation_phases,0.75], arrowprops=Dict("arrowstyle"=>"->"), fontsize = 12)
            d += 1
        end
        c += (d - 1)
    end
    c = 1
    while c <= length(sorted_sig_smaller_ideal_phases)
        desired_phases = sorted_sig_smaller_ideal_phases[c]
        desired_gene_indices = findall(in(desired_phases), plot_ideals)
        sig_desired_gene_indices = desired_gene_indices
        # sig_desired_gene_indices = intersect(desired_gene_indices, findall(l_sig_p_values[kk]))
        d = 1
        while d <= length(sig_desired_gene_indices)
            desired_annotation_phases = ideal_annotation_x_vals_smaller[c + d - 1]
            annotate(GOI_Ideal[sig_desired_gene_indices[d]], xy = [desired_phases, 0.5], xytext = [desired_annotation_phases,0.75], arrowprops=Dict("arrowstyle"=>"->"), fontsize = 12)
            d += 1
        end
        c += (d - 1)
    end
    savefig(plotname * ".png", bbox_inches = "tight", dpi = 300)
    savefig(plotname * ".pdf", bbox_inches = "tight")

    close(fig)
end

plot_clock_face(plotname_arg, GOI_arg, GOI_Acrophases_arg, GOI_Ideal_arg,  Ideal_Acrophases_arg, signif_Values_arg, amp_ratio_arg)
