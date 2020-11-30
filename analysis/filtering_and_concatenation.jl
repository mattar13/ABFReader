### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ acb06ef0-042f-11eb-2b35-e7f2578cf3bd
using PlutoUI

# ╔═╡ eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
using NeuroPhys

# ╔═╡ 7aec9f70-23e6-11eb-293a-8f4e69df4f50
using Plots, DSP, StatsPlots

# ╔═╡ 8c373e20-23dc-11eb-01b3-9143dc22e796
using DataFrames, Query

# ╔═╡ e1d96250-31af-11eb-2719-0ffa95a30d85
using PyCall

# ╔═╡ e7c07a90-042e-11eb-2565-8f992ddf6aea
pyplot()

# ╔═╡ 6aa33000-0426-11eb-3757-d55b61aebc53
md"
## [1] Filtering individual ERG files

- The first thing we can do to generate ERG reports enter the location of the file in here

"

# ╔═╡ e09e64b0-0425-11eb-1a08-8f78d2ceca08
target_path = "to_filter.abf"

# ╔═╡ cc74a240-042c-11eb-257c-f969882fcc79
md"
### Data Cleaning

This notebook will be in depth in order to show the functionality of the filtering
We can clean the data using the functions
- Drift cancelling is just a polynomial fit
- Subtract baseline
- Normalize to region before stimulus
- Continuous wavelet transform filtering
"

# ╔═╡ a3d6e720-31af-11eb-2a47-85b72bb63cf9
begin
	#First open the file
	trace = extract_abf(target_path);
end

# ╔═╡ 76b8eb6e-31c9-11eb-3a02-3b5f312dc1c8
begin
	p1 = plot(layout = grid(3,1), 
		title = ["Unfiltered_data" "" ""], 
		xlabel = ["" "" "Time (s)"], 
		ylabel = ["Response (\$\\mu\$V)" "Response (\$\\mu\$V)" "Stimulus"]
	)
	for (i, ch) in enumerate(NeuroPhys.eachchannel(trace))
		plot!(p1[i], trace.t, ch |> vec, label = "", c = :blue)
	end
	p1
end

# ╔═╡ cdb8fba0-2a98-11eb-022e-bdc0a537368d
md"
### Adjusting baseline and normalization
"

# ╔═╡ f30e86ae-3275-11eb-059c-43cea4921771
begin
	drift_trace = baseline_cancel(trace; mode = :slope, region = :prestim)
	p2 = plot(layout = grid(3,1), 
		title = ["Adjusting Drift" "" ""], 
		xlabel = ["" "" "Time (s)"], 
		ylabel = ["Response (\$\\mu\$V)" "Response (\$\\mu\$V)" "Stimulus"]
	)
	for (i, ch) in enumerate(eachchannel(drift_trace))
		plot!(p2[i], trace.t, getchannel(trace,i), label = "Unfiltered", c = :red)
		plot!(p2[i], drift_trace.t, ch |> vec, label = "Filtered", c = :blue)
	end
	p2
end

# ╔═╡ 1fcf25b0-0431-11eb-0c6c-2d2204083a98
begin
	baseline_trace = baseline_cancel(drift_trace; mode = :mean, region = :prestim)
	p3 = plot(layout = grid(3,1), 
		title = ["Adjusting Baseline" "" ""], 
		xlabel = ["" "" "Time (s)"], 
		ylabel = ["Response (\$\\mu\$V)" "Response (\$\\mu\$V)" "Stimulus"]
	)
	for (i, ch) in enumerate(eachchannel(baseline_trace))
		plot!(p3[i], trace.t, getchannel(drift_trace,i), label = "Unfiltered", c = :red)
		plot!(p3[i], trace.t, ch |> vec, label = "Filtered", c = :blue)
	end
	p3
end

# ╔═╡ c1b43420-327f-11eb-0a0d-c7a01291ec64
begin
	filter_trace = lowpass_filter(baseline_trace)
	p4 = plot(layout = grid(3,1), 
		title = ["Lowpass Filter" "" ""], 
		xlabel = ["" "" "Time (s)"], 
		ylabel = ["Response (\$\\mu\$V)" "Response (\$\\mu\$V)" "Stimulus"]
	)
	for (i, ch) in enumerate(eachchannel(filter_trace))
		plot!(p4[i], trace.t, getchannel(baseline_trace,i), label = "Unfiltered", c = :red)
		plot!(p4[i], trace.t, ch, label = "Filtered", c = :blue)
	end
	p4
end

# ╔═╡ 3a38a790-3286-11eb-1c8c-fff9c46bbe7e
begin
	cwt_trace = cwt_filter(baseline_trace)
	p5 = plot(layout = grid(3,1), 
		title = ["CWT Filter" "" ""], 
		xlabel = ["" "" "Time (s)"], 
		ylabel = ["Response (\$\\mu\$V)" "Response (\$\\mu\$V)" "Stimulus"]
	)
	for (i, ch) in enumerate(eachchannel(cwt_trace))
		plot!(p5[i], trace.t, ch, label = "Filtered", c = :blue)
	end
	p5
end

# ╔═╡ 4aee4550-0431-11eb-2643-29f5e0eb19b5

begin
	#Normalization
	x_norm1, _ = NeuroPhys.normalize(x_adj1);
	x_norm2, _ = NeuroPhys.normalize(x_adj2);
	pnorm = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pnorm[1], t, -x_norm1, label = "Normalized", title = "Normalized",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	plot!(pnorm[2], t, -x_norm2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	plot!(pnorm[1], t, x_adj1, label = "Baseline subtracted", 
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	plot!(pnorm[2], t, x_adj2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	
	plot!(pnorm[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 498f2320-0434-11eb-0cc3-f977a71c5196
begin
	import NeuroPhys.fft_spectrum
	stim_points = findall(x -> x>0.0, x_stim)
	t_start = t[stim_points[1]]-0.5
	t_end = 4.0
	p1FF_1 = plot(t, x_adj1, label = "", c = :blue, 
		xlabel = "Time (s)", ylabel = "Response (mV)", 
		title = "Filtering Ch1"
		) 
	plot!(p1FF_1, t, (x_cwt1./maximum(x_cwt1)), label = "CWT filtered", 
		c = :red, lw = 2.0);
	plot!(p1FF_1, t, x_bp1, label = "BS Filtered", 
		c = :green, lw  = 2.0, xlims = (t_start,t_end), ylims = (-1.0, 0.5))
	
	scatter!(p1FF_1, t[stim_points], [repeat([-1.0], length(stim_points))],  
		marker = :square, markersize = 10.0, c = :black, label = "Light Stim"
	)

	#Spectra Analysis
	freqs1, x_fft1 = fft_spectrum(t, -x_norm1);
	freqs_F1, x_fft_F1 = fft_spectrum(t, -x_cwt1);
	freqs_BS1, x_fft_BS1 = fft_spectrum(t, -x_bp1);
	
	p1FF_2 = plot(freqs1, abs.(x_fft1), label = "Normalized", c = :blue,
		yaxis = :log, xaxis = :log,  
		xlabel = "Frequency (hz)", ylabel = "Power (mv/hz)"
		);
	
	plot!(p1FF_2, freqs_F1, abs.(x_fft_F1), label = "CWT Filtered", c = :red, 
		yaxis = :log, xaxis = :log, );
	plot!(p1FF_2, freqs_BS1, abs.(x_fft_BS1), label = "BS Filtered", 
		yaxis = :log, xaxis = :log, 
		c = :green, alpha = 0.8);

	vspan!(p1FF_2, [59, 61], alpha = 0.5, c = :red, label = "Electrical Noise")
	p1FF = plot(p1FF_1, p1FF_2, layout = grid(2,1), size = (1000, 800));

	p2FF_1 = plot(t, x_adj2, label = "", c = :blue, 
		xlabel = "Time (s)", ylabel = "Response (mV)", 
		title = "Filtering Ch2"
		);
	
	plot!(p2FF_1, t, (x_cwt2./maximum(x_cwt2)), label = "CWT filtered", c = :red, lw = 2.0);
	plot!(p2FF_1, t, x_bp2, label = "BS Filtered", c = :green, lw  = 2.0, 
		xlims = (t_start,t_end), ylims = (-1.0, 0.5))
	
	scatter!(p2FF_1, t[stim_points], [repeat([-1.0], length(stim_points))],  
		marker = :square, markersize = 10.0, c = :black, label = "Light Stim"
	)

	#Spectra Analysis
	freqs2, x_fft2 = fft_spectrum(t, -x_norm2);
	freqs_F2, x_fft_F2 = fft_spectrum(t, -x_cwt2);
	freqs_BS2, x_fft_BS2 = fft_spectrum(t, -x_bp2);
	
	
	p2FF_2 = plot(freqs2, abs.(x_fft2), label = "Normalized", 
		yaxis = :log, xaxis = :log, c = :blue, 
		xlabel = "Frequency (hz)", ylabel = "Power (mv/hz)"
		);
	
	plot!(p2FF_2, freqs_F2, abs.(x_fft_F2), label = "CWT Filtered", 
		yaxis = :log, xaxis = :log, 
		c = :red
		);
	
	plot!(p2FF_2, freqs_BS2, abs.(x_fft_BS2), label = "BS Filtered", 
		yaxis = :log, xaxis = :log, 
		c = :green, alpha = 0.8, legend = :bottomleft
		);

	vline!(p2FF_2, [60], c = :black, label = "Electrical Noise")
	p2FF = plot(p2FF_1, p2FF_2, layout = grid(2,1), size = (1000, 800));

	plot(p1FF, p2FF, layout = grid(1,2))
end

# ╔═╡ 31814662-1e1e-11eb-3f29-5bccaf4079af
md"
##### Summary

I think filtering should be used only in very noisy cases. In most cases, the averaged data should be sufficient to make a good analysis. 

In the case filtering needs to be used: 
- CWT preserved the effective times, but can cause errors in the amplitudes
- Bandpass filtering can preserve the amplitudes, but cause errors in the effective time 

Therefore depending on which metric needs to be applied, use the appropriate filtering method
"

# ╔═╡ 4d825730-1e1b-11eb-3e3a-0b1c0d22971e
md"### [2] Concatenating files
"

# ╔═╡ 7ad594de-1e1b-11eb-28ce-e18d72a90517
cat_path = "to_concatenate"

# ╔═╡ f129e1e0-1e21-11eb-060c-b7c6b7444713
paths = cat_path |> parse_abf |> sort

# ╔═╡ b57790e0-1e24-11eb-0b7a-491baff911d1
md"
ERG traces: $(length(paths))
"

# ╔═╡ 90e76340-23dd-11eb-0cf7-f12c2da517d8
t2, concat_data = concat(paths; 
	filter_func = clean_data, t_cutoff = 0.51, t_eff = 0.1
		);

# ╔═╡ cc8d9ef0-23dd-11eb-01f4-c562c2d21a14
md" 
### [3] File analysis

Since we have these files, we will show a quick file analysis. The intention then being that we will delve deeper into analysis in future files
- Intensity Response

Using the previously derived equation for relating Transferrance, Intensity and Stimulus time to number of photons, we can calculate the approximate number of photons in the current recordings. 

$ Photons = t_{stim} T (γI^2 + βI + α) $
"

# ╔═╡ f4739e80-23db-11eb-1511-59f6bfa6dac6
begin
	#Normalize data
	#ch1_norm, norm_factor1 = normalize(concat_data[:,:,1]);
	#ch2_norm, norm_factor2 = normalize(concat_data[:,:,2]);
	ch1_norm = concat_data[:,:,1];
	ch2_norm = concat_data[:,:,2];
	#Extract the response
	R_ch1 = maximum(-ch1_norm, dims = 2) |> vec
	R_ch2 = maximum(-ch2_norm, dims = 2) |> vec
	
	
	analysis = DataFrame(
		OD = Float64[], 
		LED_Intensity = Float64[], 
		Stim_time = Float64[]
	)
	for path in paths
		info = split(path, "\\")[end] |> filename_extractor
		push!(analysis, (info[1], info[2], info[3]))
	end
	analysis[!, :Transferrance] = analysis[!,:OD] .|> Transferrance
	ivar = analysis[!, [:Transferrance, :LED_Intensity, :Stim_time]] |> Array
	analysis[!, :Photons] = stimulus_model(ivar)
	
	
	analysis[!, :Ch1_R] = R_ch1
	analysis[!, :Ch2_R] = R_ch2
	analysis
end

# ╔═╡ f67381a0-23e5-11eb-2c5e-55836c165487
begin
	fig1a = plot(layout = grid(2,1))
	plot_idxs = collect(1:size(concat_data,1))
	for i in plot_idxs
		plot!(fig1a[1], t2, ch1_norm[i,:], label = "", c = :delta, 
			line_z = log(analysis[i, :Photons]), 
			xlabel = "", ylabel = "Response (\$\\mu\$V)"
		)
		plot!(fig1a[2], t2, ch2_norm[i,:], label = "", 
			c = :delta, line_z = log(analysis[i, :Photons]), 
			xlabel = "Time (s)", ylabel = "Response (\$\\mu\$V)"
		)
		stim_start = findall(x -> x == true, concat_data[i,:,3])[1]
		stim_end = findall(x -> x == true, concat_data[i,:,3])[end]
		vspan!(fig1a[1], [t2[stim_start], t2[stim_end]], c = :gray, label = "")
		vspan!(fig1a[2], [t2[stim_start], t2[stim_end]], c = :gray, label = "")
	end
	fig1b = @df analysis plot(:Photons, :Ch1_R, seriestype = :scatter)
	#@df analysis plot!(fig1b, :Photons, :Ch2_R, seriestype = :scatter)
	
	fig1 = plot(fig1a, fig1b, layout = grid(1,2))	
end

# ╔═╡ Cell order:
# ╠═acb06ef0-042f-11eb-2b35-e7f2578cf3bd
# ╠═eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
# ╠═7aec9f70-23e6-11eb-293a-8f4e69df4f50
# ╠═8c373e20-23dc-11eb-01b3-9143dc22e796
# ╠═e1d96250-31af-11eb-2719-0ffa95a30d85
# ╟─e7c07a90-042e-11eb-2565-8f992ddf6aea
# ╟─6aa33000-0426-11eb-3757-d55b61aebc53
# ╠═e09e64b0-0425-11eb-1a08-8f78d2ceca08
# ╟─cc74a240-042c-11eb-257c-f969882fcc79
# ╟─a3d6e720-31af-11eb-2a47-85b72bb63cf9
# ╠═76b8eb6e-31c9-11eb-3a02-3b5f312dc1c8
# ╟─cdb8fba0-2a98-11eb-022e-bdc0a537368d
# ╟─f30e86ae-3275-11eb-059c-43cea4921771
# ╟─1fcf25b0-0431-11eb-0c6c-2d2204083a98
# ╟─c1b43420-327f-11eb-0a0d-c7a01291ec64
# ╟─3a38a790-3286-11eb-1c8c-fff9c46bbe7e
# ╠═4aee4550-0431-11eb-2643-29f5e0eb19b5
# ╠═498f2320-0434-11eb-0cc3-f977a71c5196
# ╟─31814662-1e1e-11eb-3f29-5bccaf4079af
# ╟─4d825730-1e1b-11eb-3e3a-0b1c0d22971e
# ╟─7ad594de-1e1b-11eb-28ce-e18d72a90517
# ╠═f129e1e0-1e21-11eb-060c-b7c6b7444713
# ╟─b57790e0-1e24-11eb-0b7a-491baff911d1
# ╠═90e76340-23dd-11eb-0cf7-f12c2da517d8
# ╟─cc8d9ef0-23dd-11eb-01f4-c562c2d21a14
# ╟─f4739e80-23db-11eb-1511-59f6bfa6dac6
# ╟─f67381a0-23e5-11eb-2c5e-55836c165487
