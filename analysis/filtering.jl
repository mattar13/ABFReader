### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ acb06ef0-042f-11eb-2b35-e7f2578cf3bd
using PlutoUI

# ╔═╡ eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
using NeuroPhys, Plots

# ╔═╡ e7c07a90-042e-11eb-2565-8f992ddf6aea
pyplot()

# ╔═╡ 6aa33000-0426-11eb-3757-d55b61aebc53
md"
## Analyzing Individual ERG files

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

# ╔═╡ 5dfb2940-042e-11eb-1d71-d3d70aed94e4
begin
	#First open the file
	t, raw_data, dt = extract_abf(target_path);
	data = sum(raw_data, dims = 1)/size(raw_data,1);
	x_ch1 = data[1,:,1]; x_ch2 = data[1,:,2]; x_stim = data[1,:,3] .> 0.2;
	p1 = plot(layout = grid(3,1), xlims = (0.0, 10.0),)
	plot!(p1[1], t, x_ch1, label = "", title = "Unfiltered Data",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	plot!(p1[2], t, x_ch2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	plot!(p1[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 8e5be320-0430-11eb-2ea2-c9fbd7e40caa
begin
	#First open the file
	#Cancelling drift
	x_lin1 = NeuroPhys.drift_cancel(t, x_ch1);
	x_lin2 = NeuroPhys.drift_cancel(t, x_ch2);
	pdrift = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pdrift[1], t, x_lin1, label = "Drift Cancelled Data", 
		title = "Drift Cancelled",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	plot!(pdrift[2], t, x_lin2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	#Unfiltered data
	plot!(pdrift[1], t, x_ch1, label = "Unfiltered data",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	plot!(pdrift[2], t, x_ch2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	
	
	plot!(pdrift[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 1fcf25b0-0431-11eb-0c6c-2d2204083a98
begin
	#Baseline subtraction
	stim_idxs = findall(x -> x == true, x_stim) #Stimulus is same for both channels
	x_adj1 = NeuroPhys.subtract_baseline(x_lin1, (1, stim_idxs[1]));
	x_adj2 = NeuroPhys.subtract_baseline(x_lin2, (1, stim_idxs[1]));
	pbase = plot(layout = grid(3,1), xlims = (0.0, 10.0))
	plot!(pbase[1], t, x_adj1, label = "Baseline subtracted", 
		title = "Baseline subtracted",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	plot!(pbase[2], t, x_adj2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :blue
	)
	
	#Drift cancelled data
	plot!(pbase[1], t, x_lin1, label = "Drift Cancelled Data", 
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	plot!(pbase[2], t, x_lin2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red
	)
	
	plot!(pbase[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 4aee4550-0431-11eb-2643-29f5e0eb19b5

begin
	#Normalization
	x_norm1, norm_factor1 = NeuroPhys.normalize(x_adj1);
	x_norm2, norm_factor2 = NeuroPhys.normalize(x_adj2);
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

# ╔═╡ 7dabc5d0-0431-11eb-0ca4-dfbfbc09620d

begin
	#CWT filtering (Probably not ready for CWT filtering )
	x_cwt1, cwt1_raster = NeuroPhys.cwt_filter(x_norm1, periods = 1:9);
	x_cwt2, cwt2_raster = NeuroPhys.cwt_filter(x_norm2, periods = 1:9);
	pcwt = plot(layout = grid(3,1), xlims = (0.0, 10.0))

	#Unfiltered
	plot!(pcwt[1], t, -x_norm1, label = "Normalized",c = :blue,
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[2], t, -x_norm2, label = "", title = "",c = :blue,
		xlabel = "", ylabel = "Response (\$\\mu\$V)"
	)
	plot!(pcwt[1], t, -x_cwt1, label = "CWT Filtered", title = "Using CWT filter",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red,
	)
	plot!(pcwt[2], t, -x_cwt2, label = "", title = "",
		xlabel = "", ylabel = "Response (\$\\mu\$V)", c = :red,
	)
	plot!(pcwt[3], t, raw_data[1,:,3], label = "", title = "",
		xlabel = "Time (s)", ylabel = "Stimulus"
	)
end

# ╔═╡ 498f2320-0434-11eb-0cc3-f977a71c5196
begin
	using DSP
	import NeuroPhys.fft_spectrum
	#Lowpass filtering
	responsetype = Lowpass(25.0; fs = 1/dt); designmethod = Butterworth(8)
	x_bp1 = filt(digitalfilter(responsetype, designmethod), x_norm1);
	#Lowpass filtering
	responsetype = Lowpass(25.0; fs = 1/dt); designmethod = Butterworth(8)
	x_bp2 = filt(digitalfilter(responsetype, designmethod), x_norm2);
	
	p1FF_1 = plot(t, -x_norm1, label = "", c = :blue, 
		xlabel = "Time (s)", ylabel = "Response (mV)", 
		title = "Filtering Ch1"
		) 
	plot!(p1FF_1, t, -x_cwt1, label = "CWT filtered", 
		c = :red, lw = 2.0);
	plot!(p1FF_1, t, -x_bp1, label = "BS Filtered", 
		c = :green, lw  = 2.0, xlims = (1,3.0), ylims = (-1.0, 0.5))
	stim_points = findall(x -> x>0.0, x_stim)
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

	p2FF_1 = plot(t, -x_norm2, label = "", c = :blue, 
		xlabel = "Time (s)", ylabel = "Response (mV)", 
		title = "Filtering Ch2"
		);
	
	plot!(p2FF_1, t, -x_cwt2, label = "CWT filtered", c = :red, lw = 2.0);
	plot!(p2FF_1, t, -x_bp2, label = "BS Filtered", c = :green, lw  = 2.0, 
		xlims = (1,3.0), ylims = (-1.0, 0.5))
	
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

# ╔═╡ Cell order:
# ╠═acb06ef0-042f-11eb-2b35-e7f2578cf3bd
# ╠═eec4b7f2-0426-11eb-1f69-b3fea7ffedb1
# ╟─e7c07a90-042e-11eb-2565-8f992ddf6aea
# ╟─6aa33000-0426-11eb-3757-d55b61aebc53
# ╠═e09e64b0-0425-11eb-1a08-8f78d2ceca08
# ╟─cc74a240-042c-11eb-257c-f969882fcc79
# ╠═5dfb2940-042e-11eb-1d71-d3d70aed94e4
# ╟─8e5be320-0430-11eb-2ea2-c9fbd7e40caa
# ╟─1fcf25b0-0431-11eb-0c6c-2d2204083a98
# ╟─4aee4550-0431-11eb-2643-29f5e0eb19b5
# ╟─7dabc5d0-0431-11eb-0ca4-dfbfbc09620d
# ╟─498f2320-0434-11eb-0cc3-f977a71c5196
