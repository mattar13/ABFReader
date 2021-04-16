### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ eb956370-9ba4-11eb-002e-bd530ec32c36
using Revise

# ╔═╡ e3e3b7e5-901b-41ff-88e2-d296891bdaa2
using NeuroPhys

# ╔═╡ de4ec361-908f-4365-8fb2-ca2c3f08b930
using DataFrames, Query

# ╔═╡ 70182ef7-d60a-4073-b55b-7565e6c4ab6c
md"
#### This notebook will analyze the data root that you enter in here
"

# ╔═╡ f57ebcff-72d2-4c0d-9c39-4caaa967b4ef
analysis_file = "E:\\Data\\ERG\\Gnat\\Matt\\2020_11_16_ERG\\Mouse2_P14_KO\\NoDrugs\\525Green"

# ╔═╡ 1e0cfd94-a20b-4691-9b5e-58c5b9e38ff3
begin
	#extract the file
	data = extract_abf(analysis_file; swps = -1)
	t_pre = 0.2;
	t_post = 1.0;
	truncate_data!(data; t_pre = t_pre, t_post = t_post);
	baseline_cancel!(data, mode = :slope, region = :whole); 
	baseline_cancel!(data); #Mean mode
	
	filter_data = data#lowpass_filter(data); #Lowpass filter using a 40hz 8-pole 
	data
end

# ╔═╡ 5b273f51-e6e1-4b2f-9789-6c42c290a8ef
md"
#### Rmaxes, Rdims, Time to peak, and Integration time Printouts

"

# ╔═╡ b5863ba7-f332-43e1-a73f-20a1f6203201
begin
	#Conduct the data anlysis
	rmax_lin = [0.10, 0.30]
	rmaxes = saturated_response(filter_data)
	rdims, dim_idx = dim_response(filter_data, rmaxes; rmax_lin = rmax_lin)
	t_peak = time_to_peak(data, dim_idx)
	t_Int = integration_time(filter_data, dim_idx)
	tau_fit, tau_GOF = recovery_tau(filter_data, dim_idx)
	tau_rec = map(x -> x[2]*1000, tau_fit)
	amp, amp_gof = amplification(data, rmaxes)
	df_analysis = DataFrame(
			:Rmaxes => rmaxes.*-1000,
			:Rdims => rdims.*-1000, 
			:Time_to_peak => t_peak.*1000,
			:Integration_time => t_Int, 
			:TauRec => tau_rec
		)
end

# ╔═╡ 0f308136-d3a1-46d3-95a1-a99c6913b524
md"
#### Amplification, and Intensity Response curves
"

# ╔═╡ 41126fe9-6a3f-495f-9e03-32d980a91bb1
#IR analysis
begin
	minima = minimum(filter_data, dims = 2)[:,1,:]
	IR_analysis = DataFrame(
		:Trace => collect(1:size(data,1)),
		:Photons => 0.0, 
		
	)
	
	for (idx, fn) in enumerate(data.filename)
		nt = formatted_split(fn, format_bank)
		stim_times = data.stim_protocol[idx].timestamps
		t_stim = (stim_times[2] - stim_times[1])*1000
		IR_analysis[idx, :Photons] = f_I(
			nt[:ND]|>Float64, nt[:Intensity]|>Float64, t_stim
		)

		
	end
	for ch in 1:size(data,3)
		IR_analysis[!, Symbol("Minima_$(ch)")] = minima[:, ch] *-1000
		IR_analysis[!, Symbol("Rmax_$(ch)")] = repeat(
			[rmaxes[ch]*-1000], size(data,1)
		)
		IR_analysis[!, Symbol("Amp_$(ch)")] = amp[1, :, ch]
		IR_analysis[!, Symbol("tEff_$(ch)")] = amp[2, :, ch]
		IR_analysis[!, Symbol("R2_$(ch)")] = amp_gof[:, ch]
	end	
	IR_analysis = IR_analysis |> @orderby(_.Photons) |> DataFrame
end

# ╔═╡ f4eff6dc-1d23-45a3-8514-521048ab5abc
begin
	plt = plot(filter_data, 
		c = :black, label_stim = true, grid = false, 
		title = "Rmax, Rdim, and Time to peak")
	#Plot the saturated traces
	
	saturated_traces = findall(minima .< rmaxes')
	for I in saturated_traces
		swp = I[1]
		ch = I[2]
		plot!(plt[ch], filter_data, to_plot = (swp, ch), label ="",
			c = :green, lw = 1.0
		)
	end
	for ch in size(data,3)
		plot!(plt, filter_data, to_plot = (dim_idx[ch], ch), label = "Dim trace",
			c = :red, lw = 2.0
		)
		
		hline!(plt[ch], [rmaxes[ch]], c = :green,label = "Saturation")
		vline!(plt[ch], [t_peak[ch]], c = :magenta,lw= 2.0, label = "Time to peak")
		tR_model(x,p) = map(t -> REC(t, p[1], -1.0), x)
		xdata = data.t
		ydata = data[dim_idx[ch], :, ch] 
		norm_val = minimum(ydata)
		ydata ./= norm_val #Normalize the Rdim
		#cutoff all points below -0.5 and above -1.0
		begin_rng = findall(ydata .>= 1.0)[end]
		xdata = xdata[begin_rng:end]
		ydata = ydata[begin_rng:end]
		end_rng = findall(ydata .< 0.5)[1] 

		xdata = xdata[1:end_rng]
		ydata = -ydata[1:end_rng]
		p0 = [ydata[1], -1.0]
		fit = curve_fit(tR_model, xdata.-xdata[1], ydata, p0)
		plot!(plt[ch], 
			x -> tR_model(x-xdata[1], fit.param)*-norm_val,
			LinRange(xdata[1], xdata[end], 10000), 
			label = "TauRec fit", c = :blue, lw = 4.0
		)
	end
	plt
end

# ╔═╡ b6e6e3c9-58a5-4e5b-8a83-a516787fb870
begin
	fit_plt = plot(filter_data, 
		#xlims = (0.0, 0.1), ylims = (-rmaxes[1], 0.0),
		c = :black, label_stim = true, grid = false, 
		title = "Amplification"
	)
	
	
	# Plotting the amplification model
	time_cutoff = 0.03 #50ms after stimulus
	for swp in 1:size(data,1), ch in 1:size(data,3)
		amp_model(x, p) = map(t -> AMP(t, p[1], p[2], rmaxes[ch]), x)
		idx_end = findall(data.t .>= time_cutoff)[1]
		xdata = data.t[1:idx_end]
		ydata = data[swp,1:idx_end,ch]
		p0 = [200.0, 0.002]
		lb = [0.0, 0.0]
		ub = [Inf, 0.020]
		fit = curve_fit(amp_model, xdata, ydata, p0, lower = lb, upper = ub)
		if swp == 1 
			label = "Amplification Fit"
		else
			label = ""
		end
		plot!(fit_plt[ch], 
			x -> amp_model(x, fit.param), 
			
			xdata[1], time_cutoff, 
			c = :blue, linewidth = 2.0, label = label)
	end
	fit_plt
end

# ╔═╡ 512f83bb-9b54-4845-92dd-550c4fa4c79f
begin
	IR_plot = plot(title = "Intensity Response Curve")
	IR_model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
	for ch in size(data,3)
		#Plot the intensity response curves
		xdata = IR_analysis[:,:Photons] #set up a 
		ydata = IR_analysis[:,Symbol("Minima_$ch")]
		plot!(IR_plot, xdata, ydata,seriestype = :scatter, xaxis = :log)
		#Fit the intensity response curves
		lb = [0.01, 2.0, 0.0] 
		ub = [Inf, 4.0, Inf]
		pars = [10000.0, 2.0, rmaxes[ch]*-1000]
		fit = curve_fit(IR_model, xdata, ydata, pars, lower = lb, upper = ub)
		plot!(IR_plot,
			x -> IR_model(x, fit.param), 
			LinRange(xdata[1], xdata[end], 10000),
			label = "Fitted IR curve", c = :red, lw = 2.0,
			xlabel = "Intensity (Log Photons)", ylabel = "Response (μV)", 
		)
		vline!(IR_plot, [fit.param[1]], label = "I_1/2", c = :green)
	end
	IR_plot
end

# ╔═╡ Cell order:
# ╠═eb956370-9ba4-11eb-002e-bd530ec32c36
# ╠═e3e3b7e5-901b-41ff-88e2-d296891bdaa2
# ╠═de4ec361-908f-4365-8fb2-ca2c3f08b930
# ╟─70182ef7-d60a-4073-b55b-7565e6c4ab6c
# ╟─f57ebcff-72d2-4c0d-9c39-4caaa967b4ef
# ╟─1e0cfd94-a20b-4691-9b5e-58c5b9e38ff3
# ╟─5b273f51-e6e1-4b2f-9789-6c42c290a8ef
# ╟─b5863ba7-f332-43e1-a73f-20a1f6203201
# ╟─f4eff6dc-1d23-45a3-8514-521048ab5abc
# ╟─0f308136-d3a1-46d3-95a1-a99c6913b524
# ╟─41126fe9-6a3f-495f-9e03-32d980a91bb1
# ╟─b6e6e3c9-58a5-4e5b-8a83-a516787fb870
# ╟─512f83bb-9b54-4845-92dd-550c4fa4c79f
