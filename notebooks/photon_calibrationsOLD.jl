### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 43750630-1eec-11eb-0a2a-a5c92595dfe0
using PlutoUI

# ╔═╡ 41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
using NeuroPhys

# ╔═╡ 7fc2de72-1eea-11eb-009d-d3ed8b3a77bd
using LsqFit, Polynomials, GLM

# ╔═╡ 7fc3efe0-1eea-11eb-25e3-a15506944123
using DataFrames, Query

# ╔═╡ 7acebd60-1ee7-11eb-0cb9-1f9722a11f29
using Plots, StatsPlots

# ╔═╡ 19b55080-1ee1-11eb-0b48-7f7133b81e6b
import NeuroPhys: stim_intensity, photons

# ╔═╡ 4572cda0-1ece-11eb-33da-212c1c12cb49
#In order to make this script work, this file is needed
target_folder = "E:\\Data\\Calibrations\\525Green_Curve"

# ╔═╡ 03afd290-1ecf-11eb-2738-b5c190ada7a4
stim_files = target_folder |> parse_abf

# ╔═╡ 7cbada1e-1ee0-11eb-39bb-5d77235312ae
begin
	data_photon = DataFrame(
		OD = Int64[], 
		Transferrance = Float64[], 
		LED_Intensity = Float64[], 
		Wavelength 	  = Int64[], 
		Stim_time     = Int64[],
		Photon_Energy = Float64[],
		Photons       = Float64[]
	)
	for file in stim_files
		wavelength = split(file, "\\")[end-1] |> number_extractor
		info = split(file, "\\")[end] |> filename_extractor
		if isnothing(info)
			#There might have been an error in ND3 filtering. 
			#For now remove it
			nd, per, ts = info
			if nd == 3.0
				println("Skip this one")
			else
				t, data_array, dt = extract_abf(file; chs = -1);
				stim_t, stim_y = stim_intensity(file; chs = -1);
				photon_d = stim_y .|> photons
				for idx in 1:length(photon_d)
					push!(data_photon, 
						(
							nd, 
							Transferrance(nd|>Float64), 
							per/100, 
							wavelength,
							round(stim_t[idx]),
							(stim_y[idx]/1000)*100,
							photon_d[idx]
						)
					)
				end
			end
		end
	end
end

# ╔═╡ f4a11470-2299-11eb-1f88-2545159a0c06
md"
### Light stimulus equations
In order to calibrate our light stimulus properly we have to think about a few things

We have 4 parts of the photons being emitted. 
1) First an electronic current is passed through a photon emitted LED
2) The photons are run through a neutral density filter
3) The photon emitter is activated for a certain amount of time
4) The light has to travel a certain distance from the detector

We can use our calibration sheets to do much of the work for us
"

# ╔═╡ d83ad020-1ee1-11eb-31fc-e50fcf1030e0
head(data_photon)

# ╔═╡ 982ba89e-229d-11eb-3736-0dac2acd1338
md"
#### [1] Electronic Current and Photons emitted

We can look at the graph for the only simple way to calculate the efficiency of our diode at certain percents. 

Measuring with no Optical density filters at 1ms stimulus time (ND = 0, Stim_time = 1ms)

We can see that we can account for 98% of the variability in the relationship between the intensity percent of the LED and the number of photons being emitted from the LED
"

# ╔═╡ a9960350-229e-11eb-18a0-6525509d9176
Q1 = @from i in data_photon begin
     @where i.OD == 0
	 @where i.Stim_time == 1
     @select {i.Transferrance, i.LED_Intensity, i.Stim_time, i.Photons}
	 @collect DataFrame
end

# ╔═╡ 0dbec790-21eb-11eb-1959-412f9c022bca
begin
	line_poly1 = Polynomials.fit(Q1[!,:LED_Intensity], Q1[!,:Photons], 2)
	m21 = round(line_poly1.coeffs[3]; digits = 2)
	m1 = round(line_poly1.coeffs[2]; digits = 2)
	b1 = round(line_poly1.coeffs[1]; digits = 2)
	r21 = round(RSQ(line_poly1, Q1[!,:LED_Intensity], Q1[!,:Photons]); digits = 3)
	p1 = @df Q1 plot(:LED_Intensity, :Photons, 
		ylabel = "Photons at T1.0 stim 1.0ms",
		xlabel = "LED Intensity",
		seriestype = :scatter)
	
	plot!(p1, line_poly1, minimum(Q1[!,:LED_Intensity]), maximum(Q1[!,:LED_Intensity]), 
		label = "y = $m21 x^2 $m1 x + $b1 R2 = $r21"
	)
end

# ╔═╡ d90c44e0-229a-11eb-2064-251e46b77c5c
md"
#### [2] Neutral density filters and Percent Transferrance

The number on the side of our ND filter is the optical density (D). How does optical density relate to %transmittance (T)?

$ T = f_T(D) = 10^{-D}$
"

# ╔═╡ 41351910-22a1-11eb-2ace-f13567348a5f
T = Transferrance(0.1)

# ╔═╡ 312e8990-229f-11eb-177a-337c1d75ec08
md" 
We can apply the equation from the first section to the %transmittance function we have just made. 

y = Integrated Photons (Photon per ms)
I = stimulus intensity (%)
D = Optical Density

$ y = (βx^2 + β_2x + α) * 10^{-D} $
"

# ╔═╡ 485f37b0-22ca-11eb-16c6-9f11691d5b3b
Q2 = @from i in data_photon begin
     @where i.Stim_time == 1
	 @where i.LED_Intensity == 1
	 #@where i.OD != 2
     @select {i.Transferrance, i.LED_Intensity, i.Stim_time, i.Photons}
end

# ╔═╡ 53a812a0-237d-11eb-05d5-8b6b89581555
md"
Making a equation using both Transferrance (T) and LED Intensity (I)

$ Photons = T(γI^2 + β I + α) $
"

# ╔═╡ c8562f70-23e0-11eb-3077-5bfd6940070a
import NeuroPhys.stimulus_model

# ╔═╡ a4f3baa0-237e-11eb-07e1-b5d7faa6240f
begin
	Qi = @from i in data_photon begin
		 @where i.Stim_time == 1
		 @select {i.Transferrance, i.LED_Intensity, i.Stim_time, i.Photons}
		 @collect DataFrame
	end;
	#The model we want to curve fit
	p0 = [3.3, 45.13, 41.87]
	ivars = Array(Qi[!,[:Transferrance, :LED_Intensity, :Stim_time]]);
	dvars = Qi[!,:Photons];
	excite_fit = curve_fit(stimulus_model, ivars, dvars, p0)
	#Fit the model and predict the parameters
	p_fit = excite_fit.param
	
	Qi[!, :Prediction] = stimulus_model(ivars, p_fit);
	Qi[!, :Percent_Error] = abs.(1.0 .- Qi[!,:Photons]./Qi[!,:Prediction]).*100
	
	Trng = LinRange(0.01, 1.0, 10)
	Irng = LinRange(0.01, 1.0, 10)
	Pvals = Float64[]
	for varT in Trng
		for varI in Irng
			println(varT, varI)
			dVarP = stimulus_model([varT, varI, 1.0], p_fit)
			push!(Pvals, dVarP)
		end
	end
	
	head(Qi)
end

# ╔═╡ c4868470-1eeb-11eb-1b80-57d69a615902
md"
Vx: $(@bind vx Slider(0:90))
Vy: $(@bind vy Slider(0:90))
"

# ╔═╡ 739734a0-1ee7-11eb-074d-a1e6cd45c40b
begin
	p = @df Qi plot(:Transferrance, :LED_Intensity ,:Photons, 
		seriestype = :scatter, camera = (vx, vy), label="Known values",
		xlabel = "Transferrance", ylabel = "LED Intensity", zlabel = "Photons"
	)
	@df Qi plot!(:Transferrance, :LED_Intensity ,:Prediction, 
		seriestype = :scatter, camera = (vx, vy), label = "Predicted Values",
		xlabel = "Transferrance", ylabel = "LED Intensity", zlabel = "Photons"
	)
	plot!(Trng, Irng, Pvals, seriestype = :surface)
end

# ╔═╡ 4be1f280-21dd-11eb-2980-af3191de0ae3
md"#### [3] Stim time and number of photons

We can map the relationship between the stimulus time and the number of photons. It makes sense that this relationship might be linear, because a increase in time would increase photons/ms. 

To demonstrate this experiment we will set the baseline ND filter at 0 (T1.0) and the LED intensity at 100% (I1.0).

For the most part our equation can be calculated by a linear term at the end of our equation. 

The final equation looks like this

$ y = T(γI^2 + βI + α)t_{stim} $

The fit might need small adjustments in order to account for any variances
"

# ╔═╡ 090b1050-23a8-11eb-2c0f-69c9cffb9fb7
begin
	Q3 = @from i in data_photon begin
		@where i.OD == 0
		@where i.LED_Intensity == 1.0
		@select {i.Transferrance, i.LED_Intensity, i.Stim_time, i.Photons}
		@collect DataFrame
	end
	ivars_f = Array(data_photon[!,[:Transferrance, :LED_Intensity, :Stim_time]]);
	dvars_f = data_photon[!,:Photons];
	excite_fit_f = curve_fit(stimulus_model, ivars_f, dvars_f, p0)
	p_fit_f = excite_fit_f.param
	
	ivars3 = Array(Q3[!,[:Transferrance, :LED_Intensity, :Stim_time]]);
	Q3[!, :Prediction] = stimulus_model(ivars3, p_fit)
	Q3[!, :Percent_Error] = abs.(1.0 .- (Q3[!,:Photons]./Q3[!,:Prediction])) * 100
	Q3
end

# ╔═╡ 9a95d8a0-21de-11eb-2304-eb38aee7c819
begin
	p3 = @df Q3 plot(:Stim_time, :Photons, label = "Actual Measurements",
		ylabel = "Photons at T1.0 I1.0",
		xlabel = "Stimulus time (ms)",
		seriestype = :scatter, m = :star, markersize = 10.0
	)
	
	@df Q3 plot!(p3, :Stim_time, :Photons, label = "Model predicition",
		ylabel = "Photons at T1.0 I1.0",
		xlabel = "Stimulus time (ms)",
		seriestype = :scatter
	)
	stim_data = Q3[!, :Stim_time]
	γ = round(p_fit_f[1], digits = 2)
	β = round(p_fit_f[2], digits = 2)
	α = round(p_fit_f[3], digits = 2)
	
	plot!(x -> stimulus_model([1.0, 1.0, x], p_fit), 
		minimum(stim_data), maximum(stim_data), 
		label = "P = T($γ I^2 + $β I + $α)t_stim"
	)
end

# ╔═╡ 6ad36cc0-23ac-11eb-1144-691d7fb88a47
md"
The final fit equation is:

$ Photons = t_{stim} T (γI^2 + βI + α) $

$ Photons = t_{stim} T ($γ I^2 + $β I + $α) $
"

# ╔═╡ Cell order:
# ╠═43750630-1eec-11eb-0a2a-a5c92595dfe0
# ╠═41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
# ╠═7fc2de72-1eea-11eb-009d-d3ed8b3a77bd
# ╠═7fc3efe0-1eea-11eb-25e3-a15506944123
# ╠═19b55080-1ee1-11eb-0b48-7f7133b81e6b
# ╠═7acebd60-1ee7-11eb-0cb9-1f9722a11f29
# ╠═4572cda0-1ece-11eb-33da-212c1c12cb49
# ╟─03afd290-1ecf-11eb-2738-b5c190ada7a4
# ╠═7cbada1e-1ee0-11eb-39bb-5d77235312ae
# ╟─f4a11470-2299-11eb-1f88-2545159a0c06
# ╟─d83ad020-1ee1-11eb-31fc-e50fcf1030e0
# ╟─982ba89e-229d-11eb-3736-0dac2acd1338
# ╟─a9960350-229e-11eb-18a0-6525509d9176
# ╠═0dbec790-21eb-11eb-1959-412f9c022bca
# ╟─d90c44e0-229a-11eb-2064-251e46b77c5c
# ╠═41351910-22a1-11eb-2ace-f13567348a5f
# ╟─312e8990-229f-11eb-177a-337c1d75ec08
# ╟─485f37b0-22ca-11eb-16c6-9f11691d5b3b
# ╟─53a812a0-237d-11eb-05d5-8b6b89581555
# ╟─a4f3baa0-237e-11eb-07e1-b5d7faa6240f
# ╠═c8562f70-23e0-11eb-3077-5bfd6940070a
# ╟─c4868470-1eeb-11eb-1b80-57d69a615902
# ╟─739734a0-1ee7-11eb-074d-a1e6cd45c40b
# ╟─4be1f280-21dd-11eb-2980-af3191de0ae3
# ╟─090b1050-23a8-11eb-2c0f-69c9cffb9fb7
# ╟─9a95d8a0-21de-11eb-2304-eb38aee7c819
# ╟─6ad36cc0-23ac-11eb-1144-691d7fb88a47
