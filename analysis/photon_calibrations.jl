### A Pluto.jl notebook ###
# v0.12.7

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
target_folder = "D:\\Data\\Calibrations\\525Green_Curve"

# ╔═╡ 03afd290-1ecf-11eb-2738-b5c190ada7a4
stim_files = target_folder |> parse_abf

# ╔═╡ ecea6fb0-21e0-11eb-0442-c7aebca4e7f8
"""
This function is for computing the R-squared
"""
function RSQ(poly::Polynomial, x, y)
	ŷ = poly.(x)
	ȳ = sum(ŷ)/length(ŷ)
	SSE = sum((y-ŷ).^2)
	SST = sum((y.-ȳ).^2)
	1-SSE/SST
end

# ╔═╡ e428a8f2-229a-11eb-032c-e5e98fc1c289
f_T(x) = 10^-x

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
		wavelength = split(file, "\\")[end-1] |> extract_numbers
		println(wavelength)
		info = split(file, "\\")[end] |> filename_extractor
		if info != nothing
			#There might have been an error in ND3 filtering. 
			#For now remove it
			nd, per, ts = info
			println(nd)
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
							f_T(nd|>Float64), 
							per, 
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
     @select {i.LED_Intensity, i.Photons}
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
		ylabel = "Photons at ND0 stim 1ms",
		xlabel = "LED Intensity (%)",
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
T = f_T(0.0)

# ╔═╡ 312e8990-229f-11eb-177a-337c1d75ec08
md" 
We can apply the equation from the first section to the %transmittance function we have just made. 

y = Integrated Photons (Photon per ms)
I = stimulus intensity (%)
D = Optical Density

$ y = (βx^2 + β_2x + α) * 10^{-D} $

However we can transform all of our Optical Density data to equal fractional transferrance, and this can be a much better representation (T)

$ y = T(βx^2 + β_2x + α)  $
"

# ╔═╡ 485f37b0-22ca-11eb-16c6-9f11691d5b3b
Q2 = @from i in data_photon begin
     @where i.Stim_time == 1
	 @where i.LED_Intensity == 1
     @select {i.Transferrance, i.Photons}
	 @collect DataFrame
end

# ╔═╡ 7e2b8a90-22e5-11eb-3918-3b0475614928


# ╔═╡ 5fa0bd40-22ca-11eb-0198-1728507b54e2
begin
	line_poly2 = Polynomials.fit(Q2[!,:Transferrance], Q2[!,:Photons], 1)
	m2 = round(line_poly2.coeffs[2]; digits = 2)
	b2 = round(line_poly2.coeffs[1]; digits = 2)
	r22 = round(RSQ(line_poly2, Q2[!,:Transferrance], Q2[!,:Photons]); digits = 3)
	
	p2 = @df Q2 plot(:Transferrance, :Photons, 
		ylabel = "Photons at 1ms 1%",
		xlabel = "Fractional Transferrance",
		seriestype = :scatter)
	
	plot!(p2, line_poly2, minimum(Q2[!,:Transferrance]), maximum(Q2[!,:Transferrance]), 
		label = "y = $m2 x + $b2 ; R2 = $r22"
	)
end

# ╔═╡ 4be1f280-21dd-11eb-2980-af3191de0ae3
md"#### [3] Stim time and number of photons

We can map the relationship between the stimulus time and the number of photons. It makes sense that this relationship might be linear, because a increase in time would increase photons/ms. 

To demonstrate this experiment we will set the baseline ND filter at 0 and the LED intensity at 100%. We will use the R-squared equation to show how well of a fit linear equations are. 

For the most part our equation can be calculated by a linear equation. The photon density can just be multiplied by a coefficient and the stimulus time
"

# ╔═╡ dbac6710-21dd-11eb-3d20-f9f8e9a298b7
Q3 = @from i in data_photon begin
     @where i.OD == 0
	 @where i.LED_Intensity == 1
     @select {i.Stim_time, i.Photons}
	 @collect DataFrame
end

# ╔═╡ 9a95d8a0-21de-11eb-2304-eb38aee7c819
begin
	line_poly3 = Polynomials.fit(Q3[!,:Stim_time], Q3[!,:Photons], 1)
	m3 = round(line_poly3.coeffs[2]; digits = 2)
	b3 = round(line_poly3.coeffs[1]; digits = 2)
	r23 = round(RSQ(line_poly3, Q3[!,:Stim_time], Q3[!,:Photons]); digits = 3)
	
	p3 = @df Q3 plot(:Stim_time, :Photons, 
		ylabel = "Photons at ND0 1%",
		xlabel = "Stimulus time (ms)",
		seriestype = :scatter)
	
	plot!(p3, line_poly3, minimum(Q3[!,:Stim_time]), maximum(Q3[!,:Stim_time]), 
		label = "y = $m3 x + $b3 ; R2 = $r23"
	)
end

# ╔═╡ 0f35b8f0-1ef2-11eb-3f6d-e70f9a20279d
md"#### Using the dataframe here we can estimate the properties of the non-linear model for prediction of photons
"

# ╔═╡ ee52e530-22e4-11eb-3db4-2131b916852e


# ╔═╡ 51443ed0-22c9-11eb-32d3-11726a77e642
mlm = lm(@formula(Photons~LED_Intensity+Transferrance+Stim_time), data_photon)

# ╔═╡ 7e4e1770-22c9-11eb-0d7f-e33aca5f9183
predict(mlm, DataFrame(Transferrance = 1.0, LED_Intensity = 100.0, Stim_time = 1))

# ╔═╡ c4868470-1eeb-11eb-1b80-57d69a615902
md"
Vx: $(@bind vx Slider(0:90))
Vy: $(@bind vy Slider(0:90))
"

# ╔═╡ 739734a0-1ee7-11eb-074d-a1e6cd45c40b
begin
	Q4 = @from i in data_photon begin
		 @where i.Stim_time == 1
		 @select {i.Transferrance, i.LED_Intensity, i.Photons}
		 @collect DataFrame
	end 
	p = @df Q4 plot(:Transferrance, :LED_Intensity ,:Photons, 
		seriestype = :surface, camera = (vx, vy), 
		colorbar_title = "Photons", 
		xlabel = "Transferrance", ylabel = "Stimulus Time", zlabel = "Photons"
	)	
end

# ╔═╡ Cell order:
# ╠═43750630-1eec-11eb-0a2a-a5c92595dfe0
# ╠═41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
# ╠═7fc2de72-1eea-11eb-009d-d3ed8b3a77bd
# ╠═7fc3efe0-1eea-11eb-25e3-a15506944123
# ╠═19b55080-1ee1-11eb-0b48-7f7133b81e6b
# ╠═7acebd60-1ee7-11eb-0cb9-1f9722a11f29
# ╠═4572cda0-1ece-11eb-33da-212c1c12cb49
# ╟─03afd290-1ecf-11eb-2738-b5c190ada7a4
# ╠═ecea6fb0-21e0-11eb-0442-c7aebca4e7f8
# ╠═e428a8f2-229a-11eb-032c-e5e98fc1c289
# ╠═7cbada1e-1ee0-11eb-39bb-5d77235312ae
# ╟─f4a11470-2299-11eb-1f88-2545159a0c06
# ╠═d83ad020-1ee1-11eb-31fc-e50fcf1030e0
# ╟─982ba89e-229d-11eb-3736-0dac2acd1338
# ╟─a9960350-229e-11eb-18a0-6525509d9176
# ╟─0dbec790-21eb-11eb-1959-412f9c022bca
# ╟─d90c44e0-229a-11eb-2064-251e46b77c5c
# ╠═41351910-22a1-11eb-2ace-f13567348a5f
# ╟─312e8990-229f-11eb-177a-337c1d75ec08
# ╠═485f37b0-22ca-11eb-16c6-9f11691d5b3b
# ╠═7e2b8a90-22e5-11eb-3918-3b0475614928
# ╠═5fa0bd40-22ca-11eb-0198-1728507b54e2
# ╟─4be1f280-21dd-11eb-2980-af3191de0ae3
# ╟─dbac6710-21dd-11eb-3d20-f9f8e9a298b7
# ╟─9a95d8a0-21de-11eb-2304-eb38aee7c819
# ╠═0f35b8f0-1ef2-11eb-3f6d-e70f9a20279d
# ╠═ee52e530-22e4-11eb-3db4-2131b916852e
# ╟─51443ed0-22c9-11eb-32d3-11726a77e642
# ╟─7e4e1770-22c9-11eb-0d7f-e33aca5f9183
# ╟─c4868470-1eeb-11eb-1b80-57d69a615902
# ╠═739734a0-1ee7-11eb-074d-a1e6cd45c40b
