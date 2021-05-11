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

# ╔═╡ 95b3647c-be5b-4a2b-a45e-ab6646de7b3d
using Revise, PlutoUI

# ╔═╡ 41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
using NeuroPhys

# ╔═╡ 7fc3efe0-1eea-11eb-25e3-a15506944123
using DataFrames, Query

# ╔═╡ 7acebd60-1ee7-11eb-0cb9-1f9722a11f29
using Plots

# ╔═╡ 4572cda0-1ece-11eb-33da-212c1c12cb49
#In order to make this script work, this file is needed
target_folder = "E:\\Data\\Calibrations"

# ╔═╡ 03afd290-1ecf-11eb-2738-b5c190ada7a4
stim_files = target_folder |> parse_abf

# ╔═╡ 3e8c4504-328d-461a-8873-afdcf6221e15
begin
	data_photon = DataFrame(
		Wavelength = Int64[], ND = Int64[], Intensity = Float64[], 
		Joules = Float64[], Frequency = Float64[], Photons = Float64[]
	)
	for (idx, file) in enumerate(stim_files)
		format = formatted_split(file, stim_format)
		#Extract the data from the calibration file
		data = extract_abf(file, chs = ["Opt3", "IN 7"])
		#calculate photon energy (joules)
		joules = sum(data, dims = 2)[1,1,1] * data.dt * 1000
		frequency = 299792458 / (format.Wavelength * 10e-9)
		photons = (joules/(6.626e-34*frequency))*10e-8
		println(frequency)
		push!(data_photon, (format..., joules, frequency, photons))
	end
	data_photon
end

# ╔═╡ 4c73fb07-2f80-48bf-bc1a-5177704bc9c8
begin
	qGREEN = data_photon |> @filter(_.Wavelength == 525) |> DataFrame
	pGREEN = plot(qGREEN.Intensity, qGREEN.Photons, 
		xlabel = "Stimulus intensity", ylabel = "Photons per micron",
		st = :scatter, c = colormatch(525), legend = :topleft
		
	)
	qUV = data_photon |> @filter(_.Wavelength == 365) |> DataFrame
	pUV = plot(qUV.Intensity, qUV.Photons, 
		xlabel = "Stimulus intensity", ylabel = "Photons per micron",
		st = :scatter, c = colormatch(440), legend = :topleft
	)
	plot(pGREEN, pUV, layout = grid(1,2))
end

# ╔═╡ 1afd7a70-b1fd-49bd-9b8f-b5b81ab63354
md"
#### Neutral density filters and Percent Transferrance

The number on the side of our ND filter is the optical density (D). How does optical density relate to %transmittance (T)?

$ T = f_T(D) = 10^{-D}$
- ND0 = 100%   transferrance or multiply photons by 1
- ND1 = 10%    transferrance or multiply photons by 0.1
- ND2 = 1%     transferrance or multiply photons by 0.01
- ND3 = 0.1%   transferrance or multiply photons by 0.001
- ND4 = 0.01%  transferrance or multiply photons by 0.0001
- ND5 = 0.001% transferrance or multiply photons by 0.00001

#### Stim time and number of photons

The number of photons is integrated. This means that we can easily double the amount of photons appearing at 1ms. 
"

# ╔═╡ f73cfdd5-40e5-4d75-995b-44e88d739585
import NeuroPhys: IR #import the IR model

# ╔═╡ 87d6f3b9-96c6-428d-a9a1-633d4ca5e837
md"
#### Finding the perfect stimuli

We can strategically pick stimulus points to try to map the IR curve and get a confidence value

##### Parameters
Ih = $(@bind Ih NumberField(0.0:1e5, default = 11000))
n =  $(@bind N NumberField(1:4, default = 3))

number of points = $(@bind n NumberField(1:10, default = 5))
"

# ╔═╡ e404fd99-d613-42cc-a58d-6b9a869203c7
begin
	lo_phot = Ih
	hi_phot = Ih*4
	plot(i -> IR(i, Ih, N), 1.0, 10e4, 
		xlabel = "Log Photons", ylabel = "Response (Normalized)", 
		label = "", c = :black, lw = 4.0, xaxis = :log
	)
	#Pick out several intensity values evenly spaced
	ideal_range = LinRange(lo_phot, hi_phot, n)
	plot!(ideal_range, i -> IR(i, Ih, N), 
		legend = :topleft, c = :red,
		st = :scatter, label = "Ideal Test Flashes")
	plot!([Ih], [IR(Ih, Ih, N)], st = :scatter, marker = :square, label = "I_half")
end	

# ╔═╡ 0065c026-10a6-4659-b4fa-9768da1bc4b5
md"#### Ideal flashes to use:"

# ╔═╡ 755f4574-1755-43a6-be7d-82b7f1fd68e5
ideal_range

# ╔═╡ fe07a8f0-0899-4274-b99f-89cc93fae213
md"#### Ideal responses:"

# ╔═╡ 8a473d8a-d899-4d34-b7d4-07ab6d3260cf
map(i -> IR(i, Ih, N), ideal_range)

# ╔═╡ df0ca5b9-d7d8-4f51-9a75-2937bb267c7d
begin
	Is = [5540.0, 445177, 890355, 1780710.0, 3561421.0]
	Rs = [3.48, 26.14, 38.42, 41.701, 46.2878]
end

# ╔═╡ 092b4d0f-105e-4903-940a-1fe6dbfaf5e9
begin
	xdata = Is
	ydata = Rs
	lb = [0.0, 1.0, 0.0]
	ub = [Inf, 4.0, Inf]
	model_pars = [Ih, N, 1.0]
	model(x, p) = map(I -> IR(I, p[1], p[2]) * p[3], x)
	fit = curve_fit(model, xdata, ydata, model_pars, lower = lb, upper = ub)
	#Calculate the goodness of fit of the model
	SSE = sum(fit.resid.^2)
	ȳ = sum(model(xdata, fit.param))/length(xdata)
	SST = sum((ydata .- ȳ).^2)
	GOF = 1- SSE/SST
	
	#Plot the results
	plot(xdata, ydata, st = :scatter, xaxis = :log, 
		c = :red
	)
	plot!(x -> model(x, fit.param), 1.0, 10e6,  
		c = :black, lw = 3.0,
		label = "IR R2 = $(round(GOF, digits = 2))", legend = :topleft)
	
end

# ╔═╡ 4c915ae2-d080-4c3d-8f9b-239282c3a29b
md"
#### Fit results
- I_half = $(fit.param[1])
- N = $(fit.param[2])
- Rmax = $(fit.param[3])
R-Squared = $(GOF)
"

# ╔═╡ Cell order:
# ╠═95b3647c-be5b-4a2b-a45e-ab6646de7b3d
# ╠═41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
# ╠═7fc3efe0-1eea-11eb-25e3-a15506944123
# ╠═7acebd60-1ee7-11eb-0cb9-1f9722a11f29
# ╠═4572cda0-1ece-11eb-33da-212c1c12cb49
# ╟─03afd290-1ecf-11eb-2738-b5c190ada7a4
# ╟─3e8c4504-328d-461a-8873-afdcf6221e15
# ╟─4c73fb07-2f80-48bf-bc1a-5177704bc9c8
# ╟─1afd7a70-b1fd-49bd-9b8f-b5b81ab63354
# ╠═f73cfdd5-40e5-4d75-995b-44e88d739585
# ╟─87d6f3b9-96c6-428d-a9a1-633d4ca5e837
# ╟─e404fd99-d613-42cc-a58d-6b9a869203c7
# ╟─0065c026-10a6-4659-b4fa-9768da1bc4b5
# ╟─755f4574-1755-43a6-be7d-82b7f1fd68e5
# ╟─fe07a8f0-0899-4274-b99f-89cc93fae213
# ╟─8a473d8a-d899-4d34-b7d4-07ab6d3260cf
# ╠═df0ca5b9-d7d8-4f51-9a75-2937bb267c7d
# ╟─092b4d0f-105e-4903-940a-1fe6dbfaf5e9
# ╟─4c915ae2-d080-4c3d-8f9b-239282c3a29b
