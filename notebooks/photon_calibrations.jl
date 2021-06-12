### A Pluto.jl notebook ###
# v0.14.7

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
using DataFrames, Query, XLSX

# ╔═╡ 7acebd60-1ee7-11eb-0cb9-1f9722a11f29
using Plots, StatsPlots

# ╔═╡ 3e8c4504-328d-461a-8873-afdcf6221e15
begin
	calibration_file = "E:\\Data\\Calibrations\\photon_lookup.xlsx"
	target_folder = "D:\\Calibrations"
	
	if isfile(calibration_file)
		println("File exists")
		rm(calibration_file)
	end
	
	stim_format = (
		"\\", ~, ~, 
		("_", :Wavelength, ~),
		("_", :ND, ~, (".", :Intensity, ~))
	)
	stim_files = target_folder |> parse_abf
	data_photon = DataFrame(
		Wavelength = Int64[], ND = Int64[], Intensity = Float64[], 
		stim_time = Int64[],
		Joules = Float64[], Frequency = Float64[], Photons = Float64[]
	)
	for (idx, file) in enumerate(stim_files)
		format = formatted_split(file, stim_format)
		println(format)
		#Extract the data from the calibration file
		data = extract_abf(file, chs = ["Opt3", "IN 7"])
		#calculate photon energy (joules)
		joules = sum(data, dims = 2)[1,1,1] * data.dt / 0.001 
		frequency = 299792458 / (format.Wavelength*1e-9)
		#the factors are all over the place
		if format.ND == 0 || format.ND == 3
			joules *= 0.1
		elseif format.ND == 1 || format.ND == 2 || format.ND == 4
			joules *= 0.01
		else
		end
		photons = (joules/(6.626e-34*frequency))*1e-8
		#Extrapolate for ND filters and stim time
		for t in 1:4
			push!(data_photon, 
				(
						format.Wavelength, format.ND, format.Intensity+1, t,
						joules*t, frequency, photons*t
					)
				)
			end
		end
	data_photon
	XLSX.writetable(calibration_file, 
		collect(DataFrames.eachcol(data_photon)), DataFrames.names(data_photon)
	)
	data_photon
end

# ╔═╡ 9a6982f6-b953-45f5-8ccb-b61115e9d6da
6.626e-34

# ╔═╡ 75660e71-e90f-444e-bc1a-32087e087118
data_photon |> 
	@filter(_.Wavelength == 525) |>  
	@filter(_.stim_time == 1) |>
	@filter(_.Intensity == 1) |>
	DataFrame

# ╔═╡ f9f446ee-b022-40b9-a2d3-dd0c73eae816
photon_lookup(525, 0, 1, 1, calibration_file)

# ╔═╡ 4c73fb07-2f80-48bf-bc1a-5177704bc9c8
begin
	qGREEN = data_photon |> 
		@filter(_.Wavelength == 525) |>
		#@filter(_.ND == 1) |>
		@filter(_.stim_time == 1) |> 
		DataFrame

	pGREEN = @df qGREEN plot(
		:ND, :Intensity, :Photons, st = :surf, c = colormatch(525), cbar = false
	)
	
	qUV = data_photon |> 
		@filter(_.Wavelength == 365) |>
		#@filter(_.ND == 1) |>
		@filter(_.stim_time == 1) |>  
	DataFrame
	pUV = @df qUV plot(
		:ND, :Intensity, :Photons, st = :surface, c = colormatch(440), cbar = false
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
Ih = $(@bind Ih NumberField(0.0:1e5, default = 400))
n =  $(@bind N NumberField(1:4, default = 3))

number of points = $(@bind n NumberField(1:10, default = 10))

Dimmest = 10^ $(@bind lo NumberField(1:10, default = round(log10(Ih)-1)))

Brightest = 10^ $(@bind hi NumberField(1:10, default = round(log10(Ih)+1)))

"

# ╔═╡ 50598ab8-6fbc-48af-80c5-f22bdcd8520f
ideal_range = exp10.(LinRange(lo, hi, n))

# ╔═╡ e404fd99-d613-42cc-a58d-6b9a869203c7
begin
	plot(i -> IR(i, Ih, N), exp10(lo-1), exp10(hi+1), 
		xlabel = "Log Photons", ylabel = "Response (Normalized)", 
		label = "", c = :black, lw = 4.0, xaxis = :log
	)
	#Pick out several intensity values evenly spaces
	plot!(ideal_range, i -> IR(i, Ih, N), 
		legend = :topleft, c = :red,
		st = :scatter, label = "Ideal Test Flashes")
	plot!([Ih], [IR(Ih, Ih, N)], st = :scatter, marker = :square, label = "I_half")
end	

# ╔═╡ 0065c026-10a6-4659-b4fa-9768da1bc4b5
md"#### Ideal flashes to use:"

# ╔═╡ 9d81b556-460f-49db-8d89-9f33513a1061
begin
	Flashes_Protocol = DataFrame(
		Wavelength = zeros(length(ideal_range)*2),
		ND = zeros(length(ideal_range)*2), 
		Intensity = zeros(length(ideal_range)*2),
		stim_time = zeros(length(ideal_range)*2),
		Target_Photons = [ideal_range... ,ideal_range...],
		Actual_Photons = zeros(length(ideal_range)*2),
		Predicted_Norm_Resp = map(i -> IR(i, Ih, N), [ideal_range... ,ideal_range...])
	)
	for (idx, row) in enumerate(eachrow(Flashes_Protocol))
		println(idx)
		photon_std = 0.1
		photon_min = row.Target_Photons - (row.Target_Photons * photon_std)
		photon_max = row.Target_Photons + (row.Target_Photons * photon_std)
		println("$photon_min < $photon_max")
		wv = idx > length(ideal_range) ? 525 : 365
		qi = data_photon |> 
			@filter(_.Wavelength == wv) |>
			@filter(photon_min < _.Photons < photon_max) |> 
			#@filter(_.Intensity < 50.0) |>
			@map({
				_.Wavelength,
				_.ND, 
				_.Intensity, 
				_.Photons,	
				_.stim_time,
				diff = abs(_.Photons - row.Target_Photons)}) |> 
			@orderby(_.diff) |> 
		DataFrame
		println(qi)
		Flashes_Protocol[idx, :Wavelength] = Int64(wv)
		if !isempty(qi)
			Flashes_Protocol[idx, :ND] = qi.ND[1]|>Int64
			Flashes_Protocol[idx, :Intensity] = qi.Intensity[1]|>Int64
			Flashes_Protocol[idx, :stim_time] = qi.stim_time[1]|>Int64
			Flashes_Protocol[idx, :Actual_Photons] = qi.Photons[1]
		else
			#println(row)
		end		
	end
	Flashes_Protocol
end

# ╔═╡ 7096cc7e-f27f-460e-956c-5c930ae69902
q_search = data_photon |> 
@filter(_.ND == 5 && _.Intensity == 1.0 && _.stim_time == 1.0) |> 
DataFrame

# ╔═╡ df0ca5b9-d7d8-4f51-9a75-2937bb267c7d
begin
	Is = ideal_range
	Rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
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
	plot!(x -> model(x, fit.param), exp10(lo-1), exp10(hi+1),  
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
# ╠═3e8c4504-328d-461a-8873-afdcf6221e15
# ╠═9a6982f6-b953-45f5-8ccb-b61115e9d6da
# ╠═75660e71-e90f-444e-bc1a-32087e087118
# ╠═f9f446ee-b022-40b9-a2d3-dd0c73eae816
# ╠═4c73fb07-2f80-48bf-bc1a-5177704bc9c8
# ╟─1afd7a70-b1fd-49bd-9b8f-b5b81ab63354
# ╠═f73cfdd5-40e5-4d75-995b-44e88d739585
# ╟─87d6f3b9-96c6-428d-a9a1-633d4ca5e837
# ╟─50598ab8-6fbc-48af-80c5-f22bdcd8520f
# ╟─e404fd99-d613-42cc-a58d-6b9a869203c7
# ╟─0065c026-10a6-4659-b4fa-9768da1bc4b5
# ╠═9d81b556-460f-49db-8d89-9f33513a1061
# ╠═7096cc7e-f27f-460e-956c-5c930ae69902
# ╠═df0ca5b9-d7d8-4f51-9a75-2937bb267c7d
# ╟─092b4d0f-105e-4903-940a-1fe6dbfaf5e9
# ╟─4c915ae2-d080-4c3d-8f9b-239282c3a29b
