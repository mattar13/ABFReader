### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 95b3647c-be5b-4a2b-a45e-ab6646de7b3d
using Revise

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


# ╔═╡ Cell order:
# ╠═95b3647c-be5b-4a2b-a45e-ab6646de7b3d
# ╠═41e9ea6e-1ecd-11eb-3dfc-ab8db76a04bc
# ╠═7fc3efe0-1eea-11eb-25e3-a15506944123
# ╠═7acebd60-1ee7-11eb-0cb9-1f9722a11f29
# ╠═4572cda0-1ece-11eb-33da-212c1c12cb49
# ╟─03afd290-1ecf-11eb-2738-b5c190ada7a4
# ╠═3e8c4504-328d-461a-8873-afdcf6221e15
# ╟─4c73fb07-2f80-48bf-bc1a-5177704bc9c8
# ╠═1afd7a70-b1fd-49bd-9b8f-b5b81ab63354
# ╠═f73cfdd5-40e5-4d75-995b-44e88d739585
